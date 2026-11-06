"""
Fase 2 (Universal): OptKnock real con straindesign.
Funciona en cualquier modelo SBML (E. coli, B. subtilis, S. cerevisiae, etc.)
usando solvers open-source (GLPK/SCIP).
"""
import argparse
import sys
import logging
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import straindesign as sd
from cobra.io import read_sbml_model
import metabo_engineer as me
from common import production_metrics, wt_baseline, save_results, Timer

logging.getLogger().setLevel(logging.ERROR)

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--model", required=True)
    p.add_argument("--target", default="EX_succ_e")
    p.add_argument("--biomass", default=None)
    p.add_argument("--k", type=int, default=2)
    p.add_argument("--anaerobic", action="store_true")
    p.add_argument("--time-limit", type=int, default=600)
    p.add_argument("--max-solutions", type=int, default=5)
    p.add_argument("--solver", default="glpk")
    p.add_argument("--default-medium", action="store_true", help="Ignorar configurar_ambiente")
    args = p.parse_args()

    print(f"Cargando {args.model}...")
    model = read_sbml_model(args.model)
    model.solver = 'glpk' 
    
    # 1. Configuración del medio
    if not args.default_medium:
        ambiente = "anaerobio" if args.anaerobic else "aerobio"
        model = me.configurar_ambiente(model, ambiente)
        model.solver = 'glpk'
    else:
        print("  Usando medio por defecto del modelo (--default-medium).")

    # 2. FIX UNIVERSAL: Resolver colisiones entre nombres de genes e IDs de reacciones
    rxn_ids = set(model.reactions.list_attr("id"))
    renamed = 0
    for gene in model.genes:
        if gene.name in rxn_ids:
            gene.name = f"{gene.name}_gene"
            renamed += 1
    if renamed:
        print(f"  Renombrados {renamed} genes que colisionaban con reacciones.")

    # 3. Detectar biomasa y target
    if args.biomass is None:
        args.biomass = next((r.id for r in model.reactions if "biomass" in r.id.lower()), None)
        if not args.biomass:
            print("❌ No se encontró reacción de biomasa.")
            return

    if args.target in model.reactions:
        model.reactions.get_by_id(args.target).bounds = (0, 1000.0)
    else:
        print(f"❌ Target {args.target} no encontrado en el modelo.")
        return

    print(f"  Biomass:    {args.biomass}")
    print(f"  Target:     {args.target}")
    print(f"  k:          {args.k}")
    print(f"  Solver:     {args.solver}")

    print("\n--- Baseline WT ---")
    wt_growth = model.slim_optimize()
    if wt_growth is None or wt_growth < 0.01:
        print(f"  ⚠️ WT max growth: {wt_growth} 1/h")
        min_growth = 0.0
    else:
        print(f"  WT max growth: {wt_growth:.4f} 1/h")
        wt = wt_baseline(model, args.target, args.biomass)
        print(f"  WT max@90%:   {wt['at_90pct_growth']['max']:.4f}")
        print(f"  WT min@100%:  {wt['at_100pct_growth']['min']:.4f}")
        min_growth = wt_growth * 0.1

    print(f"\n=== OptKnock real (k={args.k}, compress=True) ===")
    
    # 4. Crear el módulo con Plan B (si la restricción falla, se quita)
    if min_growth > 0:
        print(f"  Restricción de crecimiento mínimo: >= {min_growth:.4f} 1/h")
    
    try:
        module = sd.SDModule(
            model,
            'optknock',
            inner_objective=args.biomass,
            outer_objective=args.target,
            constraints=f"{args.biomass} >= {min_growth}" if min_growth > 0 else None
        )
    except Exception as e:
        print(f"⚠️ Restricción de crecimiento inviable. Reintentando sin restricción...")
        module = sd.SDModule(
            model,
            'optknock',
            inner_objective=args.biomass,
            outer_objective=args.target,
            constraints=None
        )

    gene_name_to_id = {g.name: g.id for g in model.genes}

    # 5. Resolver
    with Timer() as t:
        try:
            sols = sd.compute_strain_designs(
                model,
                sd_modules=[module],
                max_solutions=args.max_solutions,
                max_cost=args.k,
                solution_approach='best',
                solver=args.solver,
                time_limit=args.time_limit,
                gene_kos=True,
                compress=True
            )
        except Exception as e:
            print(f"❌ Error en el solver: {e}")
            return

    print(f"  Tiempo: {t.elapsed:.1f}s")
    print(f"  Soluciones crudas: {len(sols.gene_sd)}")

    # 6. Evaluar y filtrar soluciones letales
    rows = []
    for i, gene_dict in enumerate(sols.gene_sd):
        raw_genes = [g for g, v in gene_dict.items() if v == -1]
        
        genes = []
        for g in raw_genes:
            if g in model.genes:
                genes.append(g)
            elif g in gene_name_to_id:
                genes.append(gene_name_to_id[g])
                
        if not genes:
            continue
            
        try:
            m90 = production_metrics(model, genes, args.target, args.biomass, 0.9)
            m100 = production_metrics(model, genes, args.target, args.biomass, 1.0)
            
            if m100["growth"] < 0.01:
                continue
                
            rows.append({
                "solution_id": i,
                "genes": genes,
                "prod_max_at_90pct": m90["max"],
                "prod_min_at_100pct": m100["min"],
                "ratio_at_100pct": m100["ratio"],
                "mutant_growth_90pct": m90["growth"],
            })
        except Exception as e:
            pass

    rows_sorted = sorted(rows, key=lambda r: r["prod_max_at_90pct"], reverse=True)
    
    print("\n--- Top soluciones viables por max@90% ---")
    if not rows_sorted:
        print("  [!] No se encontraron soluciones viables (todas eran letales).")
    else:
        for r in rows_sorted[:3]:
            print(f"  Genes: {r['genes']}")
            print(f"    max@90%={r['prod_max_at_90pct']:.4f}, min@100%={r['prod_min_at_100pct']:.4f}, growth={r['mutant_growth_90pct']:.4f}")

    best = rows_sorted[0] if rows_sorted else None

    results = {
        "method": "OptKnock (compressed)",
        "solver": args.solver,
        "k": args.k,
        "wt_growth": wt_growth,
        "n_solutions_valid": len(rows),
        "solutions_top10": rows_sorted[:10],
        "best": best,
        "time_s": t.elapsed,
        "params": vars(args),
    }
    save_results("straindesign", results)

if __name__ == "__main__":
    main()