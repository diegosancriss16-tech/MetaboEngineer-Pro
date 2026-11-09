"""
Fase 2 (Optimizada): RobustKnock real con compresión de red.
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

logging.getLogger().setLevel(logging.ERROR)  # Silenciar warnings de GPR

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--model", required=True)
    p.add_argument("--target", default="EX_succ_e")
    p.add_argument("--biomass", default=None)
    p.add_argument("--k", type=int, default=2)
    p.add_argument("--anaerobic", action="store_true")
    p.add_argument("--time-limit", type=int, default=1800)
    p.add_argument("--max-solutions", type=int, default=1)
    p.add_argument("--solver", default="glpk")
    p.add_argument("--min-growth", type=float, default=0.05)
    args = p.parse_args()

    print(f"Cargando {args.model}...")
    model = read_sbml_model(args.model)
    model.solver = 'glpk' 
    
    ambiente = "anaerobio" if args.anaerobic else "aerobio"
    model = me.configurar_ambiente(model, ambiente)
    model.solver = 'glpk'

    if args.biomass is None:
        args.biomass = next(r.id for r in model.reactions if "biomass" in r.id.lower())
    if args.target in model.reactions:
        model.reactions.get_by_id(args.target).bounds = (0, 1000.0)

    print(f"  Ambiente:   {ambiente}")
    print(f"  Biomass:    {args.biomass}")
    print(f"  Target:     {args.target}")
    print(f"  k:          {args.k}")
    print(f"  Solver:     {args.solver}")

    print("\n--- Baseline WT ---")
    wt = wt_baseline(model, args.target, args.biomass)
    print(f"  WT max growth: {model.slim_optimize():.4f} 1/h")
    print(f"  WT max@90%:   {wt['at_90pct_growth']['max']:.4f}")
    print(f"  WT min@100%:  {wt['at_100pct_growth']['min']:.4f}")

    print(f"\n=== RobustKnock (k={args.k}, compress=True) ===")
    
    # Diccionario para mapear nombres a IDs de genes (corrección del error get_by_name)
    gene_name_to_id = {g.name: g.id for g in model.genes}

    # 1. SDModule SIN el parámetro compress
    module = sd.SDModule(
        model,
        'robustknock',
        inner_objective=args.biomass,
        outer_objective=args.target,
        constraints=f"{args.biomass} >= {args.min_growth}"
    )

    with Timer() as t:
        try:
            # 2. compress=True va AQUÍ, en compute_strain_designs
            sols = sd.compute_strain_designs(
                model,
                sd_modules=[module],
                max_solutions=args.max_solutions,
                max_cost=args.k,
                solution_approach='best',
                solver=args.solver,
                time_limit=args.time_limit,
                gene_kos=True,
                compress=True  # <--- LA MAGIA ESTÁ AQUÍ
            )
        except Exception as e:
            print(f"❌ Error straindesign: {e}")
            return

    print(f"  Tiempo: {t.elapsed:.1f}s")
    print(f"  Soluciones crudas: {len(sols.gene_sd)}")

    rows = []
    for i, gene_dict in enumerate(sols.gene_sd):
        raw_genes = [g for g, v in gene_dict.items() if v == -1]
        
        genes = []
        for g in raw_genes:
            if g in model.genes:
                genes.append(g)
            elif g in gene_name_to_id:
                genes.append(gene_name_to_id[g])  # Conversión segura con el diccionario
            else:
                print(f"  [!] Gen no encontrado: {g}")
                
        if not genes:
            continue
            
        try:
            m90 = production_metrics(model, genes, args.target, args.biomass, 0.9)
            m100 = production_metrics(model, genes, args.target, args.biomass, 1.0)
            rows.append({
                "solution_id": i,
                "genes": genes,
                "prod_max_at_90pct": m90["max"],
                "prod_min_at_100pct": m100["min"],
                "ratio_at_100pct": m100["ratio"],
                "mutant_growth_90pct": m90["growth"],
            })
        except Exception as e:
            print(f"  [!] Solución {i} falló: {e}")

    rows_sorted = sorted(rows, key=lambda r: r["prod_min_at_100pct"], reverse=True)
    
    print("\n--- Top soluciones por acoplamiento (min@100%) ---")
    for r in rows_sorted[:3]:
        print(f"  Genes: {r['genes']}")
        print(f"    min@100%={r['prod_min_at_100pct']:.4f}, ratio={r['ratio_at_100pct']:.4f}, growth={r['mutant_growth_90pct']:.4f}")

    best = rows_sorted[0] if rows_sorted else None

    results = {
        "method": "RobustKnock (compressed)",
        "solver": args.solver,
        "k": args.k,
        "ambiente": ambiente,
        "wt_baseline": wt,
        "n_solutions_valid": len(rows),
        "solutions_top10": rows_sorted[:10],
        "best": best,
        "time_s": t.elapsed,
        "params": vars(args),
    }
    save_results("straindesign_robust", results)


if __name__ == "__main__":
    main()