import sys
import time
import math
import argparse
from pathlib import Path
from unittest.mock import patch
import pandas as pd

PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import metabo_engineer as me
from cobra.io import read_sbml_model

def parse_args():
    parser = argparse.ArgumentParser(description="Pipeline Universal de Diseño Asistido")
    parser.add_argument("--modelo", type=str, required=True, help="Ruta al archivo SBML")
    parser.add_argument("--target", type=str, required=True, help="Reacción objetivo")
    parser.add_argument("--biomasa", type=str, default=None, help="ID de la biomasa")
    parser.add_argument("--ambiente", type=str, default="anaerobio", choices=["aerobio", "anaerobio"])
    parser.add_argument("--seed", type=int, default=42, help="Semilla para el GA")
    parser.add_argument("--k", type=int, default=3, help="Número de genes a apagar")
    parser.add_argument("--gens", type=int, default=10, help="Generaciones del GA")
    parser.add_argument("--pop", type=int, default=20, help="Población del GA")
    parser.add_argument("--seeding", type=str, default=None, help="Genes a inyectar (ej. b3738)")
    return parser.parse_args()

class _FakePrompt:
    responses = []
    @classmethod
    def ask(cls, *a, **k): return cls.responses.pop(0)

class _FakeIntPrompt:
    responses = []
    @classmethod
    def ask(cls, *a, **k): return int(cls.responses.pop(0))

def run_ga(model, target, seed, gens, pop, k):
    captured = {}
    def _fake_export(modelo_id, quimico, metodo, genes, fitness, seed_, parametros_extra=None):
        captured.update({"genes": list(genes), "fitness": float(fitness)})
    
    _FakePrompt.responses = [target]
    _FakeIntPrompt.responses = [gens, pop, k]
    
    patches = [
        patch.object(me, 'Prompt', _FakePrompt),
        patch.object(me, 'IntPrompt', _FakeIntPrompt),
        patch.object(me, '_exportar_resultado_json', _fake_export),
        patch.object(me, 'analisis_profundo', lambda *a, **k: None),
        patch.object(me, 'generar_documento_markdown', lambda *a, **k: None),
        patch.object(me, 'get_organism_id', lambda m: "83333"),
        patch.object(me, 'traducir_gen_uniprot', lambda g, o: (g, "patched")),
    ]
    for p in patches: p.start()
    try:
        me.disenar_cepa_ia(model, "Custom", "Custom", seed=seed)
    finally:
        for p in patches: p.stop()
    return captured

# ==========================================
# FUNCIÓN DE SEGURIDAD (Evita errores con listas vacías o NaN)
# ==========================================
def safe_evaluate(m, genes, target, biomass):
    if not genes:
        with m:
            m.objective = biomass
            wt_mu = m.slim_optimize()
            
            # Blindaje contra NaN o None
            if wt_mu is None or math.isnan(wt_mu) or math.isinf(wt_mu) or wt_mu < 0.05: 
                return 0.0
                
            bio_rxn = m.reactions.get_by_id(biomass)
            bio_rxn.lower_bound = wt_mu * 0.9
            m.objective = target
            m.objective.direction = 'min'
            sol = m.optimize()
            return max(0.0, sol.objective_value) if sol.status == 'optimal' else 0.0
    return me.evaluar_bacteria(m, genes, target, biomass)

def main():
    args = parse_args()
    
    print("="*70)
    print(f" 🧬 PIPELINE DE DISEÑO ASISTIDO - {args.modelo}")
    print("="*70)

    print("\n[ETAPA 1] Construyendo entorno experimental...")
    t0 = time.time()
    
    try:
        m = read_sbml_model(args.modelo)
    except Exception as e:
        print(f"  ❌ Error cargando modelo: {e}")
        return
        
    m.solver = 'glpk'
    
    if args.biomasa is None:
        args.biomasa = next((r.id for r in m.reactions if "biomass" in r.id.lower()), None)
        if not args.biomasa:
            print("  ❌ No se encontró reacción de biomasa.")
            return
            
    m = me.configurar_ambiente(m, args.ambiente)
    m.solver = 'glpk'
    
    wt_mu = m.slim_optimize()
    print(f"  ✓ Modelo: {len(m.reactions)} rxns | Biomasa: {args.biomasa}")
    
    # BLINDAJE: Detectar NaN, infinito o crecimiento nulo
    if wt_mu is None or math.isnan(wt_mu) or math.isinf(wt_mu) or wt_mu < 0.01:
        print(f"  ❌ REGLA DE PARADA: Crecimiento WT inválido o nulo ({wt_mu}). Abortando.")
        return
        
    print(f"  ✓ Crecimiento WT en {args.ambiente}: {wt_mu:.4f} h⁻¹")
    print("  ✅ CHECKPOINT PASADO.\n")

    # ==========================================
    # ETAPA 2: EXPLORACIÓN
    # ==========================================
    if args.seeding:
        print(f"[ETAPA 2] Usando Seeding Manual (Genes: {args.seeding})")
        ga_genes = [g.strip() for g in args.seeding.split(",") if g.strip()]
        ga_fitness = safe_evaluate(m.copy(), ga_genes, args.target, args.biomasa)
    else:
        print(f"[ETAPA 2] Explorando con el GA (Seed={args.seed}, k={args.k})...")
        res = run_ga(m, args.target, args.seed, args.gens, args.pop, args.k)
        ga_genes = res.get('genes', [])
        ga_fitness = res.get('fitness', -1e6)

    print(f"  ✓ Genes propuestos: {ga_genes}")
    print(f"  ✓ Fitness inicial: {ga_fitness:.4f}\n")

    if ga_fitness <= -1e5:
        print("❌ VEREDICTO: El diseño inicial es LETAL (mata a la bacteria).")
        print("   Sugerencia: Cambia la semilla o usa otro conjunto de genes en --seeding.")
        return

    # ==========================================
    # ETAPA 3: REFINAMIENTO (LEAVE-ONE-OUT)
    # ==========================================
    print("[ETAPA 3] Refinando diseño (Leave-One-Out)...")
    core_genes = list(ga_genes)
    wt_fitness = safe_evaluate(m.copy(), [], args.target, args.biomasa)
    
    for gen in list(ga_genes):
        if not core_genes: break
        test_genes = [g for g in core_genes if g != gen]
        f_minus = safe_evaluate(m.copy(), test_genes, args.target, args.biomasa)
        
        delta = f_minus - ga_fitness
        if delta > -0.001: 
            core_genes.remove(gen)
            print(f"  ✂️  Pasajero eliminado: {gen}")
        else:
            print(f"  📌 Gen importante:     {gen}")

    print(f"  ✓ Diseño minimalista: {core_genes}\n")

    # ==========================================
    # ETAPA 4: VERIFICACIÓN
    # ==========================================
    print("[ETAPA 4] Verificación biológica...")
    
    if not core_genes:
        print("  ⚠️ El GA no encontró ningún gen útil. El diseño es igual al WT.")
        final_fitness = wt_fitness
        mut_mu = wt_mu
    else:
        final_fitness = safe_evaluate(m.copy(), core_genes, args.target, args.biomasa)
        m_mut = m.copy()
        for g in core_genes:
            if g in m_mut.genes: m_mut.genes.get_by_id(g).knock_out()
        
        mut_mu = m_mut.slim_optimize()
        print(f"  ✓ Fitness final: {final_fitness:.4f}")
        
        if mut_mu is None or math.isnan(mut_mu):
            print("  ⚠️ Crecimiento mutante: inválido (nan)")
        else:
            print(f"  ✓ Crecimiento mutante: {mut_mu:.4f} h⁻¹")

    # ==========================================
    # ETAPA 5: DOCUMENTACIÓN
    # ==========================================
    print("\n[ETAPA 5] Reporte Final")
    print("-" * 50)
    
    if not core_genes:
        print("❌ VEREDICTO: 'El GA se quedó corto'. Usar OptKnock exacto o cambiar semilla.")
        return

    nombres = {}
    for g in core_genes:
        try:
            nombre_real = m.genes.get_by_id(g).name
            nombres[g] = nombre_real if nombre_real else g
        except KeyError:
            nombres[g] = g

    df = pd.DataFrame({
        "Gen KO": [nombres[g] for g in core_genes],
        "ID": core_genes,
    })
    
    print(df.to_string(index=False))
    print("-" * 50)
    
    if final_fitness > 1.0:
        print("🏆 VEREDICTO: 'Diseño listo para el lab'. Acoplamiento validado.")
    elif final_fitness > wt_fitness:
        print("⚠️ VEREDICTO: 'Interesante pero frágil'. Requiere validación.")
    else:
        print("❌ VEREDICTO: 'El GA se quedó corto'.")
    
    print(f"\n[Tiempo total: {time.time()-t0:.1f} segundos]")

if __name__ == "__main__":
    main()