import sys
import time
from pathlib import Path
from unittest.mock import patch
import pandas as pd

PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import metabo_engineer as me
from cobra.io import read_sbml_model
from cobra.flux_analysis import flux_variability_analysis

# ==========================================
# CONFIGURACIÓN DEL PIPELINE
# ==========================================
MODELO = "iML1515.xml"
BIOMASA = "BIOMASS_Ec_iML1515_core_75p37M"
TARGET = "EX_succ_e"
AMBIENTE = "anaerobio"
SEED = 42

class _FakePrompt:
    responses = []
    @classmethod
    def ask(cls, *a, **k): return cls.responses.pop(0)

class _FakeIntPrompt:
    responses = []
    @classmethod
    def ask(cls, *a, **k): return int(cls.responses.pop(0))

def run_ga(model, target, seed):
    captured = {}
    def _fake_export(modelo_id, quimico, metodo, genes, fitness, seed_, parametros_extra=None):
        captured.update({"genes": list(genes), "fitness": float(fitness)})
    
    _FakePrompt.responses = [target]
    _FakeIntPrompt.responses = [10, 20, 3] # gens=10, pop=20, mut=3 (rápido para demo)
    
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
        me.disenar_cepa_ia(model, "E. coli", "iML1515", seed=seed)
    finally:
        for p in patches: p.stop()
    return captured

def main():
    print("="*70)
    print(" 🧬 PIPELINE DE DISEÑO ASISTIDO - MetaboEngineer Pro")
    print("="*70)

    # ==========================================
    # ETAPA 1: CONSTRUIR Y VALIDAR ENTORNO
    # ==========================================
    print("\n[ETAPA 1] Construyendo entorno experimental...")
    t0 = time.time()
    m = read_sbml_model(MODELO)
    m.solver = 'glpk'
    m = me.configurar_ambiente(m, AMBIENTE)
    m.solver = 'glpk'
    
    wt_mu = m.slim_optimize()
    print(f"  ✓ Modelo cargado: {len(m.reactions)} reacciones.")
    print(f"  ✓ Crecimiento WT en {AMBIENTE}: {wt_mu:.4f} h⁻¹")
    
    if not (0.15 < wt_mu < 0.17):
        print("  ❌ REGLA DE PARADA: El crecimiento WT no coincide con la literatura (~0.157). Abortando.")
        return
    print("  ✅ CHECKPOINT PASADO. El entorno es válido.\n")

    # ==========================================
    # ETAPA 2: EXPLORACIÓN (GA)
    # ==========================================
    print(f"[ETAPA 2] Explorando con el GA (Seed={SEED}, k=3)... (esto tarda ~3 min)")
    res = run_ga(m, TARGET, SEED)
    ga_genes = res.get('genes', [])
    ga_fitness = res.get('fitness', -1e6)
    print(f"  ✓ GA propuso {len(ga_genes)} genes: {ga_genes}")
    print(f"  ✓ Fitness inicial: {ga_fitness:.4f}\n")

    # ==========================================
    # ETAPA 3: REFINAMIENTO (LEAVE-ONE-OUT)
    # ==========================================
    print("[ETAPA 3] Refinando diseño (Leave-One-Out)...")
    core_genes = list(ga_genes)
    
    for gen in list(ga_genes):
        test_genes = [g for g in core_genes if g != gen]
        # Evaluamos sin el gen (si la lista queda vacía, el fitness es WT)
        f_minus = me.evaluar_bacteria(m.copy(), test_genes, TARGET, BIOMASA) if test_genes else me.evaluar_bacteria(m.copy(), [], TARGET, BIOMASA)
        
        delta = f_minus - ga_fitness
        if delta > -0.001: # Si al quitarlo el fitness NO baja, es un pasajero
            core_genes.remove(gen)
            print(f"  ✂️  Pasajero eliminado: {gen} (Fitness sin él: {f_minus:.4f})")
        else:
            print(f"  📌 Gen importante:     {gen} (Fitness sin él cae a {f_minus:.4f})")

    print(f"  ✓ Diseño minimalista: {core_genes}\n")

    # ==========================================
    # ETAPA 4: VERIFICACIÓN (CONTROLES)
    # ==========================================
    print("[ETAPA 4] Verificación biológica del diseño mínimo...")
    final_fitness = me.evaluar_bacteria(m.copy(), core_genes, TARGET, BIOMASA)
    
    print(f"  ✓ Fitness del diseño final: {final_fitness:.4f}")
    
    print("  Verificando viabilidad del mutante...")
    m_mut = m.copy()
    for g in core_genes:
        if g in m_mut.genes: m_mut.genes.get_by_id(g).knock_out()
    mut_mu = m_mut.slim_optimize()
    
    if mut_mu < 0.05:
        print(f"  ⚠️ ADVERTENCIA: Cepa moribunda (μ={mut_mu:.4f})")
    else:
        print(f"  ✅ Cepa viable (μ={mut_mu:.4f} h⁻¹)\n")

    # ==========================================
    # ETAPA 5: DOCUMENTACIÓN
    # ==========================================
    print("[ETAPA 5] Reporte Final")
    print("-" * 50)
    
    # Traducir nombres para el reporte
    nombres = {}
    for g in core_genes:
        nombre_real = m.genes.get_by_id(g).name
        nombres[g] = nombre_real if nombre_real else g

    df = pd.DataFrame({
        "Gen KO": [nombres[g] for g in core_genes],
        "ID BiGG": core_genes,
    })
    
    print(df.to_string(index=False))
    print("-" * 50)
    
    if final_fitness > 1.0:
        print("🏆 VEREDICTO: 'Diseño listo para el lab'. Acoplamiento validado.")
    elif final_fitness > wt_mu:
        print("⚠️ VEREDICTO: 'Interesante pero frágil'. Requiere validación termodinámica.")
    else:
        print("❌ VEREDICTO: 'El GA se quedó corto'. Usar OptKnock exacto o Seeding.")
    
    print(f"\n[Tiempo total de cómputo: {time.time()-t0:.1f} segundos]")

if __name__ == "__main__":
    main()
