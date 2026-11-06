import sys
import numpy as np
from pathlib import Path
from unittest.mock import patch

PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import metabo_engineer as me
from cobra.io import read_sbml_model

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
    def _fake_export(modelo_id, quimico, metodo, genes, fitness, seed, parametros_extra=None):
        captured.update({"genes": list(genes), "fitness": float(fitness)})
    
    patches = [
        patch.object(me, 'Prompt', _FakePrompt),
        patch.object(me, 'IntPrompt', _FakeIntPrompt),
        patch.object(me, '_exportar_resultado_json', _fake_export),
        patch.object(me, 'analisis_profundo', lambda *a, **k: None),
        patch.object(me, 'generar_documento_markdown', lambda *a, **k: None),
        patch.object(me, 'get_organism_id', lambda m: "83333"),
        patch.object(me, 'traducir_gen_uniprot', lambda g, o: (g, "patched")),
    ]
    
    _FakePrompt.responses = [target]
    _FakeIntPrompt.responses = [10, 20, 3] # gens, pop, mut
    
    for p in patches: p.start()
    try:
        me.disenar_cepa_ia(model, "E. coli", "iML1515", seed=seed)
    finally:
        for p in patches: p.stop()
    return captured

def main():
    print("Cargando modelo...")
    model = read_sbml_model("iML1515.xml")
    model.solver = 'glpk'  # <--- FIX: forzar GLPK
    
    model = me.configurar_ambiente(model, "anaerobio")
    model.solver = 'glpk'  # <--- FIX: configurar_ambiente devuelve copia que hereda Gurobi
    
    semillas = [42, 100, 2024]
    resultados = []
    
    for s in semillas:
        print(f"\n--- Corriendo GA con semilla {s} (FVA-min) ---")
        try:
            res = run_ga(model, "EX_succ_e", seed=s)
            print(f"  Genes: {res.get('genes')}")
            print(f"  Fitness (FVA-min): {res.get('fitness')}")
            resultados.append(res.get('fitness', -1e6))
        except Exception as e:
            print(f"  ❌ Error en semilla {s}: {e}")
            resultados.append(-1e6)
        
    print("\n=== RESULTADOS BEST-OF-N ===")
    validos = [r for r in resultados if r > -1e5]
    if validos:
        print(f"  Fitness media: {np.mean(validos):.4f}")
        print(f"  Desviación estándar: {np.std(validos):.4f}")
        print(f"  Mejor fitness: {max(validos):.4f}")
        print(f"  Fitness individuales: {[f'{v:.4f}' for v in validos]}")
    else:
        print("  ❌ Ninguna corrida produjo resultados válidos.")

if __name__ == "__main__":
    main()
