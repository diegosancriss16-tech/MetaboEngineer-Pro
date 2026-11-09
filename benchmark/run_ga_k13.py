import sys
import json
from pathlib import Path
from unittest.mock import patch

PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import metabo_engineer as me
from cobra.io import read_sbml_model

# Simulador de inputs secuenciales para el CLI
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
    
    # Limpiamos y preparamos las respuestas (Pop=50, Gens=20, k=13)
    _FakePrompt.responses = [target]
    _FakeIntPrompt.responses = [50, 20, 13] 
    
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
    print("Cargando modelo...")
    model = read_sbml_model("iML1515.xml")
    model.solver = 'glpk'
    model = me.configurar_ambiente(model, "anaerobio")
    model.solver = 'glpk'

    cache_file = Path("benchmark/results/essential_cache.json")
    cache_file.parent.mkdir(exist_ok=True)
    
    if cache_file.exists():
        print("✅ Usando caché de esencialidad...")
        with open(cache_file, 'r') as f:
            genes_no_esenciales = json.load(f)
            
        def fake_obtener_genes(modelo, biomasa):
            return genes_no_esenciales
            
        patch_esencialidad = patch.object(me, 'obtener_genes_no_esenciales', fake_obtener_genes)
        patch_esencialidad.start()
    else:
        print("Calculando filtro de esencialidad (esto tardará ~20 min)...")
        # Guardamos el caché para futuras corridas
        genes = me.obtener_genes_no_esenciales(model, next(r.id for r in model.reactions if "biomass" in r.id.lower()))
        with open(cache_file, 'w') as f:
            json.dump(genes, f)
        patch_esencialidad = None

    seeds = [42, 100] 
    
    print(f"\n=== LANZANDO GA CON K=13 (pop=50, gens=20) ===")
    for s in seeds:
        print(f"\n--- Corriendo semilla {s} ---")
        res = run_ga(model, "EX_succ_e", seed=s)
        print(f"Genes KO: {res.get('genes')}")
        print(f"Fitness: {res.get('fitness')}")
        
        out_file = Path("benchmark/results/ga_k13_results.jsonl")
        with open(out_file, "a") as f:
            f.write(json.dumps({"seed": s, "genes": res.get('genes'), "fitness": res.get('fitness')}) + "\n")

    if patch_esencialidad:
        patch_esencialidad.stop()

if __name__ == "__main__":
    main()
