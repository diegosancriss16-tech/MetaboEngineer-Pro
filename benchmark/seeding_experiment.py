import sys
import math
from pathlib import Path
from unittest.mock import patch

PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import metabo_engineer as me
from cobra.io import read_sbml_model

# Familias de la cepa KJ060 (Jantama 2008) + isozimas
FAMILIAS_KJ060 = {
    "lactato":  ["b1380"],
    "etanol":   ["b1241", "b1478", "b0351"],
    "acetato":  ["b2296", "b1849", "b3115"],
    "formato":  ["b0904", "b2492"],
    "pfl":      ["b0903", "b3114", "b3951", "b3952"],
}

def construir_semillas_kj060(fracciones=(1.0, 0.75, 0.5, 0.25)):
    camino = [g for familia in FAMILIAS_KJ060.values() for g in familia]
    semillas, vistas = [], set()
    for f in fracciones:
        n = max(1, math.ceil(f * len(camino)))
        semilla = camino[:min(n, len(camino))]
        clave = tuple(semilla)
        if clave not in vistas:
            vistas.add(clave)
            semillas.append(list(semilla))
    return semillas

class _FakePrompt:
    responses = []
    @classmethod
    def ask(cls, *a, **k): return cls.responses.pop(0)

class _FakeIntPrompt:
    responses = []
    @classmethod
    def ask(cls, *a, **k): return int(cls.responses.pop(0))

def run_ga_with_seeds(model, target, seed, seeds_to_inject):
    captured = {}
    def _fake_export(modelo_id, quimico, metodo, genes, fitness, seed_, parametros_extra=None):
        captured.update({"genes": list(genes), "fitness": float(fitness)})
    
    _FakePrompt.responses = [target]
    _FakeIntPrompt.responses = [20, 50, 13] # gens=20, pop=50, k=13
    
    # Parcheamos la función _evolucionar para inyectar semillas
    original_evolucionar = me._evolucionar
    def _evolucionar_con_semillas(pool_genes, k, generaciones, pop_size, fitness_fn, **kwargs):
        return original_evolucionar(
            pool_genes, k, generaciones, pop_size, fitness_fn, 
            semillas=seeds_to_inject, **kwargs
        )

    patches = [
        patch.object(me, 'Prompt', _FakePrompt),
        patch.object(me, 'IntPrompt', _FakeIntPrompt),
        patch.object(me, '_exportar_resultado_json', _fake_export),
        patch.object(me, 'analisis_profundo', lambda *a, **k: None),
        patch.object(me, 'generar_documento_markdown', lambda *a, **k: None),
        patch.object(me, 'get_organism_id', lambda m: "83333"),
        patch.object(me, 'traducir_gen_uniprot', lambda g, o: (g, "patched")),
        patch.object(me, '_evolucionar', _evolucionar_con_semillas),
    ]
    for p in patches: p.start()
    try:
        me.disenar_cepa_ia(model, "E. coli", "iML1515", seed=seed)
    finally:
        for p in patches: p.stop()
    return captured

def main():
    print("="*70)
    print(" 🧬 EXPERIMENTO DE SEEDING - ¿Se puede romper el valle de aptitud?")
    print("="*70)

    m = read_sbml_model("iML1515.xml")
    m.solver = 'glpk'
    m = me.configurar_ambiente(m, "anaerobio")
    m.solver = 'glpk'

    semillas = construir_semillas_kj060()
    print(f"Generando {len(semillas)} semillas (escalones del valle)...")
    for i, s in enumerate(semillas):
        print(f"  Semilla {i+1}: {len(s)} genes")

    seeds = [42, 100]
    for s in seeds:
        print(f"\n--- Corriendo GA Sembrado (Semilla original={s}) ---")
        res = run_ga_with_seeds(m, "EX_succ_e", seed=s, seeds_to_inject=semillas)
        print(f"  Genes KO: {res.get('genes')}")
        print(f"  Fitness: {res.get('fitness')}")

if __name__ == "__main__":
    main()
