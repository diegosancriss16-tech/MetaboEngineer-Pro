import sys
import numpy as np
from pathlib import Path
from collections import Counter
from unittest.mock import patch

PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import metabo_engineer as me
from cobra.io import read_sbml_model

# ============================================================
# PARÁMETROS PRE-REGISTRADOS (no tocar tras correr)
# ============================================================
FRACTION   = 1.0
UMBRAL_BIO = 0.05   # DEBE coincidir con el umbral interno del GA
SEEDS      = [42, 100, 2024, 7, 13]
GENS, POP, K = 20, 50, 5     # VERIFICAR: el 3er prompt debe ser k (no mutación)
KJ060 = ['ldhA', 'adhE', 'ackA', 'focA', 'pflB']   # Jantama et al. 2008

def fitness_fva_min(modelo, genes_ko, rxn_objetivo, rxn_biomasa,
                    fraction=FRACTION, umbral=UMBRAL_BIO, verbose=False):
    modelo_copia = modelo.copy()
    for gen_id in genes_ko:
        if gen_id in modelo_copia.genes:
            modelo_copia.genes.get_by_id(gen_id).knock_out()
    try:
        modelo_copia.objective = rxn_biomasa
        sol_bio = modelo_copia.optimize()
        if sol_bio.status != 'optimal' or sol_bio.objective_value < umbral:
            return -1e6
        max_bio = sol_bio.objective_value
        modelo_copia.reactions.get_by_id(rxn_biomasa).lower_bound = max_bio * fraction
        modelo_copia.objective = rxn_objetivo
        modelo_copia.objective.direction = 'min'
        sol_prod = modelo_copia.optimize()
        if sol_prod.status == 'optimal':
            return max(0.0, sol_prod.objective_value)
        if verbose: print(f"  [warn] LP producto no óptimo: {sol_prod.status}")
    except Exception as e:
        if verbose: print(f"  [warn] excepción: {e}")
    return -1e6

def reacciones_knocked(modelo, genes_ko):
    """Reacciones cuyos bounds cambian al aplicar los KOs (nivel funcional)."""
    ref = {r.id: r.bounds for r in modelo.reactions}
    m = modelo.copy()
    for g in genes_ko:
        if g in m.genes:
            m.genes.get_by_id(g).knock_out()
    return {r.id for r in m.reactions if r.bounds != ref[r.id]}

class _FakePrompt:
    responses = []
    @classmethod
    def ask(cls, *a, **k):
        if not cls.responses:
            raise RuntimeError(f"Prompt inesperado (desincronización): {a} {k}")
        return cls.responses.pop(0)

class _FakeIntPrompt(_FakePrompt):
    @classmethod
    def ask(cls, *a, **k):
        return int(super().ask(*a, **k))

def run_ga(model, target, seed):
    captured = {}
    def _fake_export(modelo_id, quimico, metodo, genes, fitness, seed_, parametros_extra=None):
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
    _FakeIntPrompt.responses = [GENS, POP, K]
    for p in patches: p.start()
    try:
        me.disenar_cepa_ia(model, "E. coli", "iML1515", seed=seed)
    finally:
        for p in patches: p.stop()
    # Detección de desincronización: sobraron prompts sin consumir
    assert not _FakeIntPrompt.responses, \
        f"Prompts sin consumir: {_FakeIntPrompt.responses} — la secuencia de la CLI cambió"
    return captured

def main():
    print("Cargando modelo...")
    model = read_sbml_model("iML1515.xml")
    model = me.configurar_ambiente(model, "anaerobio")
    model.solver = 'glpk'

    target = "EX_succ_e"
    biomass = next(r.id for r in model.reactions if "biomass" in r.id.lower())

    # --- Mapeo con verificación explícita de KJ060 ---
    conteo = Counter(g.name for g in model.genes if g.name)
    dupes = [n for n, c in conteo.items() if c > 1]
    id_por_nombre = {g.name: g.id for g in model.genes if g.name}
    kj_ids, faltan = [], []
    for n in KJ060:
        if n in id_por_nombre: kj_ids.append(id_por_nombre[n])
        else: faltan.append(n)
    assert not faltan, f"Genes de KJ060 ausentes en iML1515: {faltan}"
    if dupes:
        print(f"[warn] nombres de gen duplicados en el modelo: {dupes}")

    print(f"KJ060 -> {dict(zip(KJ060, kj_ids))}")

    # --- Baselines ---
    f_wt    = fitness_fva_min(model, [], target, biomass)
    f_atpb  = fitness_fva_min(model, ['b3731'], target, biomass)
    f_kj060 = fitness_fva_min(model, kj_ids, target, biomass, verbose=True)
    rx_kj   = reacciones_knocked(model, kj_ids)

    print(f"\nBaselines  WT={f_wt:.4f}  atpB={f_atpb:.4f}  KJ060={f_kj060:.4f}")
    print(f"KJ060 bloquea {len(rx_kj)} reacciones")

    # --- Multi-semilla ---
    resultados = []
    print(f"\nGA (k={K}, pop={POP}, gens={GENS}) en {len(SEEDS)} semillas\n")
    for s in SEEDS:
        print(f"--- semilla {s} ---")
        res = run_ga(model, target, s)
        genes, f_int = res.get('genes', []), res.get('fitness', -1e6)
        # AUDITORÍA DE PARIDAD: el fitness interno vs la MISMA función local
        f_rec = fitness_fva_min(model, genes, target, biomass, verbose=True)
        delta = abs(f_int - f_rec) if f_rec > -1e5 and f_int > -1e5 else float('nan')
        ov_g  = set(genes) & set(kj_ids)
        ov_r  = reacciones_knocked(model, genes) & rx_kj
        resultados.append(dict(seed=s, genes=genes, f_int=f_int, f_rec=f_rec,
                               delta=delta, ov_g=len(ov_g), ov_r=len(ov_r)))
        print(f"  genes={genes}\n  fitness_int={f_int:.6f}  fitness_recalculado={f_rec:.6f}  delta={delta}")
        print(f"  overlap_genes={len(ov_g)}/5  overlap_reacciones={len(ov_r)}/{len(rx_kj)}")

    # --- Resumen con criterios pre-registrados ---
    print("\n" + "=" * 60)
    print("CRITERIOS PRE-REGISTRADOS")
    print("=" * 60)
    paridad = all(r['delta'] < 1e-6 for r in resultados)
    mejor = max(resultados, key=lambda r: r['f_rec'])
    n_mejor = sum(1 for r in resultados if r['f_rec'] > f_kj060)
    print(f"A (paridad fitness interno==local):   {'PASS' if paridad else 'FAIL'}")
    print(f"B (recuperación): mejor semilla overlap_gen={mejor['ov_g']}/5, "
          f"overlap_reaccion={mejor['ov_r']}/{len(rx_kj)}")
    print(f"C (calidad): {n_mejor}/{len(resultados)} semillas superan a KJ060 "
          f"(mejor fitness={mejor['f_rec']:.4f} vs KJ060={f_kj060:.4f})")
    if not paridad:
        print("\n>>> A FALLÓ: el GA interno y la función local NO miden lo mismo.")
        print(">>> El experimento NO es interpretable hasta corregir la discrepancia.")

if __name__ == "__main__":
    main()
