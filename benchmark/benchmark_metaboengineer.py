"""
Fase 1: benchmark interno GA vs OptKnock-exacto (ambos tuyos).

Uso:
    python benchmark/benchmark_metaboengineer.py --method both --target EX_ac_e --model iML1515.xml
    python benchmark/benchmark_metaboengineer.py --method both --anaerobic --target EX_succ_e --k 1 --model iML1515.xml
"""
import argparse
import sys
from pathlib import Path
from unittest.mock import patch

PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import metabo_engineer as me
from common import (
    load_iml1515, production_metrics, wt_baseline, save_results, Timer,
)


class _FakePrompt:
    responses = []
    @classmethod
    def ask(cls, *a, **k):
        return cls.responses.pop(0)


class _FakeIntPrompt:
    responses = []
    @classmethod
    def ask(cls, *a, **k):
        return int(cls.responses.pop(0))


def _make_patches(captured):
    def _fake_export(modelo_id, quimico, metodo, genes, fitness, seed, parametros_extra=None):
        captured.update({
            "metodo_interno": metodo,
            "genes": list(genes),
            "fitness_interno": float(fitness),
        })
    return [
        patch.object(me, 'Prompt', _FakePrompt),
        patch.object(me, 'IntPrompt', _FakeIntPrompt),
        patch.object(me, '_exportar_resultado_json', _fake_export),
        patch.object(me, 'analisis_profundo', lambda *a, **k: None),
        patch.object(me, 'generar_documento_markdown', lambda *a, **k: None),
        patch.object(me, 'get_organism_id', lambda modelo_id: "83333"),
        patch.object(me, 'traducir_gen_uniprot', lambda gen_id, org: (gen_id, "patched")),
    ]


def _run_with_patches(fn, prompt_responses, int_responses):
    captured = {}
    _FakePrompt.responses = list(prompt_responses)
    _FakeIntPrompt.responses = list(int_responses)
    patches = _make_patches(captured)
    for p in patches:
        p.start()
    try:
        fn()
    finally:
        for p in patches:
            p.stop()
    return captured


def run_ga(model, quimico, gens, pop_size, mutaciones, seed):
    return _run_with_patches(
        lambda: me.disenar_cepa_ia(model, "E. coli", "iML1515", seed=seed),
        prompt_responses=[quimico],
        int_responses=[gens, pop_size, mutaciones],
    )


def run_optknock(model, quimico, max_ko, seed):
    return _run_with_patches(
        lambda: me.optknock_exacto_milp(model, "E. coli", "iML1515", seed=seed),
        prompt_responses=[quimico],
        int_responses=[max_ko],
    )


def _attach_metrics(res, model, args, label):
    if not res or not res.get("genes"):
        print(f"  [!] {label}: sin resultado")
        return res or {}
    print(f"  Genes: {res['genes']}")
    print(f"  Fitness interno (tu GA/OptKnock): {res.get('fitness_interno')}")
    m_90 = production_metrics(model, res["genes"], args.target, args.biomass, 0.9)
    m_100 = production_metrics(model, res["genes"], args.target, args.biomass, 1.0)
    print(f"  [comparable GA]  max@{args.target} a 90% growth = {m_90['max']:.4f}  "
          f"(mutante growth={m_90['growth']:.4f})")
    print(f"  [comparable OK]  min@{args.target} a 100% growth = {m_100['min']:.4f}  "
          f"(ratio={m_100['ratio']:.4f})")
    res["prod_max_at_90pct"] = m_90["max"]
    res["prod_min_at_100pct"] = m_100["min"]
    res["mutant_growth_90pct"] = m_90["growth"]
    res["ratio_at_100pct"] = m_100["ratio"]
    return res


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--model", default=None)
    p.add_argument("--target", default="EX_succ_e")
    p.add_argument("--biomass", default=None)
    p.add_argument("--k", type=int, default=1)
    p.add_argument("--gens", type=int, default=10)
    p.add_argument("--pop", type=int, default=20)
    p.add_argument("--mut", type=int, default=3)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--anaerobic", action="store_true")
    p.add_argument("--method", choices=["ga", "optknock", "both"], default="both")
    args = p.parse_args()

    print("Cargando modelo...")
    model = load_iml1515(args.model)
    ambiente = "anaerobio" if args.anaerobic else "aerobio"
    print(f"Configurando medio: {ambiente}")
    # FIX: asignar el retorno (configurar_ambiente devuelve una COPIA)
    model = me.configurar_ambiente(model, ambiente)

    if args.biomass is None:
        args.biomass = next(r.id for r in model.reactions if "biomass" in r.id.lower())
    if args.target in model.reactions:
        model.reactions.get_by_id(args.target).bounds = (0, 1000.0)

    print(f"  Biomass: {args.biomass}")
    print(f"  Target:  {args.target}")
    print(f"  Genes:   {len(model.genes)}")

    print("\n--- Baseline WT ---")
    wt = wt_baseline(model, args.target, args.biomass)
    print(f"  WT max growth: {model.slim_optimize():.4f} 1/h")
    print(f"  WT {args.target} a 90% growth:  max={wt['at_90pct_growth']['max']:.4f}  "
          f"(baseline para GA)")
    print(f"  WT {args.target} a 100% growth: min={wt['at_100pct_growth']['min']:.4f}  "
          f"(baseline para OptKnock)")
    if wt["at_90pct_growth"]["max"] < 0.01:
        print(f"  ⚠️  WT no produce {args.target} a 90% growth. Cambia --target o --anaerobic.")

    results = {"wt_baseline": wt, "ambiente": ambiente}

    if args.method in ("ga", "both"):
        print(f"\n=== GA (gens={args.gens}, pop={args.pop}, mut={args.mut}, seed={args.seed}) ===")
        with Timer() as t:
            res = run_ga(model, args.target, args.gens, args.pop, args.mut, args.seed)
        res = _attach_metrics(res, model, args, "GA")
        res["time_s"] = t.elapsed
        res["params"] = vars(args)
        print(f"  Tiempo total: {t.elapsed:.1f}s")
        results["ga"] = res

    if args.method in ("optknock", "both"):
        print(f"\n=== OptKnock exacto (k={args.k}) — puede tardar ===")
        if args.k >= 2:
            print("  ⚠️  k>=2 es fuerza bruta (27k+ combos). Con iML1515 puede tardar HORAS.")
            print("  ⚠️  Recomendación: usa k=1 en Fase 1, k>=2 requiere MILP real (Fase 2).")
        with Timer() as t:
            res = run_optknock(model, args.target, args.k, args.seed)
        res = _attach_metrics(res, model, args, "OptKnock")
        res["time_s"] = t.elapsed
        res["params"] = vars(args)
        print(f"  Tiempo total: {t.elapsed:.1f}s")
        results["optknock_exacto"] = res

    save_results("metaboengineer", results)
    print("\n✅ Listo. Ejecuta: python benchmark/compare_results.py")


if __name__ == "__main__":
    main()