"""Lee benchmark/results/*.json y genera tabla comparativa."""
import json
from pathlib import Path
import pandas as pd

RESULTS = Path(__file__).parent / "results"


def _load(name):
    p = RESULTS / f"{name}.json"
    return json.loads(p.read_text()) if p.exists() else None


def _row(label, r):
    return {
        "Método": label,
        "Genes KO": ", ".join(r.get("genes") or []) or "—",
        "N KO": len(r.get("genes") or []),
        "Fitness interno": r.get("fitness_interno"),
        "Max prod @90% (GA-comparable)": r.get("prod_max_at_90pct"),
        "Min prod @100% (OK-comparable)": r.get("prod_min_at_100pct"),
        "Ratio @100%": r.get("ratio_at_100pct"),
        "Crec. mutante (1/h)": r.get("mutant_growth_90pct"),
        "Tiempo (s)": r.get("time_s"),
    }


def main():
    rows = []
    me_res = _load("metaboengineer")
    if not me_res:
        print("⚠️  Sin metaboengineer.json — corre primero benchmark_metaboengineer.py")
        return

    wt = me_res.get("wt_baseline", {})
    if wt:
        rows.append({
            "Método": f"WT baseline ({me_res.get('ambiente', '?')})",
            "Genes KO": "—",
            "N KO": 0,
            "Fitness interno": None,
            "Max prod @90% (GA-comparable)": wt.get("at_90pct_growth", {}).get("max"),
            "Min prod @100% (OK-comparable)": wt.get("at_100pct_growth", {}).get("min"),
            "Ratio @100%": wt.get("at_100pct_growth", {}).get("ratio"),
            "Crec. mutante (1/h)": wt.get("at_90pct_growth", {}).get("growth"),
            "Tiempo (s)": 0,
        })

    if "ga" in me_res:
        rows.append(_row("GA (MetaboEngineer)", me_res["ga"]))
    if "optknock_exacto" in me_res:
        rows.append(_row("OptKnock exacto (bruto)", me_res["optknock_exacto"]))

    df = pd.DataFrame(rows)
    print("\n=== TABLA COMPARATIVA ===\n")
    print(df.to_string(index=False))
    out = RESULTS / "comparison.csv"
    df.to_csv(out, index=False)
    print(f"\n✅ Guardado: {out}")


if __name__ == "__main__":
    main()