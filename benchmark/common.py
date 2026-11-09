"""Utilidades compartidas por todos los benchmarks."""
import json
import time
from pathlib import Path

import numpy as np
from cobra.io import read_sbml_model
from cobra.flux_analysis import flux_variability_analysis
from cobra.manipulation import knock_out_model_genes

RESULTS_DIR = Path(__file__).parent / "results"


def load_iml1515(model_path=None):
    if model_path and Path(model_path).exists():
        return read_sbml_model(model_path)
    try:
        from cobra.io import load_model
        return load_model("iML1515")
    except Exception as e:
        raise RuntimeError(
            "No se pudo cargar iML1515 desde BiGG. Descarga iML1515.xml y pásalo con --model."
        ) from e


def _safe_float(x, default=0.0):
    if x is None:
        return default
    try:
        v = float(x)
        return v if np.isfinite(v) else default
    except (TypeError, ValueError):
        return default


def production_metrics(model, knockouts, product_rxn, biomass_rxn, fraction=0.9):
    """Métricas al fracción de crecimiento indicada.
    fracción=0.9 → comparable con el GA. fracción=1.0 → comparable con OptKnock."""
    with model:
        knock_out_model_genes(model, knockouts)
        max_growth = model.slim_optimize()
        if max_growth is None or not np.isfinite(max_growth) or max_growth < 1e-6:
            return {"max": 0.0, "min": 0.0, "growth": 0.0, "ratio": 0.0}

        max_growth = float(max_growth)
        bio = model.reactions.get_by_id(biomass_rxn)
        new_bound = fraction * max_growth

        # Si la bound que vamos a asignar es nan/inf, abortar
        if not np.isfinite(new_bound):
            return {"max": 0.0, "min": 0.0, "growth": 0.0, "ratio": 0.0}

        bio.lower_bound = new_bound
        bio.upper_bound = new_bound

        try:
            fva = flux_variability_analysis(
                model, reaction_list=[product_rxn], fraction_of_optimum=0.0
            )
            mx = _safe_float(fva.loc[product_rxn, "maximum"])
            mn = _safe_float(fva.loc[product_rxn, "minimum"])
        except Exception:
            mx = mn = 0.0

        return {
            "max": mx,
            "min": mn,
            "growth": max_growth,
            "ratio": mx / max_growth if max_growth > 0 else 0.0,
        }


def wt_baseline(model, product_rxn, biomass_rxn):
    return {
        "at_90pct_growth": production_metrics(model, [], product_rxn, biomass_rxn, 0.9),
        "at_100pct_growth": production_metrics(model, [], product_rxn, biomass_rxn, 1.0),
    }


def save_results(name, results):
    RESULTS_DIR.mkdir(exist_ok=True)
    path = RESULTS_DIR / f"{name}.json"
    with open(path, "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"✅ Guardado: {path}")
    return path


class Timer:
    def __enter__(self):
        self.start = time.time()
        return self
    def __exit__(self, *args):
        self.elapsed = time.time() - self.start