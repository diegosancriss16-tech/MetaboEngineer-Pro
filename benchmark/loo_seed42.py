# -*- coding: utf-8 -*-
"""Leave-one-out sobre la solucion seed=42 del GA k=13 (fitness 1.6444)."""
import sys, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import cobra
import metabo_engineer as me

BIOMASA = "BIOMASS_Ec_iML1515_core_75p37M"
SOL42 = ["b2521","b3708","b1396","b2867","b4099","b1227","b1850","b3738",
         "b1912","b2787","b2093","b1517","b1695"]
ATPB = "b3738"

crudo = me.cargar_modelo_universal("iML1515")
m = crudo[0] if isinstance(crudo, (tuple, list)) else crudo
m.solver = "glpk"
r = me.configurar_ambiente(m, "anaerobio")
m = r if isinstance(r, cobra.Model) else m

def fit(genes):
    return me.evaluar_bacteria(m.copy(), genes, "EX_succ_e", BIOMASA)

base_full = fit(SOL42)
base_atpb = fit([ATPB])
print(f"Referencias: 13 KOs (seed42) = {base_full:.4f} | atpB solo = {base_atpb:.4f}\n")

print(f"{'Gen quitado':<12}{'Fitness':>9}   {'Efecto'}")
print("-" * 48)
efectos = []
for g in SOL42:
    restantes = [x for x in SOL42 if x != g]
    v = fit(restantes)
    delta = base_full - v
    tag = "NEUTRAL (pasajero)" if delta < 1e-4 else f"CONTRIBUYE (-{delta:.4f})"
    efectos.append((g, v, delta))
    print(f"{g:<12}{v:>9.4f}   {tag}")

contribuyen = [g for g, v, d in efectos if d >= 1e-4]
print("-" * 48)
print(f"Genes que sostienen la mejora: {contribuyen if contribuyen else 'NINGUNO - todo es atpB'}")

out = {"fitness_completo": base_full, "fitness_atpB_solo": base_atpb,
       "genes_contribuyentes": contribuyen,
       "loo": [{"gen": g, "fitness_sin_el": round(v, 6), "delta": round(d, 6)} for g, v, d in efectos]}
Path("benchmark/results/loo_seed42.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
print("\nGuardado: benchmark/results/loo_seed42.json")
