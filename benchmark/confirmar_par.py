import sys; from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import cobra, metabo_engineer as me
m = me.cargar_modelo_universal("iML1515")
m = m[0] if isinstance(m, (tuple, list)) else m
m.solver = "glpk"
r = me.configurar_ambiente(m, "anaerobio")
m = r if isinstance(r, cobra.Model) else m
B = "BIOMASS_Ec_iML1515_core_75p37M"
print(f"atpB + ydiO:            {me.evaluar_bacteria(m.copy(), ['b3738','b1695'], 'EX_succ_e', B):.4f}")
print(f"atpB + ydiO + tnaA:      {me.evaluar_bacteria(m.copy(), ['b3738','b1695','b3708'], 'EX_succ_e', B):.4f}")
print(f"ydiO solo (control):     {me.evaluar_bacteria(m.copy(), ['b1695'], 'EX_succ_e', B):.4f}")
