import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))
import metabo_engineer as me
from cobra.io import read_sbml_model

def fitness_fva_min(modelo, genes_ko, rxn_obj, rxn_bio, fraction=1.0):
    m = modelo.copy()
    for g in genes_ko:
        if g in m.genes:
            m.genes.get_by_id(g).knock_out()
    try:
        m.objective = rxn_bio
        s = m.optimize()
        if s.status != 'optimal' or s.objective_value < 0.05: return -1e6
        m.reactions.get_by_id(rxn_bio).lower_bound = s.objective_value * fraction
        m.objective = rxn_obj
        m.objective.direction = 'min'
        p = m.optimize()
        return max(0.0, p.objective_value) if p.status == 'optimal' else -1e6
    except Exception:
        return -1e6

def main():
    model = read_sbml_model("iML1515.xml")
    model = me.configurar_ambiente(model, "anaerobio")
    model.solver = 'glpk'
    target, biomass = "EX_succ_e", next(r.id for r in model.reactions if "biomass" in r.id.lower())

    NOMBRES = {'b3736':'atpF','b4238':'nrdD','b2093':'gatB','b2417':'crr','b2149':'mglA'}
    FULL = ['b3736','b4238','b2093','b2417','b2149']  # diseño semilla 100
    
    print(f"Completo: {fitness_fva_min(model, FULL, target, biomass):.4f}\n")
    print("=== LEAVE-ONE-OUT (quien caiga al baseline ES el que aporta) ===")
    
    for g in FULL:
        resto = [x for x in FULL if x != g]
        f = fitness_fva_min(model, resto, target, biomass)
        print(f"  sin {NOMBRES[g]:6s}: {f:.4f}")
        
    print(f"\nBaseline atpF solo: {fitness_fva_min(model, ['b3736'], target, biomass):.4f}")

if __name__ == "__main__":
    main()
