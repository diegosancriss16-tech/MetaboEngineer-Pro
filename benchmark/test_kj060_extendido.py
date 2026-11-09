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

    FAMILIAS = {
        'lactato  (ldhA)':        ['b1380'],
        'etanol   (adhE+adhP+mhpF)': ['b1241','b1478','b0351'],
        'acetato  (ackA+tdcD+purT)': ['b2296','b1849','b3115'],
        'formato  (focA+focB)':    ['b0904','b2492'],
        'PFL      (pflB+yfiD+alt)': ['b0903','b3114','b3951','b3952'],
    }
    print("=== KJ060 EXTENDIDO: familias completas de isozimas ===\n")
    print(f"{'Diseño':38s}{'n_KO':>4s}{'fitness':>10s}")
    print("-"*52)
    acum = []
    for nombre, genes in FAMILIAS.items():
        acum += genes
        f = fitness_fva_min(model, acum, target, biomass)
        print(f"{nombre:38s}{len(acum):4d}{f:10.4f}")
    print(f"\nBaseline: WT=0.0526, atpB=1.5562")

if __name__ == "__main__":
    main()
