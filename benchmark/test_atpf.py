import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import metabo_engineer as me
from cobra.io import read_sbml_model

def fitness_fva_min(modelo, genes_ko, gen_objetivo, gen_biomasa, fraction=1.0):
    modelo_copia = modelo.copy()
    for gen_id in genes_ko:
        if gen_id in modelo_copia.genes:
            modelo_copia.genes.get_by_id(gen_id).knock_out()
    try:
        modelo_copia.objective = gen_biomasa
        sol_bio = modelo_copia.optimize()
        if sol_bio.status != 'optimal' or sol_bio.objective_value < 0.05:
            return -1e6
        max_bio = sol_bio.objective_value
        modelo_copia.reactions.get_by_id(gen_biomasa).lower_bound = max_bio * fraction
        modelo_copia.objective = gen_objetivo
        modelo_copia.objective.direction = 'min'
        sol_prod = modelo_copia.optimize()
        if sol_prod.status == 'optimal':
            return max(0.0, sol_prod.objective_value)
    except Exception:
        pass
    return -1e6

def main():
    model = read_sbml_model("iML1515.xml")
    model = me.configurar_ambiente(model, "anaerobio")
    model.solver = 'glpk'
    
    target = "EX_succ_e"
    biomass = next(r.id for r in model.reactions if "biomass" in r.id.lower())
    
    print("=== TEST DE ACOPLAMIENTO INDIVIDUAL ===")
    print(f"atpB (b3731) solo: {fitness_fva_min(model, ['b3731'], target, biomass):.4f}")
    print(f"atpF (b3736) solo: {fitness_fva_min(model, ['b3736'], target, biomass):.4f}")
    print(f"atpD (b3732) solo: {fitness_fva_min(model, ['b3732'], target, biomass):.4f}")

if __name__ == "__main__":
    main()
