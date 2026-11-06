import sys
from pathlib import Path
PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

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
        bio_rxn = modelo_copia.reactions.get_by_id(gen_biomasa)
        bio_rxn.lower_bound = max_bio * fraction
        modelo_copia.objective = gen_objetivo
        modelo_copia.objective.direction = 'min'
        sol_prod = modelo_copia.optimize()
        if sol_prod.status == 'optimal':
            return max(0.0, sol_prod.objective_value)
    except Exception:
        pass
    return -1e6

def main():
    print("Cargando modelo...")
    model = read_sbml_model("iML1515.xml")
    model.solver = 'glpk'
    model = me.configurar_ambiente(model, "anaerobio")
    model.solver = 'glpk'
    
    target = "EX_succ_e"
    biomass = next(r.id for r in model.reactions if "biomass" in r.id.lower())
    
    print(f"\n{'='*60}")
    print(f"TEST DE PASAJEROS SILENCIOSOS (fracción 1.0)")
    print(f"{'='*60}")
    
    f_atpa = fitness_fva_min(model, ['b3734'], target, biomass) # atpA solo
    f_pasajeros = fitness_fva_min(model, ['b2528', 'b2256'], target, biomass) # iscA + arnD
    f_wt = fitness_fva_min(model, [], target, biomass) # WT
    
    print(f"atpA solo (b3734)     : {f_atpa:.4f}  (Esperado: ~1.5561)")
    print(f"Pasajeros (iscA+arnD) : {f_pasajeros:.4f}  (Esperado: ~0.0526)")
    print(f"WT (sin KO)           : {f_wt:.4f}  (Baseline)")

if __name__ == "__main__":
    main()
