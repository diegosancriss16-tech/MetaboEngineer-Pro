import sys
from pathlib import Path
PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import metabo_engineer as me
from cobra.io import read_sbml_model

def fitness_fva_min(modelo, genes_ko, gen_objetivo, gen_biomasa, fraction):
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
    
    print(f"\n{'='*70}")
    print(f"TEST 1: Reproducir atpB a fracción 1.0 (100%) vs 0.9 (90%)")
    print(f"{'='*70}")
    print(f"WT @ 1.0:    {fitness_fva_min(model, [], target, biomass, 1.0):.4f}")
    print(f"WT @ 0.9:    {fitness_fva_min(model, [], target, biomass, 0.9):.4f}")
    print(f"atpB @ 1.0:  {fitness_fva_min(model, ['b3731'], target, biomass, 1.0):.4f}  (Esperado: ~1.5562)")
    print(f"atpB @ 0.9:  {fitness_fva_min(model, ['b3731'], target, biomass, 0.9):.4f}")

    print(f"\n{'='*70}")
    print(f"TEST 2 y 3: Knockouts canónicos de succinato")
    print(f"ldhA=b1383, adhE=b1241, pta=b0904, pflB=b0902")
    print(f"{'='*70}")
    canonicos = [
        (['b1383'], 'ldhA'),
        (['b1241'], 'adhE'),
        (['b0904'], 'pta'),
        (['b1383', 'b1241'], 'ldhA+adhE'),
        (['b1383', 'b1241', 'b0904'], 'ldhA+adhE+pta'),
        (['b1383', 'b1241', 'b0904', 'b0902'], 'ldhA+adhE+pta+pflB')
    ]
    for combo, label in canonicos:
        f90 = fitness_fva_min(model, combo, target, biomass, 0.9)
        f100 = fitness_fva_min(model, combo, target, biomass, 1.0)
        print(f"{label:25} @0.9: {f90:7.4f}  |  @1.0: {f100:7.4f}")

if __name__ == "__main__":
    main()
