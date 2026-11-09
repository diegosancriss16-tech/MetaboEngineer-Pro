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

    # Los 13 genes de la cepa KJ060 extendida
    KO_13 = ['b1380','b1241','b1478','b0351','b2296','b1849','b3115','b0904','b2492','b0903','b3114','b3951','b3952']

    print("=== VERIFICACIÓN 1: TASA DE CRECIMIENTO ===")
    wt_growth = model.slim_optimize()
    print(f"WT max growth: {wt_growth:.4f}")

    m_mut = model.copy()
    for g in KO_13:
        if g in m_mut.genes:
            m_mut.genes.get_by_id(g).knock_out()
    m_mut.objective = biomass
    mu_mut = m_mut.slim_optimize()
    
    if mu_mut is None or mu_mut < 0.01:
        print("❌ Mutante 13-KO: INVABLE (crecimiento casi cero)")
    else:
        print(f"Mutante 13-KO max growth: {mu_mut:.4f}  ({mu_mut/wt_growth*100:.1f}% del WT)")

    print("\n=== VERIFICACIÓN 2: KNIFE-EDGE (FVA-min a 90%) ===")
    f100 = fitness_fva_min(model, KO_13, target, biomass, 1.0)
    f90  = fitness_fva_min(model, KO_13, target, biomass, 0.9)
    
    print(f"Fitness 13-KO @ 100% growth: {f100:.4f}")
    print(f"Fitness 13-KO @ 90% growth:  {f90:.4f}")
    
    if f90 < 0.1 * f100:
        print("⚠️ El acoplamiento es extremadamente frágil (colapsa al 90%).")
    elif f90 > 1.0:
        print("✅ El acoplamiento sobrevive a la relajación del crecimiento (es robusto).")

if __name__ == "__main__":
    main()
