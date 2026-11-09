import sys
import random
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import metabo_engineer as me
from cobra.io import read_sbml_model
from cobra.manipulation import knock_out_model_genes

def fitness_fva_min(modelo, genes_ko, gen_objetivo, gen_biomasa, fraction=0.9):
    """Réplica EXACTA de la función evaluar_bacteria modificada."""
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
    except Exception as e:
        return f"ERROR: {e}"
    return -1e6

def main():
    print("Cargando modelo...")
    model = read_sbml_model("iML1515.xml")
    model.solver = 'glpk'
    model = me.configurar_ambiente(model, "anaerobio")
    model.solver = 'glpk'
    
    target = "EX_succ_e"
    biomass = next(r.id for r in model.reactions if "biomass" in r.id.lower())
    
    # Gen conocido-bueno (atpB)
    ko_atpb = ["b3731"]
    
    # WT (sin KO)
    ko_wt = []
    
    # Random 3 (sin atpB)
    genes_no_esenciales = [g.id for g in model.genes if g.id != "b3731"]
    random.seed(42)
    ko_random = random.sample(genes_no_esenciales, 3)
    
    print(f"\n{'='*60}")
    print(f"DIAGNÓSTICO FVA-MIN (biomass={biomass}, target={target}, fraction=0.9)")
    print(f"{'='*60}\n")
    
    f_atpb = fitness_fva_min(model, ko_atpb, target, biomass, 0.9)
    f_wt = fitness_fva_min(model, ko_wt, target, biomass, 0.9)
    f_rand = fitness_fva_min(model, ko_random, target, biomass, 0.9)
    
    print(f"fitness(['b3731'] atpB) : {f_atpb}")
    print(f"fitness([] WT)           : {f_wt}")
    print(f"fitness(random 3: {ko_random}): {f_rand}")
    
    print(f"\n{'='*60}")
    print("DIAGNÓSTICO:")
    print(f"{'='*60}")
    
    if isinstance(f_atpb, str):
        print(f"🐛 ERROR en atpB: {f_atpb}")
    elif abs(f_atpb - 0.0473) < 1e-3:
        print("🐛 BUG CONFIRMADO: el fitness de atpB es idéntico al plateau (0.0473).")
        print("   evaluar_bacteria NO está aplicando los knockouts correctamente.")
        print("   El 'paisaje plano' era un artefacto de medición.")
    elif f_atpb > 0.5:
        print(f"✅ El paisaje SÍ tiene señal: atpB da {f_atpb:.4f}.")
        print("   El GA está explorando mal (no es bug, es limitación del algoritmo).")
    else:
        print(f"⚠️ Zona gris: atpB da {f_atpb:.4f}. Investigar.")
    
    if abs(f_wt - 0.0473) < 1e-3:
        print(f"\n🔴 El WT (sin KO) da 0.0473 = el mismo valor que el GA.")
        print("   El GA NO HA MEJORADO NADA sobre no hacer nada.")
    
    if abs(f_rand - 0.0473) < 1e-3 and abs(f_atpb - 0.0473) < 1e-3:
        print("\n🔴 CONFIRMACIÓN DOBLE: atpB, WT y random3 dan el mismo fitness.")
        print("   El fitness es INDEPENDIENTE de los knockouts → BUG de medición.")

if __name__ == "__main__":
    main()
