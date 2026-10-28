import argparse
import warnings
warnings.filterwarnings("ignore")

import cobra
from cobra.io import load_model
import random
import sys

def ejecutar_cli(modelo_id, quimico, genes_a_apagar, generaciones=10):
    print(f"🚀 Iniciando Modo Experto CLI...")
    print(f"Cargando modelo: {modelo_id}")
    
    try:
        if modelo_id in ["textbook", "e_coli_core"]:
            modelo = load_model("textbook")
        else:
            modelo = load_model(modelo_id)
    except Exception as e:
        print(f"❌ Error al cargar el modelo: {e}")
        sys.exit(1)

    gen_biomasa = [r.id for r in modelo.reactions if "biomass" in r.id.lower()][0]
    
    if quimico not in modelo.reactions:
        print(f"❌ Error: El químico '{quimico}' no existe en este modelo.")
        sys.exit(1)

    modelo.reactions.get_by_id(quimico).bounds = (0, 1000.0)
    todos_los_genes = [g.id for g in modelo.genes]
    
    print(f"Buscando maximizar {quimico} apagando {genes_a_apagar} genes con {generaciones} generaciones...")

    # Lógica simple del algoritmo genético (sin interfaz de usuario)
    poblacion = []
    for _ in range(20):
        genes = random.sample(todos_los_genes, genes_a_apagar)
        modelo_mut = modelo.copy()
        for g in genes: modelo_mut.genes.get_by_id(g).knock_out()
        modelo_mut.reactions.get_by_id(gen_biomasa).bounds = (0.1, 1000.0)
        modelo_mut.objective = quimico
        sol = modelo_mut.optimize()
        fit = sol.objective_value if sol.status == 'optimal' else 0.0
        poblacion.append({"genes": genes, "fitness": fit})

    for _ in range(generaciones):
        poblacion.sort(key=lambda x: x["fitness"], reverse=True)
        mejores = poblacion[:5]
        nueva = mejores.copy()
        while len(nueva) < 20:
            p1, p2 = random.choice(mejores)["genes"], random.choice(mejores)["genes"]
            hijo = list(set(p1[:2] + p2[2:]))
            if len(hijo) < genes_a_apagar: hijo.append(random.choice(todos_los_genes))
            modelo_mut = modelo.copy()
            for g in hijo: modelo_mut.genes.get_by_id(g).knock_out()
            modelo_mut.reactions.get_by_id(gen_biomasa).bounds = (0.1, 1000.0)
            modelo_mut.objective = quimico
            sol = modelo_mut.optimize()
            fit = sol.objective_value if sol.status == 'optimal' else 0.0
            nueva.append({"genes": hijo, "fitness": fit})
        poblacion = nueva

    mejor = sorted(poblacion, key=lambda x: x["fitness"], reverse=True)[0]
    print(f"\n✅ Análisis Terminado.")
    print(f"Mejor producción encontrada: {mejor['fitness']:.4f} mmol/gDW/h")
    print(f"Genes a apagar: {', '.join(mejor['genes'])}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="MetaboEngineer CLI Experto")
    parser.add_argument("--modelo", type=str, default="textbook", help="ID del modelo (ej. textbook, iML1515)")
    parser.add_argument("--quimico", type=str, required=True, help="ID del químico (ej. EX_succ_e)")
    parser.add_argument("--genes", type=int, default=3, help="Número de genes a apagar")
    parser.add_argument("--generaciones", type=int, default=10, help="Generaciones de la IA")
    
    args = parser.parse_args()
    ejecutar_cli(args.modelo, args.quimico, args.genes, args.generaciones)