import warnings
warnings.filterwarnings("ignore")

import cobra
from cobra.io import load_model
import random
import time

def evaluar_bacteria(modelo_base, genes_apagados, gen_objetivo, gen_biomasa):
    """Evalúa el fitness (producción del químico) de una combinación de mutaciones."""
    # Clonamos el modelo para no modificar el original
    modelo = modelo_base.copy()
    
    # 1. Apagar los genes (knockouts)
    for gen_id in genes_apagados:
        if gen_id in modelo.genes:
            for reaccion in modelo.genes.get_by_id(gen_id).reactions:
                reaccion.bounds = (0, 0)
                
    # 2. Forzar un crecimiento mínimo (para que la célula viva)
    try:
        modelo.reactions.get_by_id(gen_biomasa).bounds = (0.1, 1000.0)
    except:
        return 0.0 
        
    # 3. Cambiar el objetivo a maximizar el químico deseado
    try:
        modelo.objective = gen_objetivo
        solucion = modelo.optimize()
        if solucion.status == 'optimal':
            return solucion.objective_value
    except:
        pass
    return 0.0

def ejecutar_algoritmo():
    print("🚀 Iniciando Motor Genético (CLI) para Diseño de Cepas\n")
    
    # --- PARÁMETROS DEL EXPERIMENTO ---
    MODELO_ID = "textbook" # Usamos el local para evitar fallos de red por ahora
    QUIMICO_OBJETIVO = "EX_succ_e" # Queremos producir Succinato
    GENERACIONES = 10
    TAMANO_POBLACION = 20
    GENES_A_MUTAR = 3
    
    print(f"1. Cargando modelo: {MODELO_ID}...")
    modelo_base = load_model(MODELO_ID)
    
    # Detectar biomasa automáticamente
    reacciones_biomasa = [r.id for r in modelo_base.reactions if "biomass" in r.id.lower()]
    gen_biomasa = reacciones_biomasa[0]
    print(f"   -> Biomasa detectada: {gen_biomasa}")
    
    # Abrir la válvula para que el químico pueda salir
    if QUIMICO_OBJETIVO in modelo_base.reactions:
        modelo_base.reactions.get_by_id(QUIMICO_OBJETIVO).bounds = (0, 1000.0)
        print(f"   -> Objetivo: Maximizar {QUIMICO_OBJETIVO}\n")
    else:
        print(f"❌ Error: {QUIMICO_OBJETIVO} no existe en este modelo.")
        return

    todos_los_genes = [g.id for g in modelo_base.genes]
    print(f"2. Iniciando evolución (Población: {TAMANO_POBLACION}, Generaciones: {GENERACIONES})\n")

    # --- POBLACIÓN INICIAL ---
    poblacion = []
    for _ in range(TAMANO_POBLACION):
        genes_mutados = random.sample(todos_los_genes, GENES_A_MUTAR)
        fitness = evaluar_bacteria(modelo_base, genes_mutados, QUIMICO_OBJETIVO, gen_biomasa)
        poblacion.append({"genes": genes_mutados, "fitness": fitness})

    # --- CICLO EVOLUTIVO ---
    for gen in range(GENERACIONES):
        # Ordenar de mejor a peor
        poblacion.sort(key=lambda x: x["fitness"], reverse=True)
        mejores = poblacion[:5] # Los 5 mejores sobreviven
        
        print(f"Generación {gen+1:02d} | Mejor producción: {mejores[0]['fitness']:.4f} mmol/gDW/h")
        
        # Reproducción y mutación
        nueva_poblacion = mejores.copy()
        while len(nueva_poblacion) < TAMANO_POBLACION:
            padre1 = random.choice(mejores)["genes"]
            padre2 = random.choice(mejores)["genes"]
            
            # Cruzamiento
            hijo = list(set(padre1[:2] + padre2[2:])) 
            if len(hijo) < GENES_A_MUTAR:
                hijo.append(random.choice(todos_los_genes))
                
            fitness_hijo = evaluar_bacteria(modelo_base, hijo, QUIMICO_OBJETIVO, gen_biomasa)
            nueva_poblacion.append({"genes": hijo, "fitness": fitness_hijo})
            
        poblacion = nueva_poblacion
        time.sleep(0.1) # Pequeña pausa para que se lean los mensajes en consola

    # --- RESULTADOS FINALES ---
    print("\n" + "="*50)
    print("🏆 RESULTADOS DE LA EVOLUCIÓN")
    print("="*50)
    poblacion.sort(key=lambda x: x["fitness"], reverse=True)
    mejor_cepa = poblacion[0]
    
    print(f"Producción Máxima de Succinato: {mejor_cepa['fitness']:.4f} mmol/gDW/h")
    print("Genes a apagar (Knockouts) para lograr esto:")
    for g in mejor_cepa["genes"]:
        nombre = modelo_base.genes.get_by_id(g).name
        print(f"  - {g}: {nombre}")
    print("="*50 + "\n")

if __name__ == "__main__":
    ejecutar_algoritmo()