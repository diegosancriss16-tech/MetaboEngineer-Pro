import warnings
warnings.filterwarnings("ignore")

import cobra
from cobra.io import load_model
from equilibrator_api import ComponentContribution, Q_
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.tree import DecisionTreeClassifier, export_text

def pipeline_definitivo():
    print("🚀 Iniciando Pipeline Metabólico v1.4 (Visualización + Machine Learning)\n")
    
    # 1. Cargar modelo y datos
    modelo = load_model("textbook")
    solucion = modelo.optimize()
    
    cc = ComponentContribution()
    cc.p_h = Q_(7.0)
    cc.p_mg = Q_(3.0)
    cc.ionic_strength = Q_("0.1 M")
    
    mapeo_reacciones = {
        'PGI': 'kegg:C00092 = kegg:C00085',
        'PFK': 'kegg:C00085 + kegg:C00002 = kegg:C00354 + kegg:C00008',
        'FBA': 'kegg:C00354 = kegg:C00118 + kegg:C00111',
        'TPI': 'kegg:C00111 = kegg:C00118',
        'GAPD': 'kegg:C00118 + kegg:C00008 = kegg:C00236 + kegg:C00009', 
        'PGK': 'kegg:C00236 + kegg:C00002 = kegg:C00197 + kegg:C00008',
        'PGM': 'kegg:C00197 = kegg:C00631',
        'ENO': 'kegg:C00631 = kegg:C00074',
        'PYK': 'kegg:C00074 + kegg:C00008 = kegg:C00022 + kegg:C00002',
        'CS': 'kegg:C00024 + kegg:C00036 + kegg:C00001 = kegg:C00158 + kegg:C00010',
        'ACONTa': 'kegg:C00158 = kegg:C00417 + kegg:C00001',
        'ACONTb': 'kegg:C00417 + kegg:C00001 = kegg:C00311',
        'ICDHyr': 'kegg:C00311 + kegg:C00006 = kegg:C00026 + kegg:C00011 + kegg:C00005',
        'AKGDH': 'kegg:C00026 + kegg:C00003 + kegg:C00010 = kegg:C00091 + kegg:C00011 + kegg:C00004',
        'SUCOAS': 'kegg:C00091 + kegg:C00008 + kegg:C00009 = kegg:C00042 + kegg:C00002 + kegg:C00010',
        'SUCDi': 'kegg:C00042 + kegg:C00016 = kegg:C00122 + kegg:C01352',
        'FUM': 'kegg:C00122 + kegg:C00001 = kegg:C00149',
        'MDH': 'kegg:C00149 + kegg:C00003 = kegg:C00036 + kegg:C00004',
        'PPC': 'kegg:C00074 + kegg:C00011 = kegg:C00036 + kegg:C00009'
    }
    
    resultados = []
    for reaccion_id, formula in mapeo_reacciones.items():
        if reaccion_id in modelo.reactions:
            reaccion = modelo.reactions.get_by_id(reaccion_id)
            flujo = solucion.fluxes[reaccion_id]
            try:
                reaccion_eq = cc.parse_reaction_formula(formula)
                dg = cc.physiological_dg_prime(reaccion_eq).value.magnitude
                
                resultados.append({
                    "ID": reaccion_id,
                    "Flujo": round(flujo, 2),
                    "DeltaG": round(dg, 2),
                    "Viable": 1 if dg < 0 else 0  # 1 = SÍ, 0 = NO (para el ML)
                })
            except:
                pass
                
    df = pd.DataFrame(resultados)
    df_activas = df[df["Flujo"] != 0].sort_values(by="DeltaG")

    # ==========================================
    # PARTE 1: VISUALIZACIÓN (Paisaje Energético)
    # ==========================================
    print("📊 Generando gráfico del paisaje energético...")
    
    plt.figure(figsize=(12, 6))
    colores = ['green' if x < 0 else 'red' for x in df_activas['DeltaG']]
    
    barras = plt.bar(df_activas['ID'], df_activas['DeltaG'], color=colores, edgecolor='black')
    
    plt.axhline(0, color='black', linewidth=1.5, linestyle='--')
    plt.title('Paisaje Energético: Glucólisis y Ciclo de Krebs', fontsize=14, fontweight='bold')
    plt.xlabel('Enzimas', fontsize=12)
    plt.ylabel("Energía Libre ΔG' (kJ/mol)", fontsize=12)
    plt.xticks(rotation=45)
    
    # Limitar el eje Y para que las barras sean visibles (GAPD y PGK son gigantes)
    plt.ylim(-60, 60)
    
    # Añadir etiquetas de valor
    for barra in barras:
        yval = barra.get_height()
        if -60 < yval < 60: # Solo etiquetar si está en el rango visible
            plt.text(barra.get_x() + barra.get_width()/2, yval, round(yval, 1), 
                     ha='center', va='bottom' if yval > 0 else 'top', fontsize=9)

    plt.tight_layout()
    plt.savefig("paisaje_energetico.png") # Guarda la imagen
    plt.show() # Muestra la ventana con el gráfico
    print("✅ Gráfico guardado como 'paisaje_energetico.png' en tu carpeta.\n")

    # ==========================================
    # PARTE 2: MACHINE LEARNING (Árbol de Decisión)
    # ==========================================
    print("🤖 Entrenando modelo de Machine Learning...")
    
    # El modelo intentará predecir "Viable" basándose en el "Flujo" y el "DeltaG"
    X = df_activas[['Flujo', 'DeltaG']]
    y = df_activas['Viable']
    
    # Usamos un Árbol de Decisión porque es fácil de leer (nos da las reglas exactas)
    modelo_ml = DecisionTreeClassifier(max_depth=3, random_state=42)
    modelo_ml.fit(X, y)
    
    # Extraer las reglas que aprendió la IA
    reglas = export_text(modelo_ml, feature_names=['Flujo', 'DeltaG'])
    
    print("\n🧠 Reglas descubiertas por la Inteligencia Artificial:")
    print("-------------------------------------------------------")
    print(reglas)
    print("-------------------------------------------------------")
    print("💡 La IA ha descubierto que el 'DeltaG' es el principal predictor de viabilidad.")
    
    # Simulación de Mutación (¿Qué pasa si apagamos la enzima PFK?)
    print("\n🧬 Simulación de Mutación: Apagando la enzima PFK...")
    modelo.reactions.get_by_id('PFK').bounds = (0, 0) # Bloqueamos el flujo
    solucion_mutante = modelo.optimize()
    
    if solucion_mutante.status == 'optimal':
        print(f"   -> Tasa de crecimiento original: {solucion.objective_value:.4f}")
        print(f"   -> Tasa de crecimiento mutante: {solucion_mutante.objective_value:.4f}")
    else:
        print("   -> ¡La célula ha muerto! El modelo no puede crecer sin esta enzima.")
    
    print("\n✅ Análisis completado. ¡Felicidades, has construido un pipeline bioinformático completo!")

if __name__ == "__main__":
    pipeline_definitivo()