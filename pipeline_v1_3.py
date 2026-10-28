import warnings
warnings.filterwarnings("ignore")

import cobra
from cobra.io import load_model
from equilibrator_api import ComponentContribution, Q_
import pandas as pd

def analizar_glucolisis_y_krebs():
    print("🚀 Iniciando Pipeline Metabólico v1.3 (Glucólisis + Ciclo de Krebs)\n")
    
    # 1. Cargar modelo y optimizar (FBA)
    modelo = load_model("textbook")
    solucion = modelo.optimize()
    
    # 2. Configurar eQuilibrator
    cc = ComponentContribution()
    cc.p_h = Q_(7.0)
    cc.p_mg = Q_(3.0)
    cc.ionic_strength = Q_("0.1 M")
    
    # 3. Mapeo Manual Extendido: Glucólisis + Ciclo de Krebs
    mapeo_reacciones = {
        # --- GLUCÓLISIS ---
        'PGI': 'kegg:C00092 = kegg:C00085',
        'PFK': 'kegg:C00085 + kegg:C00002 = kegg:C00354 + kegg:C00008',
        'FBA': 'kegg:C00354 = kegg:C00118 + kegg:C00111',
        'TPI': 'kegg:C00111 = kegg:C00118',
        'GAPD': 'kegg:C00118 + kegg:C00008 = kegg:C00236 + kegg:C00009', 
        'PGK': 'kegg:C00236 + kegg:C00002 = kegg:C00197 + kegg:C00008',
        'PGM': 'kegg:C00197 = kegg:C00631',
        'ENO': 'kegg:C00631 = kegg:C00074',
        'PYK': 'kegg:C00074 + kegg:C00008 = kegg:C00022 + kegg:C00002',
        
        # --- CICLO DE KREBS (TCA) ---
        # Citrato sintasa: Acetil-CoA + OAA + H2O -> Citrato + CoA
        'CS': 'kegg:C00024 + kegg:C00036 + kegg:C00001 = kegg:C00158 + kegg:C00010',
        # Aconitasa (parte A): Citrato -> cis-Aconitato + H2O
        'ACONTa': 'kegg:C00158 = kegg:C00417 + kegg:C00001',
        # Aconitasa (parte B): cis-Aconitato + H2O -> Isocitrato
        'ACONTb': 'kegg:C00417 + kegg:C00001 = kegg:C00311',
        # Isocitrato deshidrogenasa: Isocitrato + NADP+ -> alpha-KG + CO2 + NADPH
        'ICDHyr': 'kegg:C00311 + kegg:C00006 = kegg:C00026 + kegg:C00011 + kegg:C00005',
        # alpha-KG deshidrogenasa: alpha-KG + NAD+ + CoA -> Succinil-CoA + CO2 + NADH
        'AKGDH': 'kegg:C00026 + kegg:C00003 + kegg:C00010 = kegg:C00091 + kegg:C00011 + kegg:C00004',
        # Succinil-CoA sintetasa: Succinil-CoA + ADP + Pi -> Succinato + ATP + CoA
        'SUCOAS': 'kegg:C00091 + kegg:C00008 + kegg:C00009 = kegg:C00042 + kegg:C00002 + kegg:C00010',
        # Succinato deshidrogenasa: Succinato + FAD -> Fumarato + FADH2
        'SUCDi': 'kegg:C00042 + kegg:C00016 = kegg:C00122 + kegg:C01352',
        # Fumarasa: Fumarato + H2O -> Malato
        'FUM': 'kegg:C00122 + kegg:C00001 = kegg:C00149',
        # Malato deshidrogenasa: Malato + NAD+ -> OAA + NADH
        'MDH': 'kegg:C00149 + kegg:C00003 = kegg:C00036 + kegg:C00004',
        
        # --- REACCIONES ANAPLERÓTICAS (Relleno del ciclo) ---
        # PEP carboxilasa: PEP + CO2 -> OAA + Pi
        'PPC': 'kegg:C00074 + kegg:C00011 = kegg:C00036 + kegg:C00009'
    }
    
    resultados = []
    print("Calculando termodinámica (Glucólisis y Krebs)...\n")
    
    for reaccion_id, formula_reaccion in mapeo_reacciones.items():
        if reaccion_id in modelo.reactions:
            reaccion_cobra = modelo.reactions.get_by_id(reaccion_id)
            flujo = solucion.fluxes[reaccion_id]
            
            try:
                reaccion_eq = cc.parse_reaction_formula(formula_reaccion)
                dg_prime = cc.physiological_dg_prime(reaccion_eq)
                viable = "✅ SÍ" if dg_prime.value.magnitude < 0 else "❌ NO"
                
                resultados.append({
                    "ID": reaccion_id,
                    "Nombre": reaccion_cobra.name[:30] + "...",
                    "Flux (FBA)": round(flujo, 4),
                    "ΔG' (kJ/mol)": round(dg_prime.value.magnitude, 2),
                    "Viable?": viable
                })
            except Exception as e:
                print(f"⚠️ Error termodinámico en {reaccion_id}: {e}")
    
    # Mostrar resultados
    df = pd.DataFrame(resultados)
    
    # Filtrar solo las reacciones que la célula está usando activamente
    df_activas = df[df["Flux (FBA)"] != 0].sort_values(by="ΔG' (kJ/mol)")
    
    print(df_activas.to_string(index=False))
    print("\n✅ Análisis completado.")

if __name__ == "__main__":
    analizar_glucolisis_y_krebs()