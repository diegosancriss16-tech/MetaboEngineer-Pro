import warnings
warnings.filterwarnings("ignore")

import os
import requests
import time
from datetime import datetime
from equilibrator_api import ComponentContribution, Q_
from rich.console import Console

from modulo_termodinamica_v2 import construir_reaccion

console = Console()

DB_CINETICA = {
    'PFK': {'kcat': 150.0, 'Km': 0.1, 'enzima': 'Fosfofructoquinasa'},
    'FBA': {'kcat': 45.0, 'Km': 0.5, 'enzima': 'Fructosa-bisfosfato aldolasa'},
    'GAPD': {'kcat': 120.0, 'Km': 0.2, 'enzima': 'Gliceraldehído-3-fosfato deshidrogenasa'},
    'PYK': {'kcat': 200.0, 'Km': 0.4, 'enzima': 'Piruvato quinasa'},
    'CS': {'kcat': 80.0, 'Km': 0.05, 'enzima': 'Citrato sintasa'},
    'ICDHyr': {'kcat': 60.0, 'Km': 0.03, 'enzima': 'Isocitrato deshidrogenasa'}
}

def buscar_info_uniprot_rhea(gen_id, organism_id="83333"):
    if not (gen_id.startswith('b') and gen_id[1:].isdigit()) and not gen_id.startswith('Y'):
        return gen_id, "Gen ya nombrado", "N/A"
    
    url = f"https://rest.uniprot.org/uniprotkb/search?query=organism_id:{organism_id}+AND+gene_exact:{gen_id}&fields=gene_names,protein_name,xref_rhea&format=json"
    try:
        time.sleep(0.2)
        res = requests.get(url, timeout=5)
        if res.status_code == 200 and res.json().get('results'):
            data = res.json()['results'][0]
            nombre = data.get('genes', [{}])[0].get('geneName', {}).get('value', 'Desconocido')
            funcion = data.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', 'Desconocida')
            rhea_ids = [x.get('id') for x in data.get('uniProtKBCrossReferences', []) if x.get('database') == 'RHEA']
            return nombre, funcion, ", ".join(rhea_ids) if rhea_ids else "No en Rhea"
    except:
        pass
    return gen_id, "No encontrado", "N/A"

def generar_documento_markdown(modelo, modelo_id, quimico, mejor_genes, fitness):
    console.print(f"\n[bold yellow]📄 Generando Reporte Científico Automático...[/bold yellow]")
    
    organism_id = "83333" if "coli" in modelo_id or modelo_id == "textbook" else "559292"
    
    genes_info = []
    for g in mejor_genes:
        nombre, funcion, rhea = buscar_info_uniprot_rhea(g, organism_id)
        genes_info.append({"id": g, "nombre": nombre, "funcion": funcion, "rhea": rhea})
    
    if not os.path.exists("reportes"): 
        os.makedirs("reportes")
        
    fecha = datetime.now().strftime("%Y-%m-%d_%H-%M")
    nombre_archivo = f"reportes/REPORTE_{modelo_id}_{quimico}_{fecha}.md"
    
    with open(nombre_archivo, "w", encoding="utf-8") as f:
        f.write(f"# 🧬 Reporte de Ingeniería Metabólica: Producción de {quimico}\n")
        f.write(f"**Fecha:** {datetime.now().strftime('%Y-%m-%d %H:%M')}\n")
        f.write(f"**Organismo Base:** `{modelo_id}`\n\n")
        
        f.write("## 1. Resumen Ejecutivo\n")
        f.write(f"La cepa diseñada tiene un rendimiento ajustado de **`{fitness:.4f} mmol/gDW/h`**.\n")
        
        f.write("## 2. Validación Genómica (UniProt & Rhea)\n")
        f.write("| Gen (ID) | Nombre | Función Proteica | Rhea ID |\n")
        f.write("| :--- | :--- | :--- | :--- |\n")
        for info in genes_info:
            f.write(f"| `{info['id']}` | **{info['nombre']}** | {info['funcion']} | `{info['rhea']}` |\n")
        
        f.write("\n## 3. Análisis Cinético (BRENDA)\n")
        f.write("| Reacción | Enzima | kcat (1/s) | Km (mM) | Eficiencia |\n")
        f.write("| :--- | :--- | :--- | :--- | :--- |\n")
        for rxn, data in DB_CINETICA.items():
            if rxn in modelo.reactions:
                f.write(f"| `{rxn}` | {data['enzima']} | {data['kcat']} | {data['Km']} | **{data['kcat']/data['Km']:.1f}** |\n")
        
        f.write("\n## 4. Análisis Termodinámico (eQuilibrator)\n")
        f.write("| Reacción | ΔG' (kJ/mol) | Viabilidad |\n")
        f.write("| :--- | :--- | :--- |\n")
        
        cc = ComponentContribution()
        cc.p_h = Q_(7.0); cc.p_mg = Q_(3.0); cc.ionic_strength = Q_("0.1 M")
        for r_id in [r.id for r in modelo.reactions if not r.id.startswith("EX_")][:10]:
            try:
                formula = construir_reaccion(modelo, modelo.reactions.get_by_id(r_id))
                dg = cc.physiological_dg_prime(cc.parse_reaction_formula(formula)).value.magnitude
                f.write(f"| `{r_id}` | {dg:.2f} | {'✅' if dg < 0 else '❌'} |\n")
            except: pass
        
        f.write("\n---\n*Generado automáticamente por MetaboEngineer Pro V4.1*\n")
    
    console.print(f"✅ [bold green]Reporte guardado en:[/bold green] [cyan]{nombre_archivo}[/cyan]")