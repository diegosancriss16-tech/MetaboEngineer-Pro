import warnings
import logging
warnings.filterwarnings("ignore")

import os
import requests
import time
import re
import json
from datetime import datetime
from equilibrator_api import ComponentContribution, Q_
from rich.console import Console

from modulo_termodinamica_v2 import construir_reaccion

console = Console()

def get_organism_id(modelo_id: str) -> str:
    if "iYO844" in modelo_id: return "224308"
    if "iCW773" in modelo_id: return "196627"
    if "iMM904" in modelo_id: return "559292"
    if "iEK1008" in modelo_id: return "1773"
    return "83333"

def buscar_info_uniprot_rhea(gen_id, organism_id="83333"):
    if not re.search(r'\d', gen_id): return gen_id, "Nombre comun", "N/A"
    url = f'https://rest.uniprot.org/uniprotkb/search?query=organism_id:{organism_id}+AND+"{gen_id}"&fields=gene_names,protein_name,xref_rhea&format=json'
    try:
        time.sleep(0.2)
        res = requests.get(url, timeout=5)
        if res.status_code == 200 and res.json().get('results'):
            data = res.json()['results'][0]
            nombre = data.get('genes', [{}])[0].get('geneName', {}).get('value', 'Desconocido')
            funcion = data.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', 'Desconocida')
            
            # LEER RHEA IDS CORRECTAMENTE
            rhea_ids = [x.get('id') for x in data.get('uniProtKBCrossReferences', []) if x.get('database') == 'RHEA']
            rhea_str = ", ".join(rhea_ids) if rhea_ids else "No en Rhea"
            
            return nombre, funcion, rhea_str
    except Exception as e:
        logging.warning(f"Consulta UniProt/Rhea fallo: {e}")
    return gen_id, "No encontrado", "N/A"

def generar_documento_markdown(modelo, modelo_id, quimico, mejor_genes, fitness):
    console.print(f"\n[cyan]INFO[/cyan] Generando Reporte Cientifico...")
    organism_id = get_organism_id(modelo_id)

    genes_info = []
    for g in mejor_genes:
        nombre, funcion, rhea = buscar_info_uniprot_rhea(g, organism_id)
        genes_info.append({"id": g, "nombre": nombre, "funcion": funcion, "rhea": rhea})

    if not os.path.exists("reportes"):
        os.makedirs("reportes")

    fecha = datetime.now().strftime("%Y-%m-%d_%H-%M")
    nombre_archivo = f"reportes/REPORTE_{modelo_id}_{quimico}_{fecha}.md"

    with open(nombre_archivo, "w", encoding="utf-8") as f:
        f.write(f"# Reporte de Ingenieria Metabolica: {quimico}\n\n")
        f.write(f"**Fecha:** {datetime.now().strftime('%Y-%m-%d %H:%M')}\n")
        f.write(f"**Organismo Base:** `{modelo_id}`\n")
        f.write(f"**Reacciones totales del modelo:** {len(modelo.reactions)}\n\n")

        f.write("## 1. Resumen Ejecutivo\n\n")
        # FIX: Se cambió "Rendimiento ajustado" por "Flujo Garantizado" para reflejar correctamente la unidad (mmol/gDW/h).
        f.write(f"- **Flujo Garantizado (Fitness):** `{fitness:.4f} mmol/gDW/h`\n")
        f.write(f"- **Genes knock-out:** `{', '.join(mejor_genes)}`\n\n")

        f.write("## 2. Validacion Genomica (UniProt & Rhea)\n\n")
        f.write("| Gen (ID) | Nombre | Funcion Proteica | Rhea ID |\n")
        f.write("| :--- | :--- | :--- | :--- |\n")
        for info in genes_info:
            f.write(f"| `{info['id']}` | **{info['nombre']}** | {info['funcion']} | `{info['rhea']}` |\n")
        f.write("\n")

        # FIX: La sección de "Cinetica Enzimatica" ha sido completamente eliminada.
        # El README establece que la cinética fue removida por falta de datos curados (kcat/Km de BRENDA).

        f.write("## 3. Analisis Termodinamico (eQuilibrator)\n\n")
        f.write("| Reaccion | DG' (kJ/mol) | Viable |\n")
        f.write("| :--- | :--- | :--- |\n")

        cc = ComponentContribution()
        cc.p_h = Q_(7.0); cc.p_mg = Q_(3.0); cc.ionic_strength = Q_("0.1 M")
        
        reacciones_a_probar = [quimico] if quimico in modelo.reactions else []
        for g in mejor_genes:
            if g in modelo.genes:
                for rxn in modelo.genes.get_by_id(g).reactions:
                    if rxn.id not in reacciones_a_probar:
                        reacciones_a_probar.append(rxn.id)

        for r_id in reacciones_a_probar[:10]:
            if r_id.startswith("CAFT") or r_id.startswith("XAN") or r_id.startswith("EX_caft"):
                f.write(f"| `{r_id}` | N/A (Sintetico) | N/A |\n")
                continue
            try:
                formula = construir_reaccion(modelo, modelo.reactions.get_by_id(r_id))
                dg = cc.physiological_dg_prime(cc.parse_reaction_formula(formula)).value.magnitude
                f.write(f"| `{r_id}` | {dg:.2f} | {'SI' if dg < 0 else 'NO'} |\n")
            except Exception as e:
                logging.warning(f"No se pudo calcular DG' de una reaccion: {e}")

        f.write("\n---\n")
        f.write("*Generado automaticamente por MetaboEngineer Pro*\n")

    console.print(f"[green]OK[/green] Reporte guardado: [cyan]{nombre_archivo}[/cyan]")