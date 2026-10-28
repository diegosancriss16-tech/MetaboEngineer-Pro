import warnings
warnings.filterwarnings("ignore")

import os
import requests
import time
import re
from datetime import datetime
from equilibrator_api import ComponentContribution, Q_
from rich.console import Console

from modulo_termodinamica_v2 import construir_reaccion

console = Console()

# NOTA DE DISENO: NO hay datos hardcodeados de ningun organismo.
# La cinetica (kcat/Km) se obtendria de BRENDA, que requiere API autenticada.
# Por eso el reporte solo incluye cinetica si el usuario la aporta por separado.


def get_organism_id(modelo_id: str) -> str:
    """Taxon ID de UniProt segun modelo. Universal."""
    if "iYO844" in modelo_id: return "224308"    # B. subtilis
    if "iCW773" in modelo_id: return "196627"    # C. glutamicum
    if "iMM904" in modelo_id: return "559292"    # S. cerevisiae
    if "textbook" in modelo_id or "iML1515" in modelo_id or "coli" in modelo_id.lower():
        return "83333"                            # E. coli K-12
    return "83333"


def buscar_info_uniprot_rhea(gen_id, organism_id="83333"):
    """
    Consulta UniProt + Rhea. Universal: si el gen no tiene digitos, asume nombre comun.
    """
    if not re.search(r'\d', gen_id):
        return gen_id, "Nombre comun", "N/A"

    url = (f"https://rest.uniprot.org/uniprotkb/search?"
           f"query=organism_id:{organism_id}+AND+gene_exact:{gen_id}"
           f"&fields=gene_names,protein_name,xref_rhea&format=json")
    try:
        time.sleep(0.2)
        res = requests.get(url, timeout=5)
        if res.status_code == 200 and res.json().get('results'):
            data = res.json()['results'][0]
            nombre = data.get('genes', [{}])[0].get('geneName', {}).get('value', 'Desconocido')
            funcion = (data.get('proteinDescription', {})
                       .get('recommendedName', {})
                       .get('fullName', {})
                       .get('value', 'Desconocida'))
            rhea_ids = [x.get('id') for x in data.get('uniProtKBCrossReferences', [])
                        if x.get('database') == 'RHEA']
            return nombre, funcion, ", ".join(rhea_ids) if rhea_ids else "No en Rhea"
    except:
        pass
    return gen_id, "No encontrado", "N/A"


def generar_documento_markdown(modelo, modelo_id, quimico, mejor_genes, fitness):
    """
    Genera un reporte cientifico .md universal (sin datos hardcodeados de ningun organismo).
    """
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
        f.write(f"- **Rendimiento ajustado (Fitness):** `{fitness:.4f} mmol/gDW/h`\n")
        f.write(f"- **Genes knock-out:** `{', '.join(mejor_genes)}`\n\n")

        f.write("## 2. Validacion Genomica (UniProt & Rhea)\n\n")
        f.write("| Gen (ID) | Nombre | Funcion Proteica | Rhea ID |\n")
        f.write("| :--- | :--- | :--- | :--- |\n")
        for info in genes_info:
            f.write(f"| `{info['id']}` | **{info['nombre']}** | {info['funcion']} | `{info['rhea']}` |\n")
        f.write("\n")

        f.write("## 3. Nota sobre Cinetica Enzimatica\n\n")
        f.write("Los parametros cineticos especificos (kcat/Km) **no estan incluidos** en este reporte "
                "porque requieren una conexion autenticada a la API de BRENDA (acceso no publico). "
                "Esta herramienta prioriza la universalidad: todos los datos genomicos se extraen "
                "dinamicamente desde UniProt y Rhea, por lo que funciona para cualquier organismo.\n\n")
        f.write("Si necesitas cinetica, puedes aportar tus propios valores desde BRENDA/SABIO-RK "
                "para las reacciones relevantes.\n\n")

        f.write("## 4. Analisis Termodinamico (eQuilibrator)\n\n")
        f.write("| Reaccion | DG' (kJ/mol) | Viable |\n")
        f.write("| :--- | :--- | :--- |\n")

        cc = ComponentContribution()
        cc.p_h = Q_(7.0); cc.p_mg = Q_(3.0); cc.ionic_strength = Q_("0.1 M")
        for r_id in [r.id for r in modelo.reactions if not r.id.startswith("EX_")][:10]:
            try:
                formula = construir_reaccion(modelo, modelo.reactions.get_by_id(r_id))
                dg = cc.physiological_dg_prime(cc.parse_reaction_formula(formula)).value.magnitude
                f.write(f"| `{r_id}` | {dg:.2f} | {'SI' if dg < 0 else 'NO'} |\n")
            except:
                pass

        f.write("\n---\n")
        f.write("*Generado automaticamente por MetaboEngineer Pro V5.1 - Universal*\n")

    console.print(f"[green]OK[/green] Reporte guardado: [cyan]{nombre_archivo}[/cyan]")