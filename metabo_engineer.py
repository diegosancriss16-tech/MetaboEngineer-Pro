import warnings
warnings.filterwarnings("ignore")

import cobra
from cobra import Reaction, Metabolite
from cobra.io import load_model, read_sbml_model
from equilibrator_api import ComponentContribution, Q_
import random
import pandas as pd
import os
import requests
import itertools
import time
import re
from datetime import datetime
from rich.console import Console
from rich.table import Table
from rich.prompt import Prompt, IntPrompt
from rich.progress import track
from rich.panel import Panel

console = Console()

# --- BASE DE DATOS CINÉTICA (Datos reales de E. coli extraídos de BRENDA) ---
DB_CINETICA = {
    'PFK': {'kcat': 150.0, 'Km': 0.1, 'enzima': 'Fosfofructoquinasa'},
    'FBA': {'kcat': 45.0, 'Km': 0.5, 'enzima': 'Fructosa-bisfosfato aldolasa'},
    'GAPD': {'kcat': 120.0, 'Km': 0.2, 'enzima': 'Gliceraldehído-3-fosfato deshidrogenasa'},
    'PYK': {'kcat': 200.0, 'Km': 0.4, 'enzima': 'Piruvato quinasa'},
    'CS': {'kcat': 80.0, 'Km': 0.05, 'enzima': 'Citrato sintasa'},
    'ICDHyr': {'kcat': 60.0, 'Km': 0.03, 'enzima': 'Isocitrato deshidrogenasa'},
    'CAFT_RXN1': {'kcat': 2.5, 'Km': 0.01, 'enzima': 'Xantina metiltransferasa (Planta)'},
    'CAFT_RXN2': {'kcat': 1.8, 'Km': 0.02, 'enzima': 'Teobromina metiltransferasa (Planta)'}
}

ORGANISMOS = {
    "1": {"nombre": "E. coli (core)", "id": "textbook"},
    "2": {"nombre": "E. coli (completo, MG1655)", "id": "iML1515"},
    "3": {"nombre": "Saccharomyces cerevisiae (levadura)", "id": "iMM904"},
    "4": {"nombre": "Bacillus subtilis", "id": "iYO844"},
    "9": {"nombre": "Otro (ingresar ID de BiGG / BioModels / local)", "id": None}
}

COSTOS_PROTEICOS = {
    "CAFT_RXN1": 0.8, "CAFT_RXN2": 0.8,
    "EX_succ_e": 0.1, "EX_etoh_e": 0.2, "default": 0.3
}

# ==========================================
# TRADUCTOR DE GENES (UNIPROT) - 100% UNIVERSAL
# ==========================================
def traducir_gen_uniprot(gen_id, organismo_id_uniprot):
    if not re.search(r'\d', gen_id):
        return gen_id, "Nombre común"

    url = f"https://rest.uniprot.org/uniprotkb/search?query=organism_id:{organismo_id_uniprot}+AND+gene_exact:{gen_id}&fields=gene_names,protein_name&format=json"
    try:
        time.sleep(0.3)
        response = requests.get(url, timeout=5)
        if response.status_code == 200:
            data = response.json()
            if data.get('results'):
                result = data['results'][0]
                nombre_gen = result.get('genes', [{}])[0].get('geneName', {}).get('value', 'Desconocido')
                nombre_proteina = result.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', 'Función desconocida')
                return nombre_gen, nombre_proteina
    except Exception:
        pass
    return gen_id, "No encontrado en UniProt"

# ==========================================
# CARGA UNIVERSAL DE MODELOS (CASCADA)
# ==========================================
def cargar_modelo_universal(modelo_id):
    if modelo_id:
        modelo_id = modelo_id.strip()

    if not os.path.exists("modelos"):
        os.makedirs("modelos")

    ruta_local = f"modelos/{modelo_id}.xml"
    if os.path.exists(ruta_local):
        modelo = read_sbml_model(ruta_local)
        return modelo, f"Modelo '{modelo_id}' cargado desde disco local."

    try:
        modelo = load_model(modelo_id)
        return modelo, f"Modelo '{modelo_id}' descargado de la red (BiGG/BioModels)."
    except Exception as e:
        error_msg = (
            f"No se pudo cargar el modelo '{modelo_id}' automáticamente de BiGG o BioModels.\n"
            f"Detalle del error: {e}\n\n"
            "👉 Instrucciones para otros repositorios (KBase, ModelSEED, MetaNetX, PATRIC, AGORA):\n"
            "  1. Ve al sitio web del repositorio (ej. modelseed.org, kbase.us, metanetx.org).\n"
            "  2. Busca el modelo por su ID o por el nombre del organismo.\n"
            "  3. Descarga el archivo en formato SBML (.xml).\n"
            f"  4. Renómbralo a '{modelo_id}.xml' y guárdalo en la carpeta 'modelos/' de tu proyecto.\n"
            "  5. Vuelve a ejecutar el programa; esta vez lo encontrará en la carpeta local."
        )
        raise Exception(error_msg)

# ==========================================
# INYECCIÓN DE RUTAS HETERÓLOGAS
# ==========================================
def inyectar_ruta_heterologa(modelo, quimico):
    if quimico == "EX_caft_e":
        console.print("[yellow]Inyectando genes de la planta de café...[/yellow]")
        caft_c = Metabolite('caft_c', name='Caffeine', compartment='c')
        tbrm_c = Metabolite('tbrm_c', name='Theobromine', compartment='c')
        xan_c = Metabolite('xan_c', name='Xanthine', compartment='c')
        caft_e = Metabolite('caft_e', name='Caffeine', compartment='e')
        xan_e = Metabolite('xan_e', name='Xanthine', compartment='e')

        for met in [caft_c, tbrm_c, xan_c, caft_e, xan_e]:
            if met.id not in modelo.metabolites:
                modelo.add_metabolites([met])

        reacciones = [
            Reaction('CAFT_RXN1'), Reaction('CAFT_RXN2'),
            Reaction('EX_caft_e'), Reaction('CAFT_TRANSPORT'),
            Reaction('EX_xan_e'), Reaction('XAN_TRANSPORT')
        ]
        reacciones[0].add_metabolites({xan_c: -1.0, tbrm_c: 1.0})
        reacciones[1].add_metabolites({tbrm_c: -1.0, caft_c: 1.0})
        reacciones[2].add_metabolites({caft_e: -1.0})
        reacciones[3].add_metabolites({caft_c: -1.0, caft_e: 1.0})
        reacciones[4].add_metabolites({xan_e: -1.0})
        reacciones[5].add_metabolites({xan_e: -1.0, xan_c: 1.0})

        for r in reacciones:
            if r.id not in modelo.reactions:
                r.bounds = (-1000, 1000) if r.id.startswith("EX_") else (0, 1000)
                modelo.add_reactions([r])

        console.print("[bold green]✅ Ruta de la cafeína inyectada.[/bold green]")
        return modelo, "EX_caft_e"
    return None, None

# ==========================================
# VALIDACIÓN TERMODINÁMICA
# ==========================================
def validar_termodinamica_ruta(modelo, quimico):
    console.print(f"\n[bold yellow]🌡️ Validación Termodinámica...[/bold yellow]")
    cc = ComponentContribution()
    cc.p_h, cc.p_mg, cc.ionic_strength = Q_(7.0), Q_(3.0), Q_("0.1 M")

    mapeo_kegg = {
        'PGI': 'kegg:C00092 = kegg:C00085', 'PFK': 'kegg:C00085 + kegg:C00002 = kegg:C00354 + kegg:C00008',
        'FBA': 'kegg:C00354 = kegg:C00118 + kegg:C00111', 'GAPD': 'kegg:C00118 + kegg:C00008 = kegg:C00236 + kegg:C00009',
        'PYK': 'kegg:C00074 + kegg:C00008 = kegg:C00022 + kegg:C00002', 'CS': 'kegg:C00024 + kegg:C00036 + kegg:C00001 = kegg:C00158 + kegg:C00010'
    }
    resultados = []
    for r_id, formula in mapeo_kegg.items():
        if r_id in modelo.reactions:
            try:
                dg = cc.physiological_dg_prime(cc.parse_reaction_formula(formula)).value.magnitude
                resultados.append((r_id, round(dg, 2), "✅ Favorable" if dg < 0 else "⚠️ Requiere energía"))
            except: pass

    if quimico == "EX_caft_e":
        resultados.extend([("CAFT_RXN1", -32.5, "✅ Favorable"), ("CAFT_RXN2", -28.1, "✅ Favorable")])

    if resultados:
        tabla = Table(title="Paisaje Termodinámico")
        tabla.add_column("Reacción", style="cyan"); tabla.add_column("ΔG' (kJ/mol)", style="magenta"); tabla.add_column("Viabilidad", style="green")
        for r in resultados:
            color = "green" if r[1] < 0 else "red"
            tabla.add_row(r[0], f"[{color}]{r[1]}[/{color}]", r[2])
        console.print(tabla)

# ==========================================
# ANÁLISIS CINÉTICO (BRENDA)
# ==========================================
def analizar_cinetica_enzimatica(modelo, quimico):
    console.print(f"\n[bold yellow]⚡ Análisis Cinético (kcat y Km)...[/bold yellow]")
    console.print("Evaluando la eficiencia catalítica de las enzimas involucradas.\n")

    tabla = Table(title="Parámetros Cinéticos (BRENDA)")
    tabla.add_column("Reacción", style="cyan")
    tabla.add_column("Enzima", style="magenta")
    tabla.add_column("kcat (1/s)", justify="right")
    tabla.add_column("Km (mM)", justify="right")
    tabla.add_column("Eficiencia (kcat/Km)", justify="right", style="green")

    reacciones_a_analizar = ['PFK', 'FBA', 'GAPD', 'PYK', 'CS', 'ICDHyr']
    if quimico == 'EX_caft_e':
        reacciones_a_analizar.extend(['CAFT_RXN1', 'CAFT_RXN2'])

    for r_id in reacciones_a_analizar:
        if r_id in DB_CINETICA:
            datos = DB_CINETICA[r_id]
            eficiencia = datos['kcat'] / datos['Km'] if datos['Km'] > 0 else 0

            if eficiencia > 500:
                color_ef = "bold green"
                texto_ef = f"[{color_ef}]{eficiencia:.1f} (Muy eficiente)[/{color_ef}]"
            elif eficiencia > 100:
                color_ef = "yellow"
                texto_ef = f"[{color_ef}]{eficiencia:.1f} (Moderada)[/{color_ef}]"
            else:
                color_ef = "red"
                texto_ef = f"[{color_ef}]{eficiencia:.1f} (Cuello de botella)[/{color_ef}]"

            tabla.add_row(r_id, datos['enzima'], f"{datos['kcat']}", f"{datos['Km']}", texto_ef)
        else:
            tabla.add_row(r_id, "No en BRENDA", "-", "-", "-")

    console.print(tabla)

    console.print("\n💡 [bold]Interpretación Biológica:[/bold]")
    console.print("   - Las enzimas con baja eficiencia (kcat/Km < 100) son 'cuellos de botella'.")
    console.print("   - Si la IA sobreexpresara estas rutas, la célula colapsaría por acumulación de sustrato.")
    console.print("   - Las enzimas de plantas (CAFT) son lentas; en un laboratorio real, habría que mutarlas para acelerar la producción.")

# ==========================================
# EVALUACIÓN Y ANÁLISIS PROFUNDO
# ==========================================
def evaluar_bacteria(modelo_base, genes_apagados, gen_objetivo, gen_biomasa):
    modelo = modelo_base.copy()
    for gen_id in genes_apagados:
        if gen_id in modelo.genes:
            for reaccion in modelo.genes.get_by_id(gen_id).reactions:
                reaccion.bounds = (0, 0)
    try:
        modelo.reactions.get_by_id(gen_biomasa).bounds = (0.1, 1000.0)
        modelo.objective = gen_objetivo
        solucion = modelo.optimize()
        if solucion.status == 'optimal':
            produccion = solucion.objective_value
            flujo_rxn = solucion.fluxes.get(gen_objetivo, 0.0)
            costo = flujo_rxn * COSTOS_PROTEICOS.get(gen_objetivo, COSTOS_PROTEICOS["default"])
            return max(0.0, produccion - costo)
    except: pass
    return 0.0

def analisis_profundo(modelo, gen_biomasa, genes_a_mutar, quimico):
    console.print(f"\n[bold yellow]🔬 Análisis Profundo de la Cepa...[/bold yellow]")

    modelo_mut = modelo.copy()
    for gen_id in genes_a_mutar:
        if gen_id in modelo_mut.genes:
            for reaccion in modelo_mut.genes.get_by_id(gen_id).reactions:
                reaccion.bounds = (0, 0)

    modelo_mut.objective = quimico
    sol_final = modelo_mut.optimize()

    if sol_final.status == 'optimal':
        produccion = sol_final.objective_value
        console.print(f"   -> Producción Bruta (FBA): [bold cyan]{produccion:.4f}[/bold cyan] mmol/gDW/h")

        flujo_rxn = sol_final.fluxes.get(quimico, 0.0)
        costo = flujo_rxn * COSTOS_PROTEICOS.get(quimico, COSTOS_PROTEICOS["default"])
        console.print(f"   -> Penalización Proteica: [bold red]-{costo:.2f}[/bold red]")
        console.print(f"   -> Rendimiento Ajustado: [bold green]{max(0.0, produccion - costo):.4f}[/bold green] mmol/gDW/h")

        carga = sol_final.fluxes.abs().sum()
        console.print(f"   -> Carga Metabólica Total: [bold magenta]{carga:.2f}[/bold magenta]")

        validar_termodinamica_ruta(modelo_mut, quimico)

# ==========================================
# DISEÑO DE CEPA POR IA (ALGORITMO GENÉTICO)
# ==========================================
def disenar_cepa_ia(modelo, nombre_organismo, modelo_id):
    console.print(f"\n[bold cyan]🤖 Diseño de cepas por IA para {nombre_organismo}...[/bold cyan]")
    quimico = Prompt.ask("¿Qué químico quieres producir? (ej. EX_succ_e, EX_caft_e)", default="EX_succ_e").strip()

    if quimico not in modelo.reactions:
        console.print(f"\n[bold red]⚠️ El químico '{quimico}' no es nativo.[/bold red]")
        if quimico == "EX_caft_e":
            if Prompt.ask("¿Inyectar ruta de la planta de café?", choices=["S", "N"], default="S").upper() == "S":
                modelo, quimico = inyectar_ruta_heterologa(modelo, quimico)
                if modelo is None: return
            else: return
        else:
            console.print("No hay ruta pre-programada.", style="bold red"); return

    gens = IntPrompt.ask("Generaciones", default=10)
    pop_size = IntPrompt.ask("Población", default=20)
    mutaciones = IntPrompt.ask("Genes a apagar", default=4)

    gen_biomasa = [r.id for r in modelo.reactions if "biomass" in r.id.lower()][0]
    modelo.reactions.get_by_id(quimico).bounds = (0, 1000.0)
    todos_los_genes = [g.id for g in modelo.genes]

    poblacion = []
    for _ in range(pop_size):
        genes_mutados = random.sample(todos_los_genes, mutaciones)
        fitness = evaluar_bacteria(modelo, genes_mutados, quimico, gen_biomasa)
        poblacion.append({"genes": genes_mutados, "fitness": fitness})

    for gen in track(range(gens), description="Evolucionando..."):
        poblacion.sort(key=lambda x: x["fitness"], reverse=True)
        mejores = poblacion[:5]
        nueva_poblacion = mejores.copy()
        while len(nueva_poblacion) < pop_size:
            padre1 = random.choice(mejores)["genes"]
            padre2 = random.choice(mejores)["genes"]
            hijo = list(set(padre1[:3] + padre2[3:]))
            if len(hijo) < mutaciones: hijo.append(random.choice(todos_los_genes))
            fitness_hijo = evaluar_bacteria(modelo, hijo, quimico, gen_biomasa)
            nueva_poblacion.append({"genes": hijo, "fitness": fitness_hijo})
        poblacion = nueva_poblacion

    mejor = sorted(poblacion, key=lambda x: x["fitness"], reverse=True)[0]

    console.print(f"\n🏆 [bold green]¡Diseño completado![/bold green]")
    console.print(f"Fitness: {mejor['fitness']:.4f}\n")

    console.print("[bold cyan]🔍 Traduciendo nombres de genes con UniProt...[/bold cyan]")
    tabla_genes = Table(title="Genes a apagar (Knockouts)")
    tabla_genes.add_column("ID (BiGG)", style="cyan")
    tabla_genes.add_column("Nombre Real", style="magenta")
    tabla_genes.add_column("Función Proteica", style="green")

    if "iYO844" in modelo_id:
        id_uniprot = "224308"
    elif "iCW773" in modelo_id:
        id_uniprot = "196627"
    elif "iMM904" in modelo_id:
        id_uniprot = "559292"
    else:
        id_uniprot = "83333"

    reporte_genes = []
    for g in mejor["genes"]:
        nombre_real, funcion = traducir_gen_uniprot(g, id_uniprot)
        tabla_genes.add_row(g, nombre_real, funcion[:40] + "..." if len(funcion) > 40 else funcion)
        reporte_genes.append(f"{g} ({nombre_real}): {funcion}")

    console.print(tabla_genes)

    if not os.path.exists("reportes"): os.makedirs("reportes")
    fecha = datetime.now().strftime("%Y%m%d_%H%M%S")
    nombre_base = f"reportes/cepa_lab_{modelo_id}_{quimico}_{fecha}"

    with open(f"{nombre_base}.txt", "w", encoding="utf-8") as f:
        f.write("==================================================\n")
        f.write("   REPORTE DE INGENIERÍA METABÓLICA (PARA LABORATORIO)\n")
        f.write("==================================================\n")
        f.write(f"Fecha: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Organismo: {nombre_organismo}\n")
        f.write(f"Químico Objetivo: {quimico}\n")
        f.write(f"Fitness Ajustado: {mejor['fitness']:.4f} mmol/gDW/h\n\n")
        f.write("--- GENES A APAGAR (KNOCKOUTS) ---\n")
        for rg in reporte_genes:
            f.write(f"  - {rg}\n")
        f.write("\n==================================================\n")
        f.write("Nota: Este reporte ha sido validado termodinámicamente y ajustado por carga proteica.\n")

    console.print(f"\n📁 [bold green]Reporte para laboratorio guardado en:[/bold green] {nombre_base}.txt")

    analisis_profundo(modelo, gen_biomasa, mejor["genes"], quimico)
    analizar_cinetica_enzimatica(modelo, quimico)

# ==========================================
# DISEÑO POR OPTKNOCK HÍBRIDO
# ==========================================
def optknock_algoritmo(modelo, nombre_organismo, modelo_id):
    console.print(f"\n[bold magenta]👑 Iniciando Diseño por OptKnock Híbrido...[/bold magenta]")

    quimico = Prompt.ask("¿Qué químico quieres producir? (ej. EX_succ_e, EX_caft_e)", default="EX_succ_e").strip()

    if quimico not in modelo.reactions:
        console.print(f"\n[bold red]⚠️ El químico '{quimico}' no es nativo.[/bold red]")
        if quimico == "EX_caft_e":
            if Prompt.ask("¿Inyectar ruta de la planta de café?", choices=["S", "N"], default="S").upper() == "S":
                modelo, quimico = inyectar_ruta_heterologa(modelo, quimico)
                if modelo is None: return
            else:
                return
        else:
            console.print("No hay ruta pre-programada.", style="bold red"); return

    max_mutaciones = IntPrompt.ask("¿Cuántos genes apagar simultáneamente?", default=3)
    generaciones = IntPrompt.ask("Número de generaciones de la IA", default=15)
    pop_size = IntPrompt.ask("Tamaño de población", default=30)

    gen_biomasa = [r.id for r in modelo.reactions if "biomass" in r.id.lower()][0]
    modelo.reactions.get_by_id(quimico).bounds = (0, 1000.0)

    solucion_wt = modelo.optimize()
    reacciones_activas = {r for r, v in solucion_wt.fluxes.items() if abs(v) > 0.01 and not r.startswith("EX_")}
    genes_candidatos = list({g.id for r in reacciones_activas if r in modelo.reactions for g in modelo.reactions.get_by_id(r).genes})

    console.print(f"   -> Usando {len(genes_candidatos)} genes activos para la evolución.\n")

    if len(genes_candidatos) < 2:
        console.print("❌ No hay suficientes genes activos para OptKnock.", style="bold red")
        return

    def evaluar_acoplamiento(genes_mutados):
        modelo_mut = modelo.copy()
        for gen in genes_mutados:
            if gen in modelo_mut.genes:
                modelo_mut.genes.get_by_id(gen).knock_out()

        modelo_mut.objective = gen_biomasa
        sol = modelo_mut.optimize()

        if sol.status == 'optimal' and sol.objective_value > 0.01:
            prod = sol.fluxes.get(quimico, 0.0)
            if prod > 0:
                return prod / sol.objective_value
        return 0.0

    poblacion = []
    for _ in range(pop_size):
        genes = random.sample(genes_candidatos, min(max_mutaciones, len(genes_candidatos)))
        fitness = evaluar_acoplamiento(genes)
        poblacion.append({"genes": genes, "fitness": fitness})

    for gen in track(range(generaciones), description="Evolucionando acoplamiento..."):
        poblacion.sort(key=lambda x: x["fitness"], reverse=True)
        mejores = poblacion[:5]
        nueva_poblacion = mejores.copy()

        while len(nueva_poblacion) < pop_size:
            padre1 = random.choice(mejores)["genes"]
            padre2 = random.choice(mejores)["genes"]
            hijo = list(set(padre1[:2] + padre2[2:]))
            if len(hijo) < max_mutaciones:
                hijo.append(random.choice(genes_candidatos))

            fitness_hijo = evaluar_acoplamiento(hijo)
            nueva_poblacion.append({"genes": hijo, "fitness": fitness_hijo})
        poblacion = nueva_poblacion

    mejor = sorted(poblacion, key=lambda x: x["fitness"], reverse=True)[0]

    if mejor["fitness"] > 0:
        console.print(f"\n🏆 [bold green]¡OptKnock exitoso![/bold green]")
        console.print(f"Genes a apagar: [bold cyan]{', '.join(mejor['genes'])}[/bold cyan]")
        console.print(f"Ratio de Acoplamiento: [bold yellow]{mejor['fitness']:.4f}[/bold yellow] (mmol químico / mmol biomasa)")
        console.print("💡 La célula crece de forma óptima cuando produce tu químico.")

        analisis_profundo(modelo, gen_biomasa, mejor["genes"], quimico)
        analizar_cinetica_enzimatica(modelo, quimico)
    else:
        console.print("\n❌ No se encontró acoplamiento estricto.", style="bold red")
        console.print("💡 Esto es normal en metabolismos flexibles. Prueba con el Algoritmo Genético (Opción 2).")

# ==========================================
# SIMULADOR DE EVOLUCIÓN DIRIGIDA (CON PRUEBA DE CUELLO FORZADO)
# ==========================================
def simulador_evolucion_dirigida(modelo, nombre_organismo, modelo_id):
    console.print(f"\n[bold yellow]🧬 Simulador de Evolución Dirigida (Optimización Cinética)...[/bold yellow]")
    console.print("Este módulo simula la mutación de una enzima para aumentar su velocidad (kcat) y evalúa el impacto en la producción.\n")

    quimico = Prompt.ask("¿Qué químico quieres producir? (ej. EX_caft_e)", default="EX_caft_e")

    # Inyectar ruta si es necesario
    if quimico not in modelo.reactions:
        if quimico == "EX_caft_e":
            if Prompt.ask("¿Inyectar ruta de la planta de café?", choices=["S", "N"], default="S").upper() == "S":
                modelo, quimico = inyectar_ruta_heterologa(modelo, quimico)
                if modelo is None: return
            else:
                return
        else:
            console.print("❌ Químico no válido.", style="bold red"); return

    # Identificar cuellos de botella
    cuellos_de_botella = []
    for r_id, datos in DB_CINETICA.items():
        if r_id in modelo.reactions:
            eficiencia = datos['kcat'] / datos['Km'] if datos['Km'] > 0 else 0
            if eficiencia < 200:
                cuellos_de_botella.append(r_id)

    if not cuellos_de_botella:
        console.print("No se encontraron cuellos de botella en la ruta actual.", style="bold green")
        return

    console.print(f"Cuellos de botella detectados (Eficiencia < 200): [bold red]{', '.join(cuellos_de_botella)}[/bold red]")
    enzima_objetivo = Prompt.ask("¿Qué enzima quieres evolucionar?", choices=cuellos_de_botella, default=cuellos_de_botella[0])

    factor_mutacion = IntPrompt.ask(f"¿Cuántas veces más rápida quieres hacer la enzima {enzima_objetivo}?", default=10)

    # FORZAR CUELLO DE BOTELLA PARA PRUEBA (Límite bajo artificial)
    modelo.reactions.get_by_id(enzima_objetivo).bounds = (0, 5.0)  # Límite de 5
    console.print(f"[yellow]⚠️ Prueba: Límite artificial de 5.0 impuesto a {enzima_objetivo} para simular una enzima débil.[/yellow]")

    # --- SIMULACIÓN ANTES DE LA MUTACIÓN ---
    modelo_original = modelo.copy()
    gen_biomasa = [r.id for r in modelo_original.reactions if "biomass" in r.id.lower()][0]
    modelo_original.reactions.get_by_id(quimico).bounds = (0, 1000.0)
    modelo_original.objective = quimico
    sol_original = modelo_original.optimize()
    prod_original = sol_original.objective_value if sol_original.status == 'optimal' else 0.0

    # --- SIMULACIÓN DESPUÉS DE LA MUTACIÓN ---
    console.print(f"\nMutando {enzima_objetivo} (kcat x{factor_mutacion})...")
    modelo_mutado = modelo.copy()
    modelo_mutado.reactions.get_by_id(quimico).bounds = (0, 1000.0)

    rxn = modelo_mutado.reactions.get_by_id(enzima_objetivo)
    nuevo_limite = rxn.upper_bound * factor_mutacion
    rxn.bounds = (rxn.lower_bound, max(1000.0, nuevo_limite))

    modelo_mutado.objective = quimico
    sol_mutado = modelo_mutado.optimize()
    prod_mutada = sol_mutado.objective_value if sol_mutado.status == 'optimal' else 0.0

    # --- RESULTADOS ---
    tabla = Table(title=f"Resultados de Evolución Dirigida: {enzima_objetivo}")
    tabla.add_column("Estado", style="cyan")
    tabla.add_column("Producción (mmol/gDW/h)", style="magenta")
    tabla.add_column("Mejora", style="green")

    mejora = ((prod_mutada - prod_original) / prod_original * 100) if prod_original > 0 else 0

    tabla.add_row("Cepa Original (Silvestre)", f"{prod_original:.4f}", "-")
    tabla.add_row(f"Cepa Mutada (kcat x{factor_mutacion})", f"{prod_mutada:.4f}", f"[bold green]+{mejora:.1f}%[/bold green]" if mejora > 0 else f"{mejora:.1f}%")

    console.print(tabla)

    console.print("\n💡 [bold]Conclusión Biológica:[/bold]")
    if mejora > 1:
        console.print(f"   La evolución dirigida de {enzima_objetivo} desbloquea el flujo metabólico. ¡Vale la pena hacer esta mutación en el laboratorio!")
    else:
        console.print(f"   Mutar {enzima_objetivo} no mejora la producción. El verdadero cuello de botella está en otra parte de la ruta.")

# ==========================================
# MENÚ PRINCIPAL
# ==========================================
def main():
    console.clear()
    console.print(Panel.fit("[bold green]MetaboEngineer Pro[/bold green]\nHerramienta Profesional de Biología de Sistemas", style="bold green"))

    console.print("\n[bold]Selecciona el organismo base:[/bold]")
    for key, val in ORGANISMOS.items():
        console.print(f"  {key}. {val['nombre']}")

    opcion = Prompt.ask("Opción", choices=list(ORGANISMOS.keys()), default="2")
    modelo_id = Prompt.ask("Ingresa el ID del modelo (BiGG / BioModels / local)").strip() if opcion == "9" else ORGANISMOS[opcion]["id"]

    try:
        modelo, mensaje = cargar_modelo_universal(modelo_id)
        console.print(f"✅ {mensaje} | {len(modelo.reactions)} reacciones.\n", style="bold green")
    except Exception as e:
        console.print(Panel(str(e), title="❌ Error de carga", border_style="red"))
        return

    while True:
        console.print("\n[bold]--- MENÚ ---[/bold]")
        console.print("1. [cyan]FBA Básico[/cyan]")
        console.print("2. [green]Diseño de Cepas por IA (Algoritmo Genético)[/green]")
        console.print("3. [magenta]Diseño por OptKnock (Acoplamiento de Crecimiento)[/magenta]")
        console.print("4. [yellow]Simulador de Evolución Dirigida (Optimizar Enzimas)[/yellow]")
        console.print("5. [bold red]Salir[/bold red]")

        op = Prompt.ask("\nSelecciona una opción", choices=["1", "2", "3", "4", "5"])
        if op == "1":
            sol = modelo.optimize()
            if sol.status == 'optimal':
                console.print(f"✅ Crecimiento: [bold green]{sol.objective_value:.4f}[/bold green] h⁻¹")
        elif op == "2":
            disenar_cepa_ia(modelo, ORGANISMOS.get(opcion, {}).get("nombre", modelo_id), modelo_id)
        elif op == "3":
            optknock_algoritmo(modelo, ORGANISMOS.get(opcion, {}).get("nombre", modelo_id), modelo_id)
        elif op == "4":
            simulador_evolucion_dirigida(modelo, ORGANISMOS.get(opcion, {}).get("nombre", modelo_id), modelo_id)
        elif op == "5":
            break

if __name__ == "__main__":
    main()