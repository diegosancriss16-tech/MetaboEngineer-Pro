import warnings
warnings.filterwarnings("ignore")

import cobra
from cobra import Reaction, Metabolite
from cobra.io import load_model, read_sbml_model
import random
import os
import requests
import time
import re
from datetime import datetime
from rich.console import Console
from rich.table import Table
from rich.prompt import Prompt, IntPrompt
from rich.progress import track
from rich.panel import Panel

console = Console()

# --- IMPORTACION PROTEGIDA DE MODULOS EXTERNOS ---
MODULOS_EXTERNOS_DISPONIBLES = False
try:
    from visualizador_universal import generar_mapa_universal
    from modulo_termodinamica_v2 import analizar_termodinamica_universal, construir_reaccion
    MODULOS_EXTERNOS_DISPONIBLES = True
except ImportError as e:
    console.print(f"[yellow]WARN[/yellow] Modulos externos no disponibles ({e}). Modo reducido.")

REPORTE_DISPONIBLE = False
try:
    from reporte_investigador import generar_documento_markdown
    REPORTE_DISPONIBLE = True
except ImportError as e:
    console.print(f"[yellow]WARN[/yellow] Modulo de reporte no disponible ({e}).")

# NOTA DE DISENO: este modulo es 100% universal.
# NO contiene datos hardcodeados de ningun organismo (como DB_CINETICA de E. coli).
# Los datos genomicos se extraen dinamicamente de UniProt/Rhea.

ORGANISMOS = {
    "1": {"nombre": "E. coli (core)", "id": "textbook"},
    "2": {"nombre": "E. coli (completo, MG1655)", "id": "iML1515"},
    "3": {"nombre": "Saccharomyces cerevisiae (levadura)", "id": "iMM904"},
    "4": {"nombre": "Bacillus subtilis", "id": "iYO844"},
    "5": {"nombre": "Otro (ingresar ID de BiGG / BioModels / local)", "id": None}
}

COSTOS_PROTEICOS = {
    "CAFT_RXN1": 0.8, "CAFT_RXN2": 0.8,
    "EX_succ_e": 0.1, "EX_etoh_e": 0.2, "default": 0.3
}


# ==========================================
# HELPER: detectar ID de UniProt por organismo
# ==========================================
def get_organism_id(modelo_id: str) -> str:
    """Devuelve el taxon ID de UniProt segun el modelo."""
    if "iYO844" in modelo_id:
        return "224308"    # B. subtilis
    if "iCW773" in modelo_id:
        return "196627"    # C. glutamicum
    if "iMM904" in modelo_id:
        return "559292"    # S. cerevisiae
    if "textbook" in modelo_id or "iML1515" in modelo_id or "coli" in modelo_id.lower():
        return "83333"     # E. coli K-12
    return "83333"         # fallback


# ==========================================
# TRADUCTOR DE GENES (UNIPROT)
# ==========================================
def traducir_gen_uniprot(gen_id, organismo_id_uniprot):
    """Consulta UniProt. Universal: si no hay digitos, asume nombre comun."""
    if not re.search(r'\d', gen_id):
        return gen_id, "Nombre comun"

    url = (f"https://rest.uniprot.org/uniprotkb/search?"
           f"query=organism_id:{organismo_id_uniprot}+AND+gene_exact:{gen_id}"
           f"&fields=gene_names,protein_name&format=json")
    try:
        time.sleep(0.3)
        response = requests.get(url, timeout=5)
        if response.status_code == 200:
            data = response.json()
            if data.get('results'):
                result = data['results'][0]
                nombre_gen = result.get('genes', [{}])[0].get('geneName', {}).get('value', 'Desconocido')
                nombre_proteina = (result.get('proteinDescription', {})
                                    .get('recommendedName', {})
                                    .get('fullName', {})
                                    .get('value', 'Funcion desconocida'))
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
        return modelo, f"Modelo '{modelo_id}' descargado de BiGG/BioModels."
    except Exception as e:
        raise Exception(
            f"No se pudo cargar '{modelo_id}'. Detalle: {e}\n"
            "Descarga el SBML manualmente y guardalo en 'modelos/<id>.xml'."
        )


# ==========================================
# INYECCION DE RUTAS HETEROLOGAS
# ==========================================
def inyectar_ruta_heterologa(modelo, quimico):
    if quimico == "EX_caft_e":
        console.print("[cyan]INFO[/cyan] Inyectando ruta de la planta de cafe...")
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

        console.print("[green]OK[/green] Ruta de cafeina inyectada.")
        return modelo, "EX_caft_e"
    return None, None


# ==========================================
# ANALISIS TERMODINAMICO UNIVERSAL DE LA CEPA
# ==========================================
def analizar_termodinamica_cepa(modelo, quimico, genes_knockout):
    """
    Reemplaza al viejo validar_termodinamica_ruta (que era E.coli-especifico).
    Analiza termodinamicamente la reaccion objetivo + las reacciones afectadas por knockouts.
    """
    if not MODULOS_EXTERNOS_DISPONIBLES:
        console.print("[yellow]WARN[/yellow] Termodinamica universal no disponible.")
        return

    reacciones_relevantes = []
    if quimico in modelo.reactions:
        reacciones_relevantes.append(quimico)
    for gen_id in genes_knockout:
        if gen_id in modelo.genes:
            for rxn in modelo.genes.get_by_id(gen_id).reactions:
                if rxn.id not in reacciones_relevantes and not rxn.id.startswith("EX_"):
                    reacciones_relevantes.append(rxn.id)

    if not reacciones_relevantes:
        console.print("[yellow]WARN[/yellow] No hay reacciones relevantes para analizar termodinamicamente.")
        return

    try:
        analizar_termodinamica_universal(modelo, reacciones_relevantes[:10])
    except Exception as e:
        console.print(f"[yellow]WARN[/yellow] Termodinamica fallo: {e}")


# ==========================================
# EVALUACION Y ANALISIS PROFUNDO
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
    console.print(f"\n[cyan]INFO[/cyan] Analisis Profundo de la Cepa...")

    modelo_mut = modelo.copy()
    for gen_id in genes_a_mutar:
        if gen_id in modelo_mut.genes:
            for reaccion in modelo_mut.genes.get_by_id(gen_id).reactions:
                reaccion.bounds = (0, 0)

    modelo_mut.objective = quimico
    sol_final = modelo_mut.optimize()

    if sol_final.status == 'optimal':
        produccion = sol_final.objective_value
        console.print(f"   -> Produccion Bruta (FBA): [bold cyan]{produccion:.4f}[/bold cyan] mmol/gDW/h")

        flujo_rxn = sol_final.fluxes.get(quimico, 0.0)
        costo = flujo_rxn * COSTOS_PROTEICOS.get(quimico, COSTOS_PROTEICOS["default"])
        console.print(f"   -> Penalizacion Proteica: [bold red]-{costo:.2f}[/bold red]")
        console.print(f"   -> Rendimiento Ajustado: [bold green]{max(0.0, produccion - costo):.4f}[/bold green] mmol/gDW/h")

        carga = sol_final.fluxes.abs().sum()
        console.print(f"   -> Carga Metabolica Total: [bold magenta]{carga:.2f}[/bold magenta]")

        # Termodinamica UNIVERSAL (antes era E.coli-especifico)
        analizar_termodinamica_cepa(modelo_mut, quimico, genes_a_mutar)


# ==========================================
# DISENO DE CEPA POR IA (ALGORITMO GENETICO)
# ==========================================
def disenar_cepa_ia(modelo, nombre_organismo, modelo_id):
    console.print(f"\n[cyan]INFO[/cyan] Diseno de cepas por IA para {nombre_organismo}...")
    quimico = Prompt.ask("Quimico objetivo (ej. EX_succ_e, EX_caft_e)", default="EX_succ_e").strip()

    if quimico not in modelo.reactions:
        console.print(f"[red]ERROR[/red] El quimico '{quimico}' no es nativo.")
        if quimico == "EX_caft_e":
            if Prompt.ask("Inyectar ruta de la planta de cafe?", choices=["S", "N"], default="S").upper() == "S":
                modelo, quimico = inyectar_ruta_heterologa(modelo, quimico)
                if modelo is None: return
            else: return
        else:
            console.print("No hay ruta pre-programada.", style="bold red"); return

    gens = IntPrompt.ask("Generaciones", default=10)
    pop_size = IntPrompt.ask("Poblacion", default=20)
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

    console.print(f"\n[green]OK[/green] Diseno completado.")
    console.print(f"Fitness: {mejor['fitness']:.4f}\n")

    # Traduccion a nombres reales
    console.print("[cyan]INFO[/cyan] Traduciendo nombres de genes con UniProt...")
    tabla_genes = Table(title="Genes a apagar (Knockouts)")
    tabla_genes.add_column("ID (BiGG)", style="cyan")
    tabla_genes.add_column("Nombre Real", style="magenta")
    tabla_genes.add_column("Funcion Proteica", style="green")

    id_uniprot = get_organism_id(modelo_id)
    reporte_genes = []
    for g in mejor["genes"]:
        nombre_real, funcion = traducir_gen_uniprot(g, id_uniprot)
        tabla_genes.add_row(g, nombre_real, funcion[:40] + "..." if len(funcion) > 40 else funcion)
        reporte_genes.append(f"{g} ({nombre_real}): {funcion}")
    console.print(tabla_genes)

    # Reporte cientifico en Markdown (reemplaza al viejo .txt inline)
    if REPORTE_DISPONIBLE:
        try:
            generar_documento_markdown(modelo, modelo_id, quimico, mejor["genes"], mejor["fitness"])
        except Exception as e:
            console.print(f"[yellow]WARN[/yellow] Reporte fallo: {e}")
    else:
        # Fallback: guardamos un .txt basico
        if not os.path.exists("reportes"): os.makedirs("reportes")
        fecha = datetime.now().strftime("%Y%m%d_%H%M%S")
        nombre_base = f"reportes/cepa_lab_{modelo_id}_{quimico}_{fecha}"
        with open(f"{nombre_base}.txt", "w", encoding="utf-8") as f:
            f.write("REPORTE DE INGENIERIA METABOLICA\n")
            f.write(f"Organismo: {nombre_organismo}\n")
            f.write(f"Quimico: {quimico}\n")
            f.write(f"Fitness: {mejor['fitness']:.4f} mmol/gDW/h\n")
            for rg in reporte_genes:
                f.write(f"  - {rg}\n")
        console.print(f"[green]OK[/green] Reporte basico guardado: {nombre_base}.txt")

    analisis_profundo(modelo, gen_biomasa, mejor["genes"], quimico)


# ==========================================
# DISENO POR OPTKNOCK (ACOPLAMIENTO)
# ==========================================
def optknock_algoritmo(modelo, nombre_organismo, modelo_id):
    console.print(f"\n[cyan]INFO[/cyan] Diseno por OptKnock Hibrido...")

    quimico = Prompt.ask("Quimico objetivo (ej. EX_succ_e, EX_caft_e)", default="EX_succ_e").strip()

    if quimico not in modelo.reactions:
        console.print(f"[red]ERROR[/red] El quimico '{quimico}' no es nativo.")
        if quimico == "EX_caft_e":
            if Prompt.ask("Inyectar ruta de la planta de cafe?", choices=["S", "N"], default="S").upper() == "S":
                modelo, quimico = inyectar_ruta_heterologa(modelo, quimico)
                if modelo is None: return
            else: return
        else:
            console.print("No hay ruta pre-programada.", style="bold red"); return

    max_mutaciones = IntPrompt.ask("Cuántos genes apagar simultaneamente?", default=3)
    generaciones = IntPrompt.ask("Generaciones de la IA", default=15)
    pop_size = IntPrompt.ask("Tamano de poblacion", default=30)

    gen_biomasa = [r.id for r in modelo.reactions if "biomass" in r.id.lower()][0]
    modelo.reactions.get_by_id(quimico).bounds = (0, 1000.0)

    solucion_wt = modelo.optimize()
    reacciones_activas = {r for r, v in solucion_wt.fluxes.items() if abs(v) > 0.01 and not r.startswith("EX_")}
    genes_candidatos = list({g.id for r in reacciones_activas if r in modelo.reactions
                             for g in modelo.reactions.get_by_id(r).genes})

    console.print(f"   -> Usando {len(genes_candidatos)} genes activos para la evolucion.\n")

    if len(genes_candidatos) < 2:
        console.print("[red]ERROR[/red] No hay suficientes genes activos.", style="bold red")
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
        nueva = mejores.copy()
        while len(nueva) < pop_size:
            p1 = random.choice(mejores)["genes"]
            p2 = random.choice(mejores)["genes"]
            hijo = list(set(p1[:2] + p2[2:]))
            if len(hijo) < max_mutaciones:
                hijo.append(random.choice(genes_candidatos))
            nueva.append({"genes": hijo, "fitness": evaluar_acoplamiento(hijo)})
        poblacion = nueva

    mejor = sorted(poblacion, key=lambda x: x["fitness"], reverse=True)[0]

    if mejor["fitness"] > 0:
        console.print(f"\n[green]OK[/green] OptKnock exitoso.")
        console.print(f"Genes a apagar: [bold cyan]{', '.join(mejor['genes'])}[/bold cyan]")
        console.print(f"Ratio de Acoplamiento: [bold yellow]{mejor['fitness']:.4f}[/bold yellow] mmol quimico / mmol biomasa")
        console.print("La celula crece de forma optima cuando produce el quimico.")
        analisis_profundo(modelo, gen_biomasa, mejor["genes"], quimico)
    else:
        console.print("\n[red]ERROR[/red] No se encontro acoplamiento estricto.")
        console.print("Sugerencia: metabolismos flexibles requieren el Algoritmo Genetico (Opcion 2).")


# ==========================================
# ANALISIS INTEGRAL
# ==========================================
def analisis_integral(modelo, nombre_organismo, modelo_id):
    console.print(f"\n[cyan]INFO[/cyan] Analisis Integral (Diseno + Termo + Mapa + Reporte)...")

    quimico = Prompt.ask("Quimico objetivo (ej. EX_succ_e, EX_caft_e)", default="EX_succ_e")

    if quimico not in modelo.reactions:
        console.print(f"[red]ERROR[/red] El quimico '{quimico}' no es nativo.")
        if quimico == "EX_caft_e":
            if Prompt.ask("Inyectar ruta de la planta de cafe?", choices=["S", "N"], default="S").upper() == "S":
                modelo, quimico = inyectar_ruta_heterologa(modelo, quimico)
                if modelo is None: return
            else: return
        else:
            console.print("No hay ruta pre-programada.", style="bold red"); return

    gens = IntPrompt.ask("Generaciones de la IA", default=10)
    pop_size = IntPrompt.ask("Poblacion", default=20)
    mutaciones = IntPrompt.ask("Genes a apagar", default=4)

    # Fase 1: IA
    console.print("\n[cyan]INFO[/cyan] Fase 1: Evolucion Genetica...")
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
        nueva = mejores.copy()
        while len(nueva) < pop_size:
            p1 = random.choice(mejores)["genes"]
            p2 = random.choice(mejores)["genes"]
            hijo = list(set(p1[:3] + p2[3:]))
            if len(hijo) < mutaciones: hijo.append(random.choice(todos_los_genes))
            nueva.append({"genes": hijo, "fitness": evaluar_bacteria(modelo, hijo, quimico, gen_biomasa)})
        poblacion = nueva

    mejor = sorted(poblacion, key=lambda x: x["fitness"], reverse=True)[0]
    console.print(f"[green]OK[/green] Mejor Cepa: {', '.join(mejor['genes'])} (Fitness: {mejor['fitness']:.4f})")

    # Fase 2: Analisis profundo (con termo universal)
    console.print("\n[cyan]INFO[/cyan] Fase 2: Analisis Fisicoquimico...")
    analisis_profundo(modelo, gen_biomasa, mejor["genes"], quimico)

    # Fase 3: Mapa visual
    if MODULOS_EXTERNOS_DISPONIBLES:
        console.print("\n[cyan]INFO[/cyan] Fase 3: Generando Mapa Visual...")
        try:
            generar_mapa_universal(modelo_id, quimico)
        except Exception as e:
            console.print(f"[yellow]WARN[/yellow] Mapa visual fallo: {e}")
    else:
        console.print("\n[yellow]WARN[/yellow] Fase 3 omitida: modulo visual no disponible.")

    # Fase 4: Reporte cientifico
    if REPORTE_DISPONIBLE:
        console.print("\n[cyan]INFO[/cyan] Fase 4: Generando Reporte Cientifico...")
        try:
            generar_documento_markdown(modelo, modelo_id, quimico, mejor["genes"], mejor["fitness"])
        except Exception as e:
            console.print(f"[yellow]WARN[/yellow] Reporte fallo: {e}")

    console.print("\n[green]OK[/green] Analisis Integral completado.")


# ==========================================
# MENU PRINCIPAL
# ==========================================
def main():
    console.clear()
    console.print(Panel.fit(
        "[bold green]MetaboEngineer Pro V5.1[/bold green]\n"
        "Suite Universal de Biologia de Sistemas",
        style="bold green"
    ))

    console.print("\n[bold]Selecciona el organismo base:[/bold]")
    for key, val in ORGANISMOS.items():
        console.print(f"  {key}. {val['nombre']}")

    opcion = Prompt.ask("Opcion", choices=list(ORGANISMOS.keys()), default="2")

    if opcion == "9":
        modelo_id = Prompt.ask("ID del modelo (BiGG / BioModels / local)").strip()
        nombre_organismo = modelo_id
    else:
        modelo_id = ORGANISMOS[opcion]["id"]
        nombre_organismo = ORGANISMOS[opcion]["nombre"]

    with console.status(f"[bold yellow]Cargando modelo {modelo_id}...[/bold yellow]"):
        try:
            modelo, mensaje = cargar_modelo_universal(modelo_id)
            console.print(f"[green]OK[/green] {mensaje} | {len(modelo.reactions)} reacciones.\n")
        except Exception as e:
            console.print(Panel(str(e), title="ERROR de carga", border_style="red"))
            return

    while True:
        console.print("\n[bold]--- MENU PRINCIPAL (V5.1) ---[/bold]")
        console.print("1. [cyan]Analisis de Flujos Basico (FBA)[/cyan]")
        console.print("2. [green]Diseno de Cepas por IA (+ Reporte Cientifico)[/green]")
        console.print("3. [magenta]Diseno por OptKnock (Acoplamiento de Crecimiento)[/magenta]")
        console.print("4. [bold blue]Analisis Integral (Diseno + Termo + Mapa + Reporte)[/bold blue]")
        console.print("5. [bold red]Salir[/bold red]")

        op = Prompt.ask("\nSelecciona una opcion", choices=["1", "2", "3", "4", "5"])

        if op == "1":
            sol = modelo.optimize()
            if sol.status == 'optimal':
                console.print(f"[green]OK[/green] Crecimiento: [bold green]{sol.objective_value:.4f}[/bold green] h-1")
        elif op == "2":
            disenar_cepa_ia(modelo, nombre_organismo, modelo_id)
        elif op == "3":
            optknock_algoritmo(modelo, nombre_organismo, modelo_id)
        elif op == "4":
            analisis_integral(modelo, nombre_organismo, modelo_id)
        elif op == "5":
            console.print("Saliendo. Hasta luego.")
            break


if __name__ == "__main__":
    main()
    