import warnings
warnings.filterwarnings("ignore")

import cobra
from cobra import Reaction, Metabolite
from cobra.io import load_model, read_sbml_model, write_sbml_model
import random
import os
import sys
import requests
import time
import re
import json
import hashlib
import logging
import pandas as pd
from requests.exceptions import RequestException
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
    console.print(f"[yellow]WARN[/yellow] Modulos externos no disponibles ({e}).")

REPORTE_DISPONIBLE = False
try:
    from reporte_investigador import generar_documento_markdown
    REPORTE_DISPONIBLE = True
except ImportError:
    pass

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

_RUTA_HISTORIAL_CSV = "reportes/historial.csv"

# ==========================================
# LOGGING A ARCHIVO
# ==========================================
def setup_logger():
    if not os.path.exists("logs"):
        os.makedirs("logs")
    log_file = f"logs/run_{datetime.now():%Y%m%d_%H%M%S}.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[logging.FileHandler(log_file, encoding="utf-8")]
    )
    logging.info("=" * 60)
    logging.info("MetaboEngineer Pro V5.9 - Session started")
    logging.info("=" * 60)
    return log_file

# ==========================================
# SEED REPRODUCIBLE
# ==========================================
def preguntar_seed():
    try:
        entrada = Prompt.ask("Semilla aleatoria (Enter = aleatoria)", default="").strip()
        if entrada == "":
            seed = random.randint(0, 999999)
            console.print(f"[cyan]INFO[/cyan] Semilla aleatoria: [bold]{seed}[/bold]")
            return seed
        seed = int(entrada)
        console.print(f"[cyan]INFO[/cyan] Semilla fija: [bold]{seed}[/bold]")
        return seed
    except (ValueError, KeyboardInterrupt):
        seed = random.randint(0, 999999)
        return seed

# ==========================================
# HELPER: detectar ID de UniProt por organismo
# ==========================================
def get_organism_id(modelo_id: str) -> str:
    if not modelo_id:
        return "83333"
    if "iYO844" in modelo_id: return "224308"
    if "iCW773" in modelo_id: return "196627"
    if "iMM904" in modelo_id: return "559292"
    if "iEK1008" in modelo_id: return "1773"
    return "83333"

# ==========================================
# CACHE UNIPROT (offline)
# ==========================================
_RUTA_CACHE_UNIPROT = "cache_uniprot.json"

def _cargar_cache_uniprot():
    if not os.path.exists(_RUTA_CACHE_UNIPROT):
        return {}
    try:
        with open(_RUTA_CACHE_UNIPROT, "r", encoding="utf-8") as f:
            return json.load(f)
    except Exception:
        return {}

def _guardar_cache_uniprot(cache):
    tmp = _RUTA_CACHE_UNIPROT + ".tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(cache, f, indent=2, ensure_ascii=False)
    os.replace(tmp, _RUTA_CACHE_UNIPROT)

def limpiar_cache_uniprot():
    if os.path.exists(_RUTA_CACHE_UNIPROT):
        os.remove(_RUTA_CACHE_UNIPROT)
        console.print(f"[green]OK[/green] Cache UniProt eliminado.")

# ==========================================
# TRADUCTOR DE GENES (UNIPROT)
# ==========================================
def traducir_gen_uniprot(gen_id, organismo_id_uniprot):
    if not re.search(r'\d', gen_id):
        return gen_id, "Nombre comun"

    cache = _cargar_cache_uniprot()
    clave = f"{organismo_id_uniprot}:{gen_id}"
    if clave in cache:
        return tuple(cache[clave])

    # Búsqueda global: busca el ID en cualquier campo (soluciona locus_tag)
    url = f'https://rest.uniprot.org/uniprotkb/search?query=organism_id:{organismo_id_uniprot}+AND+"{gen_id}"&fields=gene_names,protein_name&format=json'
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
                cache[clave] = [nombre_gen, nombre_proteina]
                _guardar_cache_uniprot(cache)
                return nombre_gen, nombre_proteina
    except Exception:
        pass
    return gen_id, "No encontrado en UniProt"

# ==========================================
# EXPORTACION A JSON/CSV
# ==========================================
def _exportar_resultado_json(modelo_id, quimico, metodo, genes, fitness, seed, parametros_extra=None):
    if not os.path.exists("reportes"): os.makedirs("reportes")
    fecha = datetime.now()
    timestamp = fecha.strftime("%Y%m%d_%H%M%S")

    resultado = {"fecha": fecha.isoformat(), "modelo": modelo_id, "quimico": quimico, "metodo": metodo, "fitness": round(fitness, 6), "genes_knockout": list(genes), "num_genes": len(genes), "seed": seed, "parametros": parametros_extra or {}}
    ruta_json = f"reportes/{metodo}_{modelo_id}_{quimico}_{timestamp}.json"
    with open(ruta_json, "w", encoding="utf-8") as f:
        json.dump(resultado, f, indent=2, ensure_ascii=False)

    fila = {"fecha": fecha.strftime("%Y-%m-%d %H:%M:%S"), "modelo": modelo_id, "quimico": quimico, "metodo": metodo, "fitness": round(fitness, 6), "num_genes": len(genes), "genes": ";".join(genes), "seed": seed}
    df = pd.DataFrame([fila])
    df.to_csv(_RUTA_HISTORIAL_CSV, mode="a", header=not os.path.exists(_RUTA_HISTORIAL_CSV), index=False)

    console.print(f"[green]OK[/green] JSON: {ruta_json}")
    console.print(f"[green]OK[/green] CSV: {_RUTA_HISTORIAL_CSV} (append)")

def mostrar_historial():
    if not os.path.exists(_RUTA_HISTORIAL_CSV): return
    df = pd.read_csv(_RUTA_HISTORIAL_CSV)
    if df.empty: return
    console.print(df.tail(20).to_string())

# ==========================================
# DESCARGA ROBUSTA E INDEXACION
# ==========================================
_HEADERS_NAVEGADOR = {"User-Agent": "Mozilla/5.0", "Accept": "application/xml,text/xml,*/*;q=0.9"}
def _bigg_online(): return True
def _descargar_xml(url, ruta_local, max_intentos=2, timeout=25): return False, "error"
def _urls_candidatas(m_id): return []

_RUTA_INDICE = "modelos/_index.json"
def _sha256_archivo(ruta):
    h = hashlib.sha256()
    with open(ruta, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""): h.update(chunk)
    return h.hexdigest()

def _cargar_indice():
    if not os.path.exists(_RUTA_INDICE): return {"version": 1, "modelos": {}}
    try:
        with open(_RUTA_INDICE, "r", encoding="utf-8") as f: return json.load(f)
    except: return {"version": 1, "modelos": {}}

def _guardar_indice(indice):
    with open(_RUTA_INDICE + ".tmp", "w", encoding="utf-8") as f: json.dump(indice, f, indent=2)
    os.replace(_RUTA_INDICE + ".tmp", _RUTA_INDICE)

def _registrar_modelo(modelo, modelo_id, ruta, fuente, sha256=None):
    indice = _cargar_indice()
    if not sha256: sha256 = _sha256_archivo(ruta)
    indice["modelos"][modelo_id] = {"archivo": os.path.basename(ruta), "sha256": sha256, "fecha": datetime.now().strftime("%Y-%m-%dT%H:%M:%S"), "reacciones": len(modelo.reactions), "genes": len(modelo.genes)}
    _guardar_indice(indice)

def listar_modelos_cacheados(): return _cargar_indice().get("modelos", {})

def _sincronizar_indice():
    if not os.path.exists("modelos"): return
    indice = _cargar_indice()
    for archivo in os.listdir("modelos"):
        if archivo.endswith(".xml") and archivo[:-4] not in indice["modelos"]:
            try: _registrar_modelo(read_sbml_model(os.path.join("modelos", archivo)), archivo[:-4], os.path.join("modelos", archivo), "manual")
            except: pass

def cargar_modelo_universal(modelo_id):
    if not modelo_id or not str(modelo_id).strip(): raise Exception("ID de modelo vacio.")
    modelo_id = str(modelo_id).strip()
    if not os.path.exists("modelos"): os.makedirs("modelos")
    ruta_local = f"modelos/{modelo_id}.xml"
    if os.path.exists(ruta_local):
        modelo = read_sbml_model(ruta_local)
        _registrar_modelo(modelo, modelo_id, ruta_local, "local")
        return modelo, f"Modelo '{modelo_id}' cargado desde disco local."
    try:
        modelo = load_model(modelo_id)
        try:
            write_sbml_model(modelo, ruta_local)
            _registrar_modelo(modelo, modelo_id, ruta_local, "cobrapy")
        except: pass
        return modelo, f"Modelo '{modelo_id}' descargado via cobrapy."
    except Exception as e:
        raise Exception(f"No se pudo cargar '{modelo_id}': {e}")

# ==========================================
# INYECCION Y ANALISIS
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
            if met.id not in modelo.metabolites: modelo.add_metabolites([met])
        reacciones = [Reaction('CAFT_RXN1'), Reaction('CAFT_RXN2'), Reaction('EX_caft_e'), Reaction('CAFT_TRANSPORT'), Reaction('EX_xan_e'), Reaction('XAN_TRANSPORT')]
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

def analizar_termodinamica_cepa(modelo, quimico, genes_knockout):
    if not MODULOS_EXTERNOS_DISPONIBLES: return
    reacciones = [quimico] if quimico in modelo.reactions else []
    for gen_id in genes_knockout:
        if gen_id in modelo.genes:
            for rxn in modelo.genes.get_by_id(gen_id).reactions:
                if rxn.id not in reacciones and not rxn.id.startswith("EX_"): reacciones.append(rxn.id)
    try: analizar_termodinamica_universal(modelo, reacciones[:10])
    except Exception as e: console.print(f"[yellow]WARN[/yellow] Termodinamica fallo: {e}")

def mostrar_cinetica_local(modelo, quimico, genes_knockout):
    if not os.path.exists("cinetica_local.json"): return
    try:
        with open("cinetica_local.json", "r", encoding="utf-8") as f: db = json.load(f)
    except: return

    reacciones = [quimico]
    for gen_id in genes_knockout:
        if gen_id in modelo.genes:
            for rxn in modelo.genes.get_by_id(gen_id).reactions:
                if rxn.id not in reacciones: reacciones.append(rxn.id)

    tabla = Table(title="Parametros Cineticos (Datos Locales)")
    tabla.add_column("Reaccion", style="cyan"); tabla.add_column("Enzima", style="magenta"); tabla.add_column("kcat (1/s)"); tabla.add_column("Km (mM)")
    
    encontrados = 0
    for r_id in reacciones:
        if r_id in db:
            d = db[r_id]
            tabla.add_row(r_id, d.get("enzima", "-"), str(d.get("kcat", "-")), str(d.get("Km", "-")))
            encontrados += 1
    if encontrados > 0: console.print(tabla)

def evaluar_bacteria(modelo_base, genes_apagados, gen_objetivo, gen_biomasa):
    modelo = modelo_base.copy()
    for gen_id in genes_apagados:
        if gen_id in modelo.genes:
            for reaccion in modelo.genes.get_by_id(gen_id).reactions: reaccion.bounds = (0, 0)
    try:
        modelo.reactions.get_by_id(gen_biomasa).bounds = (0.1, 1000.0)
        modelo.objective = gen_objetivo
        solucion = modelo.optimize()
        if solucion.status == 'optimal':
            produccion = solucion.objective_value
            costo = solucion.fluxes.get(gen_objetivo, 0.0) * COSTOS_PROTEICOS.get(gen_objetivo, COSTOS_PROTEICOS["default"])
            return max(0.0, produccion - costo)
    except: pass
    return 0.0

def analisis_profundo(modelo, gen_biomasa, genes_a_mutar, quimico):
    console.print(f"\n[cyan]INFO[/cyan] Analisis Profundo de la Cepa...")
    modelo_mut = modelo.copy()
    for gen_id in genes_a_mutar:
        if gen_id in modelo_mut.genes:
            for reaccion in modelo_mut.genes.get_by_id(gen_id).reactions: reaccion.bounds = (0, 0)
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
        analizar_termodinamica_cepa(modelo_mut, quimico, genes_a_mutar)
        mostrar_cinetica_local(modelo, quimico, genes_a_mutar)

# ==========================================
# ALGORITMOS DE DISEÑO (CON FALLBACK LOCAL)
# ==========================================
def disenar_cepa_ia(modelo, nombre_organismo, modelo_id, seed=None):
    if seed is not None: random.seed(seed)
    console.print(f"\n[cyan]INFO[/cyan] Diseno de cepas por IA para {nombre_organismo}...")
    quimico = Prompt.ask("Quimico objetivo (ej. EX_succ_e, EX_caft_e)", default="EX_succ_e").strip()

    if quimico not in modelo.reactions:
        if quimico == "EX_caft_e" and Prompt.ask("Inyectar ruta?", choices=["S", "N"], default="S").upper() == "S":
            modelo, quimico = inyectar_ruta_heterologa(modelo, quimico)
            if modelo is None: return
        else: return

    gens = IntPrompt.ask("Generaciones", default=10); pop_size = IntPrompt.ask("Poblacion", default=20); mutaciones = IntPrompt.ask("Genes a apagar", default=4)
    gen_biomasa = [r.id for r in modelo.reactions if "biomass" in r.id.lower()][0]
    modelo.reactions.get_by_id(quimico).bounds = (0, 1000.0)
    todos_los_genes = [g.id for g in modelo.genes]

    poblacion = []
    for _ in range(pop_size):
        gm = random.sample(todos_los_genes, mutaciones)
        poblacion.append({"genes": gm, "fitness": evaluar_bacteria(modelo, gm, quimico, gen_biomasa)})

    for gen in track(range(gens), description="Evolucionando..."):
        poblacion.sort(key=lambda x: x["fitness"], reverse=True)
        mejores = poblacion[:5]
        nueva = mejores.copy()
        while len(nueva) < pop_size:
            p1, p2 = random.choice(mejores)["genes"], random.choice(mejores)["genes"]
            hijo = list(set(p1[:3] + p2[3:]))
            if len(hijo) < mutaciones: hijo.append(random.choice(todos_los_genes))
            nueva.append({"genes": hijo, "fitness": evaluar_bacteria(modelo, hijo, quimico, gen_biomasa)})
        poblacion = nueva

    mejor = sorted(poblacion, key=lambda x: x["fitness"], reverse=True)[0]
    console.print(f"\n[green]OK[/green] Diseno completado. Fitness: {mejor['fitness']:.4f}\n")
    
    id_uniprot = get_organism_id(modelo_id)
    console.print("[cyan]INFO[/cyan] Traduciendo nombres de genes (Modelo local -> UniProt)...")
    tabla = Table(title="Genes a apagar (Knockouts)")
    tabla.add_column("ID (BiGG)", style="cyan"); tabla.add_column("Nombre Real", style="magenta"); tabla.add_column("Funcion Proteica", style="green")
    
    for g in mejor["genes"]:
        # SOLUCIÓN: Primero intenta sacar el nombre del modelo local de COBRA (infalible)
        nombre_local = modelo.genes.get_by_id(g).name
        if nombre_local and nombre_local != g:
            nombre_real, funcion = nombre_local, "Nombre extraido del modelo local"
        else:
            # Si el modelo no tiene el nombre, busca en UniProt
            nombre_real, funcion = traducir_gen_uniprot(g, id_uniprot)
        tabla.add_row(g, nombre_real, funcion[:50] + "..." if len(funcion) > 50 else funcion)
    console.print(tabla)

    if REPORTE_DISPONIBLE:
        try: generar_documento_markdown(modelo, modelo_id, quimico, mejor["genes"], mejor["fitness"])
        except Exception as e: console.print(f"[yellow]WARN[/yellow] Reporte fallo: {e}")

    try: _exportar_resultado_json(modelo_id, quimico, "ia", mejor["genes"], mejor["fitness"], seed, {"gens": gens})
    except: pass
    analisis_profundo(modelo, gen_biomasa, mejor["genes"], quimico)

def optknock_algoritmo(modelo, nombre_organismo, modelo_id, seed=None):
    if seed is not None: random.seed(seed)
    quimico = Prompt.ask("Quimico objetivo", default="EX_succ_e").strip()
    if quimico not in modelo.reactions: return

    max_mutaciones = IntPrompt.ask("Genes a apagar?", default=3); generaciones = IntPrompt.ask("Generaciones?", default=15); pop_size = IntPrompt.ask("Poblacion?", default=30)
    gen_biomasa = [r.id for r in modelo.reactions if "biomass" in r.id.lower()][0]
    modelo.reactions.get_by_id(quimico).bounds = (0, 1000.0)
    reacciones_activas = {r for r, v in modelo.optimize().fluxes.items() if abs(v) > 0.01 and not r.startswith("EX_")}
    genes_candidatos = list({g.id for r in reacciones_activas if r in modelo.reactions for g in modelo.reactions.get_by_id(r).genes})
    
    if len(genes_candidatos) < 2: return

    def evaluar(genes):
        m = modelo.copy()
        for gen in genes:
            if gen in m.genes: m.genes.get_by_id(gen).knock_out()
        m.objective = gen_biomasa
        sol = m.optimize()
        return (sol.fluxes.get(quimico, 0.0) / sol.objective_value) if (sol.status == 'optimal' and sol.objective_value > 0.01) else 0.0

    poblacion = []
    for _ in range(pop_size):
        gm = random.sample(genes_candidatos, min(max_mutaciones, len(genes_candidatos)))
        poblacion.append({"genes": gm, "fitness": evaluar(gm)})
        
    for gen in track(range(generaciones), description="Evolucionando..."):
        poblacion.sort(key=lambda x: x["fitness"], reverse=True)
        mejores = poblacion[:5]
        nueva = mejores.copy()
        while len(nueva) < pop_size:
            p1, p2 = random.choice(mejores)["genes"], random.choice(mejores)["genes"]
            hijo = list(set(p1[:2] + p2[2:]))
            if len(hijo) < max_mutaciones: hijo.append(random.choice(genes_candidatos))
            nueva.append({"genes": hijo, "fitness": evaluar(hijo)})
        poblacion = nueva

    mejor = sorted(poblacion, key=lambda x: x["fitness"], reverse=True)[0]
    if mejor["fitness"] > 0:
        console.print(f"\n[green]OK[/green] OptKnock exitoso. Ratio: {mejor['fitness']:.4f}")
        try: _exportar_resultado_json(modelo_id, quimico, "optknock", mejor["genes"], mejor["fitness"], seed)
        except: pass
        analisis_profundo(modelo, gen_biomasa, mejor["genes"], quimico)
    else: console.print("\n[red]ERROR[/red] No se encontro acoplamiento.")

def analisis_integral(modelo, nombre_organismo, modelo_id, seed=None):
    if seed is not None: random.seed(seed)
    quimico = Prompt.ask("Quimico objetivo", default="EX_succ_e")
    if quimico not in modelo.reactions:
        if quimico == "EX_caft_e" and Prompt.ask("Inyectar ruta?", choices=["S", "N"], default="S").upper() == "S":
            modelo, quimico = inyectar_ruta_heterologa(modelo, quimico)
            if modelo is None: return
        else: return

    gens = IntPrompt.ask("Generaciones de la IA", default=10); pop_size = IntPrompt.ask("Poblacion", default=20); mutaciones = IntPrompt.ask("Genes a apagar", default=4)
    gen_biomasa = [r.id for r in modelo.reactions if "biomass" in r.id.lower()][0]
    modelo.reactions.get_by_id(quimico).bounds = (0, 1000.0)
    todos_los_genes = [g.id for g in modelo.genes]

    poblacion = []
    for _ in range(pop_size):
        gm = random.sample(todos_los_genes, mutaciones)
        poblacion.append({"genes": gm, "fitness": evaluar_bacteria(modelo, gm, quimico, gen_biomasa)})
    for gen in track(range(gens), description="Evolucionando..."):
        poblacion.sort(key=lambda x: x["fitness"], reverse=True)
        mejores = poblacion[:5]
        nueva = mejores.copy()
        while len(nueva) < pop_size:
            p1, p2 = random.choice(mejores)["genes"], random.choice(mejores)["genes"]
            hijo = list(set(p1[:3] + p2[3:]))
            if len(hijo) < mutaciones: hijo.append(random.choice(todos_los_genes))
            nueva.append({"genes": hijo, "fitness": evaluar_bacteria(modelo, hijo, quimico, gen_biomasa)})
        poblacion = nueva

    mejor = sorted(poblacion, key=lambda x: x["fitness"], reverse=True)[0]
    console.print(f"[green]OK[/green] Mejor Cepa: {', '.join(mejor['genes'])} (Fitness: {mejor['fitness']:.4f})")
    analisis_profundo(modelo, gen_biomasa, mejor["genes"], quimico)

    if MODULOS_EXTERNOS_DISPONIBLES:
        try: generar_mapa_universal(modelo_id, quimico)
        except Exception as e: console.print(f"[yellow]WARN[/yellow] Mapa visual fallo: {e}")
    if REPORTE_DISPONIBLE:
        try: generar_documento_markdown(modelo, modelo_id, quimico, mejor["genes"], mejor["fitness"])
        except Exception as e: console.print(f"[yellow]WARN[/yellow] Reporte fallo: {e}")
    try: _exportar_resultado_json(modelo_id, quimico, "integral", mejor["genes"], mejor["fitness"], seed)
    except: pass

# ==========================================
# MENU PRINCIPAL
# ==========================================
def main():
    console.clear()
    console.print(Panel.fit("[bold green]MetaboEngineer Pro V5.9[/bold green]\nSuite Universal de Biologia de Sistemas", style="bold green"))
    setup_logger(); _sincronizar_indice()

    console.print("\n[bold]Selecciona el organismo base:[/bold]")
    for key, val in ORGANISMOS.items(): console.print(f"  {key}. {val['nombre']}")

    cacheados = listar_modelos_cacheados()
    claves_cache = []
    if cacheados:
        console.print("\n[bold cyan]Modelos en cache:[/bold cyan]")
        for i, (mid, meta) in enumerate(cacheados.items(), start=1):
            console.print(f"  \\[c{i}] {mid} — {meta['reacciones']} reacciones")
            claves_cache.append(f"c{i}")
        claves_cache.append("c0"); console.print("  \\[c0] Seleccion normal")

    try: opcion = Prompt.ask("Opcion", choices=list(ORGANISMOS.keys()) + claves_cache, default="2")
    except KeyboardInterrupt: return

    if opcion.startswith("c") and opcion != "c0":
        modelo_id = list(cacheados.keys())[int(opcion[1:]) - 1]
        nombre_organismo = modelo_id
    else:
        if opcion == "c0": opcion = Prompt.ask("Opcion", choices=list(ORGANISMOS.keys()), default="2")
        OTRO_KEY = next((k for k, v in ORGANISMOS.items() if v["id"] is None), None)
        if opcion == OTRO_KEY:
            modelo_id = Prompt.ask("ID del modelo").strip()
            if not modelo_id: return
            nombre_organismo = modelo_id
        else:
            modelo_id = ORGANISMOS[opcion]["id"]
            nombre_organismo = ORGANISMOS[opcion]["nombre"]

    seed = preguntar_seed()

    with console.status(f"[bold yellow]Cargando modelo...[/bold yellow]"):
        try:
            modelo, msg = cargar_modelo_universal(modelo_id)
            console.print(f"[green]OK[/green] {msg} | {len(modelo.reactions)} reacciones.\n")
        except Exception as e:
            console.print(Panel(str(e), title="ERROR", border_style="red")); return

    while True:
        console.print("\n[bold]--- MENU (V5.9) ---[/bold]")
        console.print("1. FBA | 2. IA | 3. OptKnock | 4. Integral | 5. Historial | 6. Limpiar Cache | 7. Salir")
        try: op = Prompt.ask("\nSelecciona una opcion", choices=["1", "2", "3", "4", "5", "6", "7"])
        except KeyboardInterrupt: break

        if op == "1":
            sol = modelo.optimize()
            if sol.status == 'optimal': console.print(f"[green]OK[/green] Crecimiento: {sol.objective_value:.4f} h-1")
        elif op == "2": disenar_cepa_ia(modelo, nombre_organismo, modelo_id, seed)
        elif op == "3": optknock_algoritmo(modelo, nombre_organismo, modelo_id, seed)
        elif op == "4": analisis_integral(modelo, nombre_organismo, modelo_id, seed)
        elif op == "5": mostrar_historial()
        elif op == "6": limpiar_cache_uniprot()
        elif op == "7": break

if __name__ == "__main__":
    main()