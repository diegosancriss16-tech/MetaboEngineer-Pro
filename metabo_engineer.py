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
    console.print(f"[yellow]WARN[/yellow] Modulos externos no disponibles ({e}). Modo reducido.")

REPORTE_DISPONIBLE = False
try:
    from reporte_investigador import generar_documento_markdown
    REPORTE_DISPONIBLE = True
except ImportError as e:
    console.print(f"[yellow]WARN[/yellow] Modulo de reporte no disponible ({e}).")


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
    logging.info("MetaboEngineer Pro V5.5 - Session started")
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
            console.print(f"[cyan]INFO[/cyan] Semilla aleatoria generada: [bold]{seed}[/bold]")
            console.print(f"[dim]Para reproducir este resultado, usa: --seed {seed}[/dim]")
            return seed
        seed = int(entrada)
        console.print(f"[cyan]INFO[/cyan] Semilla fija: [bold]{seed}[/bold]")
        return seed
    except (ValueError, KeyboardInterrupt):
        seed = random.randint(0, 999999)
        console.print(f"[yellow]WARN[/yellow] Semilla invalida, usando: {seed}")
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
        logging.info("Cache UniProt eliminado manualmente")
    else:
        console.print("[yellow]WARN[/yellow] No hay cache que borrar.")


# ==========================================
# TRADUCTOR DE GENES (UNIPROT) - CON CACHE
# ==========================================
def traducir_gen_uniprot(gen_id, organismo_id_uniprot):
    if not re.search(r'\d', gen_id):
        return gen_id, "Nombre comun"

    cache = _cargar_cache_uniprot()
    clave = f"{organismo_id_uniprot}:{gen_id}"
    if clave in cache:
        logging.debug(f"UniProt cache HIT: {clave}")
        return tuple(cache[clave])

    logging.debug(f"UniProt cache MISS: {clave} -> consultando API")
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
                cache[clave] = [nombre_gen, nombre_proteina]
                _guardar_cache_uniprot(cache)
                return nombre_gen, nombre_proteina
    except Exception as e:
        logging.warning(f"UniProt error para {gen_id}: {e}")
    return gen_id, "No encontrado en UniProt"


# ==========================================
# EXPORTACION A JSON/CSV
# ==========================================
def _exportar_resultado_json(modelo_id, quimico, metodo, genes, fitness, seed, parametros_extra=None):
    """Exporta una corrida como JSON individual + append a CSV historial."""
    if not os.path.exists("reportes"):
        os.makedirs("reportes")

    fecha = datetime.now()
    timestamp = fecha.strftime("%Y%m%d_%H%M%S")

    # JSON individual
    resultado = {
        "fecha": fecha.isoformat(),
        "modelo": modelo_id,
        "quimico": quimico,
        "metodo": metodo,
        "fitness": round(fitness, 6),
        "genes_knockout": list(genes),
        "num_genes": len(genes),
        "seed": seed,
        "parametros": parametros_extra or {},
    }

    ruta_json = f"reportes/{metodo}_{modelo_id}_{quimico}_{timestamp}.json"
    with open(ruta_json, "w", encoding="utf-8") as f:
        json.dump(resultado, f, indent=2, ensure_ascii=False)

    # Append a CSV historial
    fila = {
        "fecha": fecha.strftime("%Y-%m-%d %H:%M:%S"),
        "modelo": modelo_id,
        "quimico": quimico,
        "metodo": metodo,
        "fitness": round(fitness, 6),
        "num_genes": len(genes),
        "genes": ";".join(genes),
        "seed": seed,
    }

    existe = os.path.exists(_RUTA_HISTORIAL_CSV)
    df = pd.DataFrame([fila])
    df.to_csv(_RUTA_HISTORIAL_CSV, mode="a", header=not existe, index=False)

    console.print(f"[green]OK[/green] JSON: {ruta_json}")
    console.print(f"[green]OK[/green] CSV: {_RUTA_HISTORIAL_CSV} (append)")
    logging.info(f"Resultado exportado: {ruta_json} + {_RUTA_HISTORIAL_CSV}")

    return ruta_json


def mostrar_historial():
    """Muestra el historial acumulativo de corridas como tabla."""
    if not os.path.exists(_RUTA_HISTORIAL_CSV):
        console.print("[yellow]WARN[/yellow] No hay historial aun. Ejecuta un analisis primero.")
        return

    try:
        df = pd.read_csv(_RUTA_HISTORIAL_CSV)
    except Exception as e:
        console.print(f"[red]ERROR[/red] No se pudo leer el historial: {e}")
        return

    if df.empty:
        console.print("[yellow]WARN[/yellow] El historial esta vacio.")
        return

    tabla = Table(title=f"Historial de corridas ({len(df)} total)")
    tabla.add_column("Fecha", style="dim")
    tabla.add_column("Modelo", style="cyan")
    tabla.add_column("Quimico", style="magenta")
    tabla.add_column("Metodo", style="yellow")
    tabla.add_column("Fitness", style="green", justify="right")
    tabla.add_column("Genes", justify="right")
    tabla.add_column("Seed", justify="right")

    for _, row in df.tail(20).iterrows():
        tabla.add_row(
            str(row["fecha"]),
            str(row["modelo"]),
            str(row["quimico"]),
            str(row["metodo"]),
            f"{row['fitness']:.4f}",
            str(int(row["num_genes"])),
            str(row["seed"]),
        )

    console.print(tabla)

    # Resumen estadistico
    if len(df) >= 2:
        console.print(f"\n[bold cyan]Resumen estadistico ({len(df)} corridas):[/bold cyan]")
        stats = df.groupby("metodo")["fitness"].agg(["count", "mean", "max", "min"]).round(4)
        tabla_stats = Table(title="Fitness por metodo")
        tabla_stats.add_column("Metodo", style="yellow")
        tabla_stats.add_column("Corridas", justify="right")
        tabla_stats.add_column("Media", justify="right", style="cyan")
        tabla_stats.add_column("Max", justify="right", style="green")
        tabla_stats.add_column("Min", justify="right", style="red")
        for metodo, row in stats.iterrows():
            tabla_stats.add_row(
                str(metodo),
                str(int(row["count"])),
                f"{row['mean']:.4f}",
                f"{row['max']:.4f}",
                f"{row['min']:.4f}",
            )
        console.print(tabla_stats)

    console.print(f"\n[dim]Mostrando ultimas {min(20, len(df))} de {len(df)}. Archivo: {_RUTA_HISTORIAL_CSV}[/dim]")


# ==========================================
# DESCARGA ROBUSTA
# ==========================================
_HEADERS_NAVEGADOR = {
    "User-Agent": (
        "Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
        "AppleWebKit/537.36 (KHTML, like Gecko) "
        "Chrome/120.0.0.0 Safari/537.36"
    ),
    "Accept": "application/xml,text/xml,*/*;q=0.9",
    "Accept-Language": "en-US,en;q=0.9",
    "Connection": "keep-alive",
}


def _bigg_online(timeout=5) -> bool:
    try:
        r = requests.head("http://bigg.ucsd.edu/", timeout=timeout, allow_redirects=True)
        return r.status_code < 500
    except Exception:
        return False


def _descargar_xml(url, ruta_local, max_intentos=2, timeout=25):
    ultimo_error = None
    for intento in range(1, max_intentos + 1):
        try:
            r = requests.get(url, headers=_HEADERS_NAVEGADOR, timeout=timeout, allow_redirects=True)
            if r.status_code == 200 and len(r.content) > 1000:
                with open(ruta_local, "wb") as f:
                    f.write(r.content)
                return True, f"OK ({len(r.content)} bytes)"
            ultimo_error = f"HTTP {r.status_code} ({len(r.content)} bytes)"
        except RequestException as e:
            ultimo_error = str(e)[:120]
        if intento < max_intentos:
            time.sleep(2 * intento)
    return False, ultimo_error or "desconocido"


def _urls_candidatas(modelo_id):
    return [
        f"http://bigg.ucsd.edu/static/models/{modelo_id}.xml",
        f"https://bigg.ucsd.edu/static/models/{modelo_id}.xml",
        f"https://raw.githubusercontent.com/opencobra/cobratestdata/master/models/{modelo_id}.xml",
    ]


# ==========================================
# INDEXACION Y DEDUPLICACION DE MODELOS
# ==========================================
_RUTA_INDICE = "modelos/_index.json"


def _sha256_archivo(ruta):
    h = hashlib.sha256()
    with open(ruta, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def _cargar_indice():
    if not os.path.exists(_RUTA_INDICE):
        return {"version": 1, "modelos": {}}
    try:
        with open(_RUTA_INDICE, "r", encoding="utf-8") as f:
            data = json.load(f)
        if "modelos" not in data:
            data["modelos"] = {}
        return data
    except Exception:
        console.print("[yellow]WARN[/yellow] Indice corrupto, recreando...")
        return {"version": 1, "modelos": {}}


def _guardar_indice(indice):
    tmp = _RUTA_INDICE + ".tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(indice, f, indent=2, ensure_ascii=False)
    os.replace(tmp, _RUTA_INDICE)


def _detectar_duplicado(indice, sha_nuevo, modelo_id_actual):
    for otro_id, meta in indice["modelos"].items():
        if otro_id == modelo_id_actual:
            continue
        if meta.get("sha256") == sha_nuevo:
            return otro_id
    return None


def _registrar_modelo(modelo, modelo_id, ruta_archivo, fuente, sha256=None):
    indice = _cargar_indice()
    if sha256 is None:
        sha256 = _sha256_archivo(ruta_archivo)
    alias_de = _detectar_duplicado(indice, sha256, modelo_id)
    indice["modelos"][modelo_id] = {
        "archivo": os.path.basename(ruta_archivo),
        "sha256": sha256,
        "fecha": datetime.now().strftime("%Y-%m-%dT%H:%M:%S"),
        "reacciones": len(modelo.reactions),
        "metabolitos": len(modelo.metabolites),
        "genes": len(modelo.genes),
        "fuente": fuente,
        "tamaño_bytes": os.path.getsize(ruta_archivo),
        "alias_de": alias_de,
    }
    _guardar_indice(indice)
    if alias_de:
        console.print(f"[yellow]INFO[/yellow] '{modelo_id}' es identico a '{alias_de}' (mismo hash).")
    else:
        console.print(f"[green]OK[/green] Modelo '{modelo_id}' registrado en indice.")


def listar_modelos_cacheados():
    indice = _cargar_indice()
    return indice.get("modelos", {})


def _sincronizar_indice():
    if not os.path.exists("modelos"):
        return
    indice = _cargar_indice()
    ids_registrados = set(indice["modelos"].keys())
    for archivo in os.listdir("modelos"):
        if not archivo.endswith(".xml"):
            continue
        modelo_id = archivo[:-4]
        if modelo_id in ids_registrados:
            continue
        ruta = os.path.join("modelos", archivo)
        try:
            modelo_tmp = read_sbml_model(ruta)
            _registrar_modelo(modelo_tmp, modelo_id, ruta, fuente="manual")
        except Exception as e:
            console.print(f"[yellow]WARN[/yellow] No se pudo indexar '{archivo}': {e}")


# ==========================================
# CARGA UNIVERSAL
# ==========================================
def cargar_modelo_universal(modelo_id):
    if not modelo_id or not isinstance(modelo_id, str) or not modelo_id.strip():
        raise Exception(
            "ID de modelo vacio o invalido. Ingresa un ID valido "
            "(ej. iML1515, iYO844, iCW773) o coloca el .xml en 'modelos/'."
        )
    modelo_id = modelo_id.strip()

    if not os.path.exists("modelos"):
        os.makedirs("modelos")

    ruta_local = f"modelos/{modelo_id}.xml"
    indice = _cargar_indice()

    if os.path.exists(ruta_local):
        if modelo_id in indice["modelos"]:
            sha_esperado = indice["modelos"][modelo_id].get("sha256")
            sha_real = _sha256_archivo(ruta_local)
            if sha_esperado and sha_real != sha_esperado:
                console.print(f"[yellow]WARN[/yellow] '{modelo_id}.xml' local corrupto. Re-descargando.")
                os.remove(ruta_local)
                del indice["modelos"][modelo_id]
                _guardar_indice(indice)
            else:
                modelo = read_sbml_model(ruta_local)
                if indice["modelos"][modelo_id].get("alias_de"):
                    alias = indice["modelos"][modelo_id]["alias_de"]
                    console.print(f"[yellow]INFO[/yellow] '{modelo_id}' es alias de '{alias}'.")
                logging.info(f"Modelo cargado desde disco: {modelo_id}")
                return modelo, f"Modelo '{modelo_id}' cargado desde disco local (indexado)."
        else:
            modelo = read_sbml_model(ruta_local)
            _registrar_modelo(modelo, modelo_id, ruta_local, fuente="local_manual")
            logging.info(f"Modelo cargado desde disco (sin indice previo): {modelo_id}")
            return modelo, f"Modelo '{modelo_id}' cargado desde disco local (nuevo indexado)."

    bigg_ok = _bigg_online()
    if not bigg_ok:
        console.print("[yellow]WARN[/yellow] BiGG no responde. Probando espejos directos...")

    try:
        modelo = load_model(modelo_id)
        try:
            write_sbml_model(modelo, ruta_local)
            _registrar_modelo(modelo, modelo_id, ruta_local, fuente="cobrapy_bigg")
        except Exception:
            pass
        logging.info(f"Modelo descargado via cobrapy: {modelo_id}")
        return modelo, f"Modelo '{modelo_id}' descargado via cobrapy (BiGG/BioModels)."
    except Exception as e_orig:
        console.print(f"[yellow]WARN[/yellow] load_model fallo: {str(e_orig)[:120]}")
        console.print("[cyan]INFO[/cyan] Probando descarga directa con headers de navegador...")
        logging.warning(f"load_model fallo para {modelo_id}: {e_orig}")

    errores = []
    for url in _urls_candidatas(modelo_id):
        ok, detalle = _descargar_xml(url, ruta_local, max_intentos=2, timeout=25)
        if ok:
            try:
                modelo = read_sbml_model(ruta_local)
                _registrar_modelo(modelo, modelo_id, ruta_local, fuente="espejo_directo")
                logging.info(f"Modelo descargado via espejo: {url}")
                return modelo, f"Modelo '{modelo_id}' descargado desde espejo directo (indexado)."
            except Exception as e_read:
                errores.append(f"{url}\n      lectura SBML fallo: {str(e_read)[:100]}")
        else:
            errores.append(f"{url}\n      {detalle}")

    detalle = "\n".join(f"  - {e}" for e in errores)
    ruta_abs = os.path.abspath(ruta_local)
    logging.error(f"Fallo total al cargar {modelo_id}: {detalle}")
    raise Exception(
        f"No se pudo cargar '{modelo_id}' desde ninguna fuente.\n\n"
        f"CAUSA PROBABLE: servidores de BiGG caidos.\n\n"
        f"Intentos:\n{detalle}\n\n"
        f"========== SOLUCION MANUAL ==========\n"
        f"1. Abre: http://bigg.ucsd.edu/models/{modelo_id}\n"
        f"2. Click 'Download XML'\n"
        f"3. Guarda en: {ruta_abs}\n"
        f"4. Vuelve a ejecutar; se indexara automaticamente.\n\n"
        f"========== ALTERNATIVA RAPIDA ==========\n"
        f"Opcion 1 del menu ('textbook') no requiere internet."
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
        logging.info("Ruta de cafeina inyectada")
        return modelo, "EX_caft_e"
    return None, None


# ==========================================
# ANALISIS TERMODINAMICO
# ==========================================
def analizar_termodinamica_cepa(modelo, quimico, genes_knockout):
    if not MODULOS_EXTERNOS_DISPONIBLES:
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
        return
    try:
        analizar_termodinamica_universal(modelo, reacciones_relevantes[:10])
    except Exception as e:
        console.print(f"[yellow]WARN[/yellow] Termodinamica fallo: {e}")
        logging.warning(f"Termodinamica fallo: {e}")


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
        logging.info(f"Analisis profundo: produccion={produccion:.4f}, penalizacion={costo:.2f}, carga={carga:.2f}")
        analizar_termodinamica_cepa(modelo_mut, quimico, genes_a_mutar)


# ==========================================
# DISENO DE CEPA POR IA
# ==========================================
def disenar_cepa_ia(modelo, nombre_organismo, modelo_id, seed=None):
    if seed is not None:
        random.seed(seed)
        logging.info(f"IA iniciada con seed={seed}")

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

    logging.info(f"IA params: quimico={quimico}, gens={gens}, pop={pop_size}, mut={mutaciones}")

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
            padre1 = random.choice(mejores)["genes"]
            padre2 = random.choice(mejores)["genes"]
            hijo = list(set(padre1[:3] + padre2[3:]))
            if len(hijo) < mutaciones: hijo.append(random.choice(todos_los_genes))
            fitness_hijo = evaluar_bacteria(modelo, hijo, quimico, gen_biomasa)
            nueva.append({"genes": hijo, "fitness": fitness_hijo})
        poblacion = nueva
        logging.debug(f"Gen {gen+1}/{gens}: best={poblacion[0]['fitness']:.4f}")

    mejor = sorted(poblacion, key=lambda x: x["fitness"], reverse=True)[0]

    console.print(f"\n[green]OK[/green] Diseno completado.")
    console.print(f"Fitness: {mejor['fitness']:.4f}\n")
    logging.info(f"IA completada: fitness={mejor['fitness']:.4f}, genes={mejor['genes']}")

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

    if REPORTE_DISPONIBLE:
        try:
            generar_documento_markdown(modelo, modelo_id, quimico, mejor["genes"], mejor["fitness"])
        except Exception as e:
            console.print(f"[yellow]WARN[/yellow] Reporte fallo: {e}")
            logging.warning(f"Reporte fallo: {e}")

    # EXPORTAR JSON + CSV
    try:
        _exportar_resultado_json(
            modelo_id=modelo_id,
            quimico=quimico,
            metodo="ia",
            genes=mejor["genes"],
            fitness=mejor["fitness"],
            seed=seed,
            parametros_extra={"generaciones": gens, "poblacion": pop_size, "mutaciones": mutaciones}
        )
    except Exception as e:
        console.print(f"[yellow]WARN[/yellow] Exportacion JSON/CSV fallo: {e}")
        logging.warning(f"Exportacion fallo: {e}")

    analisis_profundo(modelo, gen_biomasa, mejor["genes"], quimico)


# ==========================================
# DISENO POR OPTKNOCK
# ==========================================
def optknock_algoritmo(modelo, nombre_organismo, modelo_id, seed=None):
    if seed is not None:
        random.seed(seed)
        logging.info(f"OptKnock iniciado con seed={seed}")

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

    max_mutaciones = IntPrompt.ask("Cuantos genes apagar simultaneamente?", default=3)
    generaciones = IntPrompt.ask("Generaciones de la IA", default=15)
    pop_size = IntPrompt.ask("Tamano de poblacion", default=30)

    logging.info(f"OptKnock params: quimico={quimico}, mut={max_mutaciones}, gen={generaciones}, pop={pop_size}")

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
        logging.info(f"OptKnock completado: ratio={mejor['fitness']:.4f}, genes={mejor['genes']}")

        # EXPORTAR JSON + CSV
        try:
            _exportar_resultado_json(
                modelo_id=modelo_id,
                quimico=quimico,
                metodo="optknock",
                genes=mejor["genes"],
                fitness=mejor["fitness"],
                seed=seed,
                parametros_extra={"max_mutaciones": max_mutaciones, "generaciones": generaciones, "poblacion": pop_size}
            )
        except Exception as e:
            console.print(f"[yellow]WARN[/yellow] Exportacion JSON/CSV fallo: {e}")
            logging.warning(f"Exportacion fallo: {e}")

        analisis_profundo(modelo, gen_biomasa, mejor["genes"], quimico)
    else:
        console.print("\n[red]ERROR[/red] No se encontro acoplamiento estricto.")
        console.print("Sugerencia: metabolismos flexibles requieren el Algoritmo Genetico (Opcion 2).")
        logging.warning("OptKnock: no se encontro acoplamiento")


# ==========================================
# ANALISIS INTEGRAL
# ==========================================
def analisis_integral(modelo, nombre_organismo, modelo_id, seed=None):
    if seed is not None:
        random.seed(seed)
        logging.info(f"Analisis Integral iniciado con seed={seed}")

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

    logging.info(f"Integral params: quimico={quimico}, gen={gens}, pop={pop_size}, mut={mutaciones}")

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
    logging.info(f"Fase 1 completada: fitness={mejor['fitness']:.4f}, genes={mejor['genes']}")

    console.print("\n[cyan]INFO[/cyan] Fase 2: Analisis Fisicoquimico...")
    analisis_profundo(modelo, gen_biomasa, mejor["genes"], quimico)

    if MODULOS_EXTERNOS_DISPONIBLES:
        console.print("\n[cyan]INFO[/cyan] Fase 3: Generando Mapa Visual...")
        try:
            generar_mapa_universal(modelo_id, quimico)
            logging.info("Mapa visual generado")
        except Exception as e:
            console.print(f"[yellow]WARN[/yellow] Mapa visual fallo: {e}")
            logging.warning(f"Mapa visual fallo: {e}")
    else:
        console.print("\n[yellow]WARN[/yellow] Fase 3 omitida: modulo visual no disponible.")

    if REPORTE_DISPONIBLE:
        console.print("\n[cyan]INFO[/cyan] Fase 4: Generando Reporte Cientifico...")
        try:
            generar_documento_markdown(modelo, modelo_id, quimico, mejor["genes"], mejor["fitness"])
            logging.info("Reporte cientifico generado")
        except Exception as e:
            console.print(f"[yellow]WARN[/yellow] Reporte fallo: {e}")
            logging.warning(f"Reporte fallo: {e}")

    # EXPORTAR JSON + CSV
    try:
        _exportar_resultado_json(
            modelo_id=modelo_id,
            quimico=quimico,
            metodo="integral",
            genes=mejor["genes"],
            fitness=mejor["fitness"],
            seed=seed,
            parametros_extra={"generaciones": gens, "poblacion": pop_size, "mutaciones": mutaciones}
        )
    except Exception as e:
        console.print(f"[yellow]WARN[/yellow] Exportacion JSON/CSV fallo: {e}")
        logging.warning(f"Exportacion fallo: {e}")

    console.print("\n[green]OK[/green] Analisis Integral completado.")


# ==========================================
# MENU PRINCIPAL
# ==========================================
def main():
    console.clear()
    console.print(Panel.fit(
        "[bold green]MetaboEngineer Pro V5.5[/bold green]\n"
        "Suite Universal de Biologia de Sistemas",
        style="bold green"
    ))

    log_file = setup_logger()
    console.print(f"[dim]Log: {log_file}[/dim]\n")

    _sincronizar_indice()

    console.print("\n[bold]Selecciona el organismo base:[/bold]")
    for key, val in ORGANISMOS.items():
        console.print(f"  {key}. {val['nombre']}")

    cacheados = listar_modelos_cacheados()
    claves_cache = []
    if cacheados:
        console.print("\n[bold cyan]Modelos ya en cache local (listos sin internet):[/bold cyan]")
        for i, (mid, meta) in enumerate(cacheados.items(), start=1):
            alias = f" (alias de {meta['alias_de']})" if meta.get('alias_de') else ""
            console.print(f"  \\[c{i}] {mid} — {meta['reacciones']} reacciones, "
                          f"{meta['genes']} genes, {meta['fecha'][:10]}{alias}")
            claves_cache.append(f"c{i}")
        claves_cache.append("c0")
        console.print("  \\[c0] Seleccion normal (descargar nuevo)")

    opciones_validas = list(ORGANISMOS.keys()) + claves_cache

    try:
        opcion = Prompt.ask("Opcion", choices=opciones_validas, default="2")
    except KeyboardInterrupt:
        console.print("\n[yellow]INFO[/yellow] Seleccion cancelada.")
        return

    if opcion.startswith("c") and opcion != "c0":
        idx = int(opcion[1:]) - 1
        modelo_id = list(cacheados.keys())[idx]
        nombre_organismo = modelo_id
        console.print(f"[cyan]INFO[/cyan] Usando modelo cacheado: {modelo_id}")
    else:
        if opcion == "c0":
            try:
                opcion = Prompt.ask("Opcion", choices=list(ORGANISMOS.keys()), default="2")
            except KeyboardInterrupt:
                console.print("\n[yellow]INFO[/yellow] Seleccion cancelada.")
                return

        OTRO_KEY = next((k for k, v in ORGANISMOS.items() if v["id"] is None), None)

        if opcion == OTRO_KEY:
            modelo_id = Prompt.ask("ID del modelo (BiGG / BioModels / local)").strip()
            if not modelo_id:
                console.print("[red]ERROR[/red] No ingresaste un ID. Abortando.", style="bold red")
                return
            nombre_organismo = modelo_id
        else:
            modelo_id = ORGANISMOS[opcion]["id"]
            nombre_organismo = ORGANISMOS[opcion]["nombre"]

    seed = preguntar_seed()
    logging.info(f"Organismo seleccionado: {modelo_id}, seed={seed}")

    with console.status(f"[bold yellow]Cargando modelo {modelo_id}...[/bold yellow]"):
        try:
            modelo, mensaje = cargar_modelo_universal(modelo_id)
            console.print(f"[green]OK[/green] {mensaje} | {len(modelo.reactions)} reacciones.\n")
        except Exception as e:
            console.print(Panel(str(e), title="ERROR de carga", border_style="red"))
            logging.error(f"Error cargando modelo: {e}")
            return

    while True:
        console.print("\n[bold]--- MENU PRINCIPAL (V5.5) ---[/bold]")
        console.print("1. [cyan]Analisis de Flujos Basico (FBA)[/cyan]")
        console.print("2. [green]Diseno de Cepas por IA (+ Reporte + JSON/CSV)[/green]")
        console.print("3. [magenta]Diseno por OptKnock (Acoplamiento de Crecimiento)[/magenta]")
        console.print("4. [bold blue]Analisis Integral (Diseno + Termo + Mapa + Reporte)[/bold blue]")
        console.print("5. [bold cyan]Ver historial de corridas (CSV)[/bold cyan]")
        console.print("6. [bold yellow]Limpiar cache UniProt[/bold yellow]")
        console.print("7. [bold red]Salir[/bold red]")

        try:
            op = Prompt.ask("\nSelecciona una opcion", choices=["1", "2", "3", "4", "5", "6", "7"])
        except KeyboardInterrupt:
            console.print("\n[yellow]INFO[/yellow] Saliendo.")
            break

        if op == "1":
            sol = modelo.optimize()
            if sol.status == 'optimal':
                console.print(f"[green]OK[/green] Crecimiento: [bold green]{sol.objective_value:.4f}[/bold green] h-1")
                logging.info(f"FBA: crecimiento={sol.objective_value:.4f}")
        elif op == "2":
            disenar_cepa_ia(modelo, nombre_organismo, modelo_id, seed=seed)
        elif op == "3":
            optknock_algoritmo(modelo, nombre_organismo, modelo_id, seed=seed)
        elif op == "4":
            analisis_integral(modelo, nombre_organismo, modelo_id, seed=seed)
        elif op == "5":
            mostrar_historial()
        elif op == "6":
            limpiar_cache_uniprot()
        elif op == "7":
            console.print("Saliendo. Hasta luego.")
            logging.info("Session ended")
            break


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        console.print("\n[yellow]INFO[/yellow] Programa interrumpido por el usuario.")
        logging.info("Interrupted by user")
    except Exception as e:
        console.print(f"\n[red]FATAL ERROR[/red] {e}")
        logging.exception("Fatal error")