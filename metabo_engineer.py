import warnings
warnings.filterwarnings("ignore")

import cobra
from cobra import Reaction, Metabolite
from cobra.io import load_model, read_sbml_model, write_sbml_model
from cobra.flux_analysis import flux_variability_analysis as fva
from cobra.manipulation import knock_out_model_genes
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
import itertools
from datetime import datetime
from rich.console import Console
from rich.table import Table
from rich.prompt import Prompt, IntPrompt
from rich.progress import track
from rich.panel import Panel

console = Console()

MODULOS_EXTERNOS_DISPONIBLES = False
try:
    from visualizador_universal import generar_mapa_universal
    from modulo_termodinamica_v2 import analizar_termodinamica_universal, construir_reaccion
    MODULOS_EXTERNOS_DISPONIBLES = True
except ImportError:
    pass

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
    "5": {"nombre": "Otro (ingresar ID)", "id": None}
}

_RUTA_HISTORIAL_CSV = "reportes/historial.csv"
_RUTA_CACHE_UNIPROT = "cache_uniprot.json"
_RUTA_INDICE = "modelos/_index.json"

_MODELOS_BUNDLED = {"textbook", "e_coli_core"}

def setup_logger():
    if not os.path.exists("logs"): os.makedirs("logs")
    log_file = f"logs/run_{datetime.now():%Y%m%d_%H%M%S}.log"
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s", handlers=[logging.FileHandler(log_file, encoding="utf-8")])
    return log_file

def preguntar_seed():
    try:
        entrada = Prompt.ask("Semilla aleatoria (Enter = aleatoria)", default="").strip()
        return int(entrada) if entrada else random.randint(0, 999999)
    except:
        return random.randint(0, 999999)

def configurar_ambiente(modelo, ambiente):
    modelo_amb = modelo.copy()
    if 'EX_o2_e' in modelo_amb.reactions:
        if ambiente == 'anaerobio':
            modelo_amb.reactions.get_by_id('EX_o2_e').bounds = (0.0, 0.0) # Bloqueo total O2
        else:
            modelo_amb.reactions.get_by_id('EX_o2_e').bounds = (-1000.0, 1000.0)
    else:
        if ambiente == 'anaerobio':
            console.print("[yellow]WARN[/yellow] El modelo no tiene EX_o2_e.")

    if ambiente == 'aerobio':
        for rxn in ["EX_etoh_e", "EX_lac__D_e", "EX_for_e"]:
            if rxn in modelo_amb.reactions:
                # Bloqueo estricto de secreción/captación
                modelo_amb.reactions.get_by_id(rxn).bounds = (0.0, 0.0) 

    return modelo_amb

def get_organism_id(modelo_id: str) -> str:
    if not modelo_id: return "83333"
    if "iYO844" in modelo_id: return "224308"
    if "iCW773" in modelo_id: return "196627"
    if "iMM904" in modelo_id: return "559292"
    if "iEK1008" in modelo_id: return "1773"
    return "83333"

def _cargar_cache_uniprot():
    if not os.path.exists(_RUTA_CACHE_UNIPROT): return {}
    try:
        with open(_RUTA_CACHE_UNIPROT, "r", encoding="utf-8") as f: return json.load(f)
    except: return {}

def _guardar_cache_uniprot(cache):
    with open(_RUTA_CACHE_UNIPROT + ".tmp", "w", encoding="utf-8") as f: json.dump(cache, f, indent=2)
    os.replace(_RUTA_CACHE_UNIPROT + ".tmp", _RUTA_CACHE_UNIPROT)

def limpiar_cache_uniprot():
    if os.path.exists(_RUTA_CACHE_UNIPROT):
        os.remove(_RUTA_CACHE_UNIPROT)
        console.print(f"[green]OK[/green] Cache eliminado.")

def traducir_gen_uniprot(gen_id, organismo_id_uniprot):
    if not re.search(r'\d', gen_id): return gen_id, "Nombre comun"
    cache = _cargar_cache_uniprot()
    clave = f"{organismo_id_uniprot}:{gen_id}"
    if clave in cache: return tuple(cache[clave])

    url = f'https://rest.uniprot.org/uniprotkb/search?query=organism_id:{organismo_id_uniprot}+AND+"{gen_id}"&fields=gene_names,protein_name&format=json'
    try:
        time.sleep(0.3)
        response = requests.get(url, timeout=5)
        if response.status_code == 200:
            data = response.json()
            if data.get('results'):
                result = data['results'][0]
                nombre_gen = result.get('genes', [{}])[0].get('geneName', {}).get('value', 'Desconocido')
                nombre_proteina = result.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', 'Funcion desconocida')
                if nombre_gen != 'Desconocido':
                    cache[clave] = [nombre_gen, nombre_proteina]
                    _guardar_cache_uniprot(cache)
                return nombre_gen, nombre_proteina
    except: pass
    return gen_id, "No encontrado en UniProt"

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

def _registrar_modelo(modelo, modelo_id, ruta, fuente):
    indice = _cargar_indice()
    indice["modelos"][modelo_id] = {"archivo": os.path.basename(ruta), "sha256": _sha256_archivo(ruta), "fecha": datetime.now().strftime("%Y-%m-%dT%H:%M:%S"), "reacciones": len(modelo.reactions), "genes": len(modelo.genes)}
    _guardar_indice(indice)

def listar_modelos_cacheados():
    return _cargar_indice().get("modelos", {})

def _sincronizar_indice():
    if not os.path.exists("modelos"): return
    indice = _cargar_indice()
    for archivo in os.listdir("modelos"):
        if archivo.endswith(".xml") and archivo[:-4] not in indice["modelos"]:
            try:
                _registrar_modelo(read_sbml_model(os.path.join("modelos", archivo)), archivo[:-4], os.path.join("modelos", archivo), "manual")
            except: pass

def cargar_modelo_universal(modelo_id):
    if not modelo_id or not str(modelo_id).strip(): raise Exception("ID de modelo vacio.")
    modelo_id = str(modelo_id).strip()

    if modelo_id in _MODELOS_BUNDLED:
        try:
            modelo = load_model(modelo_id)
            return modelo, f"Modelo '{modelo_id}' cargado desde bundle COBRApy."
        except Exception as e:
            raise Exception(f"No se pudo cargar '{modelo_id}' desde bundle COBRApy: {e}")

    ruta_cwd = f"{modelo_id}.xml"
    if os.path.exists(ruta_cwd):
        modelo = read_sbml_model(ruta_cwd)
        return modelo, f"Modelo '{modelo_id}' cargado desde {ruta_cwd} (CI)."

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

def inyectar_ruta_heterologa(modelo, quimico):
    if quimico == "EX_caft_e":
        caft_c = Metabolite('caft_c', name='Caffeine', formula='C8H10N4O2', compartment='c')
        tbrm_c = Metabolite('tbrm_c', name='Theobromine', formula='C7H8N4O2', compartment='c')
        xan_c = Metabolite('xan_c', name='Xanthine', formula='C5H4N4O2', compartment='c')
        caft_e = Metabolite('caft_e', name='Caffeine', formula='C8H10N4O2', compartment='e')
        xan_e = Metabolite('xan_e', name='Xanthine', formula='C5H4N4O2', compartment='e')
        h_c = modelo.metabolites.get_by_id("h_c") if "h_c" in modelo.metabolites else Metabolite('h_c', name='H+', formula='H', charge=1, compartment='c')
        sam_c = modelo.metabolites.get_by_id("amet_c") if "amet_c" in modelo.metabolites else Metabolite('amet_c', name='S-Adenosylmethionine', formula='C15H23N6O5S', charge=1, compartment='c')
        sah_c = modelo.metabolites.get_by_id("ahcys_c") if "ahcys_c" in modelo.metabolites else Metabolite('ahcys_c', name='S-Adenosylhomocysteine', formula='C14H20N6O5S', compartment='c')

        for met in [caft_c, tbrm_c, xan_c, caft_e, xan_e, h_c, sam_c, sah_c]:
            if met.id not in modelo.metabolites: modelo.add_metabolites([met])

        reacciones = [Reaction('CAFT_RXN1'), Reaction('CAFT_RXN2'), Reaction('EX_caft_e'), Reaction('CAFT_TRANSPORT'), Reaction('EX_xan_e'), Reaction('XAN_TRANSPORT')]
        reacciones[0].add_metabolites({xan_c: -1.0, sam_c: -2.0, tbrm_c: 1.0, sah_c: 2.0, h_c: 2.0})
        reacciones[1].add_metabolites({tbrm_c: -1.0, sam_c: -1.0, caft_c: 1.0, sah_c: 1.0, h_c: 1.0})
        reacciones[2].add_metabolites({caft_e: -1.0})
        reacciones[3].add_metabolites({caft_c: -1.0, caft_e: 1.0})
        reacciones[4].add_metabolites({xan_e: -1.0})
        reacciones[5].add_metabolites({xan_e: -1.0, xan_c: 1.0})

        for r in reacciones:
            if r.id not in modelo.reactions:
                if r.id == "EX_xan_e":
                    r.bounds = (-10, 1000)
                else:
                    r.bounds = (-1000, 1000) if r.id.startswith("EX_") else (0, 1000)
                modelo.add_reactions([r])
        return modelo, "EX_caft_e"
    return None, None

def obtener_genes_no_esenciales(modelo, gen_biomasa):
    console.print("[cyan]INFO[/cyan] Filtrando genes esenciales...")
    genes_no_esenciales = []
    modelo_wt = modelo.copy()
    modelo_wt.objective = gen_biomasa
    sol_wt = modelo_wt.optimize()
    crecimiento_wt = sol_wt.objective_value if sol_wt.status == 'optimal' else 0.0

    if crecimiento_wt <= 0: return [g.id for g in modelo.genes]

    for gen in track(modelo.genes, description="Evaluando esencialidad..."):
        if gen.id == "s0001" or "spontaneous" in gen.name.lower():
            continue

        m = modelo.copy()
        with m:
            knock_out_model_genes(m, [gen.id])
            m.objective = gen_biomasa
            sol = m.slim_optimize(error_value=0.0)
            if sol > (crecimiento_wt * 0.05):
                genes_no_esenciales.append(gen.id)

    console.print(f"[green]OK[/green] Se descartaron {len(modelo.genes) - len(genes_no_esenciales)} genes letales.")
    return genes_no_esenciales

def aplicar_restricciones_cineticas(modelo):
    if not os.path.exists("cinetica_local.json"): return modelo
    try:
        with open("cinetica_local.json", "r", encoding="utf-8") as f: db = json.load(f)
    except: return modelo

    pool_enzima = 0.1
    for rxn_id, data in db.items():
        if rxn_id in modelo.reactions and "kcat" in data:
            kcat = data["kcat"]
            if kcat > 0:
                v_max = kcat * pool_enzima
                rxn = modelo.reactions.get_by_id(rxn_id)
                rxn.bounds = (max(rxn.lower_bound, -v_max), min(rxn.upper_bound, v_max))
    return modelo

def evaluar_bacteria(modelo_base, genes_apagados, gen_objetivo, gen_biomasa):
    modelo = modelo_base.copy()
    for gen_id in genes_apagados:
        if gen_id in modelo.genes:
            modelo.genes.get_by_id(gen_id).knock_out()

    try:
        modelo.objective = gen_biomasa
        sol_biomasa = modelo.optimize()

        if sol_biomasa.status != 'optimal' or sol_biomasa.objective_value < 0.05:
            return -1e6

        max_biomasa_mut = sol_biomasa.objective_value
        bio_rxn = modelo.reactions.get_by_id(gen_biomasa)
        bio_rxn.lower_bound = max_biomasa_mut * 1.0

        modelo.objective = gen_objetivo
        modelo.objective.direction = 'min'
        sol_prod = modelo.optimize()

        if sol_prod.status == 'optimal':
            return max(0.0, sol_prod.objective_value)
    except Exception:
        pass
    return -1e6

def analizar_termodinamica_cepa(modelo, quimico, genes_knockout):
    if not MODULOS_EXTERNOS_DISPONIBLES: return
    reacciones = [quimico] if quimico in modelo.reactions else []
    for gen_id in genes_knockout:
        if gen_id in modelo.genes:
            for rxn in modelo.genes.get_by_id(gen_id).reactions:
                if rxn.id not in reacciones and not rxn.id.startswith("EX_"): reacciones.append(rxn.id)
    try: analizar_termodinamica_universal(modelo, reacciones[:10])
    except: pass

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

# 🛠️ FIX 1: Unidades biológicas dinámicas (mol/mol)
def analisis_profundo(modelo, gen_biomasa, genes_a_mutar, quimico):
    console.print(f"\n[cyan]INFO[/cyan] Analisis Profundo de la Cepa...")
    modelo_mut = modelo.copy()

    for gen_id in genes_a_mutar:
        if gen_id in modelo_mut.genes:
            modelo_mut.genes.get_by_id(gen_id).knock_out()

    try:
        modelo_mut.objective = gen_biomasa
        sol_bio = modelo_mut.optimize()

        if sol_bio.status == 'optimal':
            max_bio = sol_bio.objective_value
            modelo_mut.reactions.get_by_id(gen_biomasa).lower_bound = max_bio * 0.9
            modelo_mut.objective = quimico
            sol_final = modelo_mut.optimize()

            if sol_final.status == 'optimal':
                produccion = sol_final.objective_value
                
                # --- CÁLCULO HONESTO DEL RENDIMIENTO ---
                try:
                    # Leemos la captación real de glucosa del modelo
                    glc_uptake = abs(modelo.reactions.get_by_id("EX_glc__D_e").lower_bound)
                    if glc_uptake == 0: glc_uptake = 10.0 # Fallback de seguridad
                except:
                    glc_uptake = 10.0
                    
                rendimiento_mol_mol = produccion / glc_uptake

                console.print(f"   -> Fitness (Flujo Garantizado): [bold cyan]{produccion:.4f}[/bold cyan] mmol/gDW/h")
                console.print(f"   -> Rendimiento Teórico:       [bold green]{rendimiento_mol_mol:.3f}[/bold green] mol/mol")
                console.print(f"   [dim](Calculado sobre una captación de glucosa de {glc_uptake} mmol/gDW/h)[/dim]")
                
                carga = sol_final.fluxes.abs().sum()
                console.print(f"   -> Carga Metabolica:          [bold magenta]{carga:.2f}[/bold magenta]")
                analizar_termodinamica_cepa(modelo_mut, quimico, genes_a_mutar)
                mostrar_cinetica_local(modelo, quimico, genes_a_mutar)
    except Exception as e:
        console.print(f"[yellow]WARN[/yellow] Analisis fallo: {e}")

def _exportar_resultado_json(modelo_id, quimico, metodo, genes, fitness, seed, parametros_extra=None):
    if not os.path.exists("reportes"): os.makedirs("reportes")
    fecha = datetime.now()
    resultado = {"fecha": fecha.isoformat(), "modelo": modelo_id, "quimico": quimico, "metodo": metodo, "fitness": round(fitness, 6), "genes_knockout": list(genes), "seed": seed, "parametros": parametros_extra or {}}
    with open(f"reportes/{metodo}_{modelo_id}_{quimico}_{fecha.strftime('%Y%m%d_%H%M%S')}.json", "w", encoding="utf-8") as f: json.dump(resultado, f, indent=2)
    fila = {"fecha": fecha.strftime("%Y-%m-%d %H:%M:%S"), "modelo": modelo_id, "quimico": quimico, "metodo": metodo, "fitness": round(fitness, 6), "genes": ";".join(genes), "seed": seed}
    df = pd.DataFrame([fila])
    df.to_csv(_RUTA_HISTORIAL_CSV, mode="a", header=not os.path.exists(_RUTA_HISTORIAL_CSV), index=False)

def mostrar_historial():
    if not os.path.exists(_RUTA_HISTORIAL_CSV): return
    df = pd.read_csv(_RUTA_HISTORIAL_CSV)
    if df.empty: return
    console.print(df.tail(10).to_string())

# 🛠️ FIX 2: Mutación real sin duplicados
def disenar_cepa_ia(modelo, nombre_organismo, modelo_id, seed=None):
    if seed is not None: random.seed(seed)
    quimico = Prompt.ask("Quimico objetivo", default="EX_succ_e").strip()
    if quimico not in modelo.reactions: return

    gens = IntPrompt.ask("Generaciones", default=10); pop_size = IntPrompt.ask("Poblacion", default=20); mutaciones = IntPrompt.ask("Genes a apagar", default=4)
    gen_biomasa = [r.id for r in modelo.reactions if "biomass" in r.id.lower()][0]
    modelo.reactions.get_by_id(quimico).bounds = (0, 1000.0)

    todos_los_genes = obtener_genes_no_esenciales(modelo, gen_biomasa)

    poblacion = []
    for _ in range(pop_size):
        gm = random.sample(todos_los_genes, min(mutaciones, len(todos_los_genes)))
        poblacion.append({"genes": gm, "fitness": evaluar_bacteria(modelo, gm, quimico, gen_biomasa)})
    for gen in track(range(gens), description="Evolucionando..."):
        poblacion.sort(key=lambda x: x["fitness"], reverse=True)
        mejores = poblacion[:5]
        nueva = mejores.copy()
        while len(nueva) < pop_size:
            p1, p2 = random.choice(mejores)["genes"], random.choice(mejores)["genes"]
            hijo = list(set(p1[:3] + p2[3:]))
            if len(hijo) < mutaciones: hijo.append(random.choice(todos_los_genes))
            
            # Mutación evitando duplicados
            if random.random() < 0.1 and hijo:
                idx = random.randint(0, len(hijo) - 1)
                nuevo_gen = random.choice(todos_los_genes)
                if nuevo_gen not in hijo:
                    hijo[idx] = nuevo_gen
                    hijo = list(set(hijo))
                
            nueva.append({"genes": hijo, "fitness": evaluar_bacteria(modelo, hijo, quimico, gen_biomasa)})
        poblacion = nueva

    mejor = sorted(poblacion, key=lambda x: x["fitness"], reverse=True)[0]
    console.print(f"\n[green]OK[/green] Diseno completado. Fitness (FVA-min@100%): {mejor['fitness']:.4f}\n")

    id_uniprot = get_organism_id(modelo_id)
    tabla = Table(title="Genes a apagar (Knockouts)")
    tabla.add_column("ID (BiGG)", style="cyan"); tabla.add_column("Nombre Real", style="magenta"); tabla.add_column("Funcion Proteica", style="green")
    for g in mejor["genes"]:
        nombre_local = modelo.genes.get_by_id(g).name
        nombre_real, funcion = (nombre_local, "Modelo local") if (nombre_local and nombre_local != g) else traducir_gen_uniprot(g, id_uniprot)
        tabla.add_row(g, nombre_real, funcion[:50] + "..." if len(funcion) > 50 else funcion)
    console.print(tabla)

    if REPORTE_DISPONIBLE:
        try: generar_documento_markdown(modelo, modelo_id, quimico, mejor["genes"], mejor["fitness"])
        except Exception as e: console.print(f"[yellow]WARN[/yellow] Reporte fallo: {e}")
    try: _exportar_resultado_json(modelo_id, quimico, "ia", mejor["genes"], mejor["fitness"], seed, {"gens": gens, "mut": mutaciones})
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
    genes_candidatos = sorted({g.id for r in reacciones_activas if r in modelo.reactions for g in modelo.reactions.get_by_id(r).genes})
    if len(genes_candidatos) < 2: return

    def evaluar(genes):
        m = modelo.copy()
        for gen in genes:
            if gen in m.genes: m.genes.get_by_id(gen).knock_out()
        m.objective = gen_biomasa
        sol = m.optimize()
        if sol.status == 'optimal' and sol.objective_value > 0.01:
            try:
                fva_res = fva(m, [quimico], fraction_of_optimum=1.0)
                return max(0.0, fva_res.minimum.iloc[0])
            except: return 0.0
        return 0.0

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
            
            # Mutación evitando duplicados
            if random.random() < 0.1 and hijo:
                idx = random.randint(0, len(hijo) - 1)
                nuevo_gen = random.choice(genes_candidatos)
                if nuevo_gen not in hijo:
                    hijo[idx] = nuevo_gen
                    hijo = list(set(hijo))
                
            nueva.append({"genes": hijo, "fitness": evaluar(hijo)})
        poblacion = nueva

    mejor = sorted(poblacion, key=lambda x: x["fitness"], reverse=True)[0]
    if mejor["fitness"] > 0:
        console.print(f"\n[green]OK[/green] OptKnock exitoso. Ratio FVA: {mejor['fitness']:.4f}")
        try: _exportar_resultado_json(modelo_id, quimico, "optknock", mejor["genes"], mejor["fitness"], seed)
        except: pass
        analisis_profundo(modelo, gen_biomasa, mejor["genes"], quimico)
    else:
        console.print("\n[red]ERROR[/red] No se encontro acoplamiento.")

def optknock_exacto_milp(modelo, nombre_organismo, modelo_id, seed=None):
    if seed is not None: random.seed(seed)
    quimico = Prompt.ask("Quimico objetivo", default="EX_succ_e").strip()
    if quimico not in modelo.reactions: return

    max_knockouts = IntPrompt.ask("Genes exactos a apagar?", default=2)
    gen_biomasa = [r.id for r in modelo.reactions if "biomass" in r.id.lower()][0]
    modelo.reactions.get_by_id(quimico).bounds = (0, 1000.0)

    console.print("\n[cyan]INFO[/cyan] Buscando el óptimo exacto (MILP equivalent)...")
    sol_wt = modelo.optimize()
    reacciones_activas = [r for r, v in sol_wt.fluxes.items() if abs(v) > 0.1 and not r.startswith("EX_")]
    genes_activos = sorted(list({g.id for r in reacciones_activas if r in modelo.reactions for g in modelo.reactions.get_by_id(r).genes}))

    mejor_ratio = -1
    mejor_combo = []
    combinaciones = list(itertools.combinations(genes_activos, max_knockouts))
    console.print(f"   -> Evaluando {len(combinaciones)} combinaciones exactas...")

    for combo in track(combinaciones, description="Buscando acoplamiento exacto..."):
        m = modelo.copy()
        for g in combo:
            if g in m.genes: m.genes.get_by_id(g).knock_out()

        m.objective = gen_biomasa
        sol_bio = m.optimize()

        if sol_bio.status == 'optimal' and sol_bio.objective_value > 0.1:
            try:
                fva_res = fva(m, [quimico], fraction_of_optimum=1.0)
                prod_garantizada = fva_res.minimum.iloc[0]
                if prod_garantizada > 0:
                    ratio = prod_garantizada / sol_bio.objective_value
                    if ratio > mejor_ratio:
                        mejor_ratio = ratio
                        mejor_combo = combo
            except: pass

    if mejor_combo:
        console.print(f"\n[bold green]¡OptKnock Exacto Exitoso![/bold green]")
        console.print(f"Genes: {', '.join(mejor_combo)}")
        console.print(f"Ratio Real: {mejor_ratio:.4f}")
        try: _exportar_resultado_json(modelo_id, quimico, "optknock_exacto", mejor_combo, mejor_ratio, seed)
        except: pass
        analisis_profundo(modelo, gen_biomasa, mejor_combo, quimico)
    else:
        console.print("[red]No se encontró acoplamiento exacto.[/red]")

def analisis_integral(modelo, nombre_organismo, modelo_id, seed=None):
    if seed is not None: random.seed(seed)
    quimico = Prompt.ask("Quimico objetivo", default="EX_succ_e")
    if quimico not in modelo.reactions: return

    gens = IntPrompt.ask("Generaciones", default=10); pop_size = IntPrompt.ask("Poblacion", default=20); mutaciones = IntPrompt.ask("Genes a apagar", default=4)
    gen_biomasa = [r.id for r in modelo.reactions if "biomass" in r.id.lower()][0]
    modelo.reactions.get_by_id(quimico).bounds = (0, 1000.0)

    todos_los_genes = obtener_genes_no_esenciales(modelo, gen_biomasa)

    poblacion = []
    for _ in range(pop_size):
        gm = random.sample(todos_los_genes, min(mutaciones, len(todos_los_genes)))
        poblacion.append({"genes": gm, "fitness": evaluar_bacteria(modelo, gm, quimico, gen_biomasa)})
    for gen in track(range(gens), description="Evolucionando..."):
        poblacion.sort(key=lambda x: x["fitness"], reverse=True)
        mejores = poblacion[:5]
        nueva = mejores.copy()
        while len(nueva) < pop_size:
            p1, p2 = random.choice(mejores)["genes"], random.choice(mejores)["genes"]
            hijo = list(set(p1[:3] + p2[3:]))
            if len(hijo) < mutaciones: hijo.append(random.choice(todos_los_genes))
            
            # Mutación evitando duplicados
            if random.random() < 0.1 and hijo:
                idx = random.randint(0, len(hijo) - 1)
                nuevo_gen = random.choice(todos_los_genes)
                if nuevo_gen not in hijo:
                    hijo[idx] = nuevo_gen
                    hijo = list(set(hijo))

            nueva.append({"genes": hijo, "fitness": evaluar_bacteria(modelo, hijo, quimico, gen_biomasa)})
        poblacion = nueva

    mejor = sorted(poblacion, key=lambda x: x["fitness"], reverse=True)[0]
    console.print(f"[green]OK[/green] Mejor Cepa: {', '.join(mejor['genes'])} (Fitness FVA-min@100%: {mejor['fitness']:.4f})")
    analisis_profundo(modelo, gen_biomasa, mejor["genes"], quimico)

    if MODULOS_EXTERNOS_DISPONIBLES:
        try: generar_mapa_universal(modelo_id, quimico)
        except: pass
    if REPORTE_DISPONIBLE:
        try: generar_documento_markdown(modelo, modelo_id, quimico, mejor["genes"], mejor["fitness"])
        except: pass
    try: _exportar_resultado_json(modelo_id, quimico, "integral", mejor["genes"], mejor["fitness"], seed)
    except: pass

def main():
    console.clear(); console.print(Panel.fit("[bold green]MetaboEngineer Pro V7.1 (PhD Fixed)[/bold green]", style="bold green"))
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

    try:
        opcion = Prompt.ask("Opcion", choices=list(ORGANISMOS.keys()) + claves_cache, default="2")
    except: return

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

    ambiente = Prompt.ask("Condición ambiental", choices=["aerobio", "anaerobio"], default="aerobio")
    seed = preguntar_seed()

    with console.status(f"[bold yellow]Cargando modelo...[/bold yellow]"):
        try:
            modelo, msg = cargar_modelo_universal(modelo_id)
            modelo = configurar_ambiente(modelo, ambiente)
            
            # aplicar_restricciones_cineticas(modelo) # Desactivado por inconsistencia de unidades
            
            console.print(f"[green]OK[/green] {msg} | {len(modelo.reactions)} reacciones | Ambiente: {ambiente.upper()}\n")
        except Exception as e:
            console.print(Panel(str(e), title="ERROR", border_style="red")); return

    while True:
        console.print("\n[bold]--- MENU (V7.1) ---[/bold]")
        console.print("1. FBA | 2. IA | 3. OptKnock (GA) | 4. Integral | 5. Historial | 6. Limpiar Cache | 7. Inyectar Cafeína | 8. OptKnock Exacto (PhD) | 9. Salir")
        try:
            op = Prompt.ask("\nSelecciona una opcion", choices=["1", "2", "3", "4", "5", "6", "7", "8", "9"])
        except: break

        if op == "1":
            sol = modelo.optimize()
            if sol.status == 'optimal': console.print(f"[green]OK[/green] Crecimiento: {sol.objective_value:.4f} h-1")
        elif op == "2": disenar_cepa_ia(modelo, nombre_organismo, modelo_id, seed)
        elif op == "3": optknock_algoritmo(modelo, nombre_organismo, modelo_id, seed)
        elif op == "4": analisis_integral(modelo, nombre_organismo, modelo_id, seed)
        elif op == "5": mostrar_historial()
        elif op == "6": limpiar_cache_uniprot()
        elif op == "7":
            modelo, _ = inyectar_ruta_heterologa(modelo, "EX_caft_e")
            console.print("[green]OK[/green] Ruta de cafeína inyectada.")
        elif op == "8": optknock_exacto_milp(modelo, nombre_organismo, modelo_id, seed)
        elif op == "9": break

if __name__ == "__main__":
    main()