#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
MetaboEngineer Pro - CLI Experto
=================================
Interfaz de línea de comandos no interactiva para automatizacion.
Ideal para servidores, HPC (SLURM/PBS) y pipelines CI/CD.

Ejemplos de uso:
----------------
# Corrida basica con seed reproducible
python cli_experto.py --modelo textbook --quimico EX_succ_e --metodo ia --seed 42

# Con parametros personalizados
python cli_experto.py -m iML1515 -q EX_succ_e --metodo optknock --seed 42 \\
                       --generaciones 20 --poblacion 30 --mutaciones 4

# Inyectar ruta de cafeina automaticamente
python cli_experto.py -m textbook -q EX_caft_e --metodo ia --inyectar --seed 42

# Batch (multiples seeds para estadistica)
for seed in 1 2 3 4 5; do
    python cli_experto.py -m textbook -q EX_succ_e --metodo ia --seed $seed
done
"""

import warnings
warnings.filterwarnings("ignore")

import argparse
import logging
import os
import sys
import random
import traceback
from datetime import datetime

# Imports del core
from metabo_engineer import (
    cargar_modelo_universal,
    inyectar_ruta_heterologa,
    evaluar_bacteria,
    get_organism_id,
    traducir_gen_uniprot,
    analisis_profundo,
    analizar_termodinamica_cepa,
    _exportar_resultado_json,
    REPORTE_DISPONIBLE,
    MODULOS_EXTERNOS_DISPONIBLES,
)

if REPORTE_DISPONIBLE:
    from reporte_investigador import generar_documento_markdown


# ==========================================
# ARGPARSE
# ==========================================
def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="cli_experto.py",
        description="MetaboEngineer Pro - CLI no interactivo para automatizacion.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Ejemplos:
  # Corrida basica
  python cli_experto.py -m textbook -q EX_succ_e --metodo ia --seed 42

  # OptKnock con parametros personalizados
  python cli_experto.py -m iML1515 -q EX_succ_e --metodo optknock --seed 42 \\
                        --generaciones 20 --poblacion 30 --mutaciones 4

  # Inyectar ruta de cafeina
  python cli_experto.py -m textbook -q EX_caft_e --metodo ia --inyectar --seed 42

  # Batch (multiples seeds)
  for s in 1 2 3; do python cli_experto.py -m textbook -q EX_succ_e --metodo ia --seed $s; done
        """,
    )

    # --- Requeridos ---
    parser.add_argument(
        "-m", "--modelo", required=True,
        help="ID del modelo (ej. textbook, iML1515, iYO844) o ruta a un .xml"
    )
    parser.add_argument(
        "-q", "--quimico", required=True,
        help="Reaccion objetivo (ej. EX_succ_e, EX_caft_e)"
    )

    # --- Metodo ---
    parser.add_argument(
        "--metodo", choices=["ia", "optknock", "integral"], default="ia",
        help="Algoritmo a usar: ia (algoritmo genetico), optknock (acoplamiento), integral (todo). Default: ia"
    )

    # --- Parametros del algoritmo ---
    parser.add_argument(
        "--generaciones", type=int, default=10,
        help="Numero de generaciones de la IA (default: 10)"
    )
    parser.add_argument(
        "--poblacion", type=int, default=20,
        help="Tamano de la poblacion (default: 20)"
    )
    parser.add_argument(
        "--mutaciones", type=int, default=4,
        help="Genes a apagar simultaneamente (default: 4)"
    )

    # --- Reproducibilidad ---
    parser.add_argument(
        "--seed", type=int, default=None,
        help="Semilla aleatoria (int). Si se omite, se usa una aleatoria. Guardala para reproducir."
    )

    # --- Biologia sintetica ---
    parser.add_argument(
        "--inyectar", action="store_true",
        help="Inyectar ruta heterologa si el quimico no esta en el modelo (ej. EX_caft_e)"
    )

    # --- Salida ---
    parser.add_argument(
        "-o", "--output", default="reportes",
        help="Directorio de salida para JSON/CSV/reporte (default: reportes)"
    )
    parser.add_argument(
        "--no-reporte", action="store_true",
        help="No generar reporte .md (mas rapido en batch)"
    )
    parser.add_argument(
        "--no-termino", action="store_true",
        help="No ejecutar analisis termodinamico (mas rapido en batch)"
    )

    # --- Logging ---
    parser.add_argument(
        "--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"], default="INFO",
        help="Nivel de logging (default: INFO)"
    )

    # --- Version ---
    parser.add_argument(
        "--version", action="version", version="MetaboEngineer Pro CLI v5.5"
    )

    return parser


# ==========================================
# LOGGING
# ==========================================
def setup_logging(log_level_str: str) -> str:
    os.makedirs("logs", exist_ok=True)
    log_file = f"logs/cli_{datetime.now():%Y%m%d_%H%M%S}.log"
    logging.basicConfig(
        level=getattr(logging, log_level_str),
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[logging.FileHandler(log_file, encoding="utf-8")],
    )
    return log_file


# ==========================================
# METODOS
# ==========================================
def metodo_ia(modelo, quimico, args):
    """Diseno de cepa por algoritmo genetico."""
    gen_biomasa = [r.id for r in modelo.reactions if "biomass" in r.id.lower()][0]
    modelo.reactions.get_by_id(quimico).bounds = (0, 1000.0)
    todos_los_genes = [g.id for g in modelo.genes]

    poblacion = []
    for _ in range(args.poblacion):
        genes = random.sample(todos_los_genes, args.mutaciones)
        fit = evaluar_bacteria(modelo, genes, quimico, gen_biomasa)
        poblacion.append({"genes": genes, "fitness": fit})

    for gen in range(args.generaciones):
        poblacion.sort(key=lambda x: x["fitness"], reverse=True)
        mejores = poblacion[:5]
        nueva = mejores.copy()
        while len(nueva) < args.poblacion:
            p1 = random.choice(mejores)["genes"]
            p2 = random.choice(mejores)["genes"]
            hijo = list(set(p1[:3] + p2[3:]))
            if len(hijo) < args.mutaciones:
                hijo.append(random.choice(todos_los_genes))
            nueva.append({"genes": hijo, "fitness": evaluar_bacteria(modelo, hijo, quimico, gen_biomasa)})
        poblacion = nueva
        logging.debug(f"Gen {gen+1}/{args.generaciones}: best={poblacion[0]['fitness']:.4f}")

    mejor = sorted(poblacion, key=lambda x: x["fitness"], reverse=True)[0]
    logging.info(f"IA completada: fitness={mejor['fitness']:.4f}, genes={mejor['genes']}")
    return mejor


def metodo_optknock(modelo, quimico, args):
    """OptKnock: acoplamiento crecimiento-produccion."""
    gen_biomasa = [r.id for r in modelo.reactions if "biomass" in r.id.lower()][0]
    modelo.reactions.get_by_id(quimico).bounds = (0, 1000.0)

    sol_wt = modelo.optimize()
    activas = {r for r, v in sol_wt.fluxes.items() if abs(v) > 0.01 and not r.startswith("EX_")}
    candidatos = list({
        g.id for r in activas if r in modelo.reactions
        for g in modelo.reactions.get_by_id(r).genes
    })

    if len(candidatos) < 2:
        raise RuntimeError(f"Solo {len(candidatos)} genes activos (minimo 2).")

    def evaluar(genes):
        m = modelo.copy()
        for g in genes:
            if g in m.genes:
                m.genes.get_by_id(g).knock_out()
        m.objective = gen_biomasa
        s = m.optimize()
        if s.status == "optimal" and s.objective_value > 0.01:
            p = s.fluxes.get(quimico, 0.0)
            if p > 0:
                return p / s.objective_value
        return 0.0

    poblacion = []
    for _ in range(args.poblacion):
        genes = random.sample(candidatos, min(args.mutaciones, len(candidatos)))
        poblacion.append({"genes": genes, "fitness": evaluar(genes)})

    for gen in range(args.generaciones):
        poblacion.sort(key=lambda x: x["fitness"], reverse=True)
        mejores = poblacion[:5]
        nueva = mejores.copy()
        while len(nueva) < args.poblacion:
            p1 = random.choice(mejores)["genes"]
            p2 = random.choice(mejores)["genes"]
            hijo = list(set(p1[:2] + p2[2:]))
            if len(hijo) < args.mutaciones:
                hijo.append(random.choice(candidatos))
            nueva.append({"genes": hijo, "fitness": evaluar(hijo)})
        poblacion = nueva
        logging.debug(f"Gen {gen+1}/{args.generaciones}: best={poblacion[0]['fitness']:.4f}")

    mejor = sorted(poblacion, key=lambda x: x["fitness"], reverse=True)[0]
    logging.info(f"OptKnock completado: ratio={mejor['fitness']:.4f}, genes={mejor['genes']}")
    return mejor


# ==========================================
# MAIN
# ==========================================
def main():
    parser = build_parser()
    args = parser.parse_args()

    log_file = setup_logging(args.log_level)
    logging.info("=" * 60)
    logging.info("MetaboEngineer Pro CLI v5.5 - Session started")
    logging.info(f"Args: {vars(args)}")
    logging.info("=" * 60)

    # --- Validar directorio de salida ---
    os.makedirs(args.output, exist_ok=True)

    # --- Seed ---
    if args.seed is None:
        args.seed = random.randint(0, 999999)
        logging.info(f"Seed aleatoria: {args.seed}")
        print(f"[INFO] Seed generada: {args.seed}")
    else:
        logging.info(f"Seed fija: {args.seed}")
    random.seed(args.seed)

    # --- Cargar modelo ---
    print(f"[INFO] Cargando modelo '{args.modelo}'...")
    try:
        modelo, mensaje = cargar_modelo_universal(args.modelo)
        logging.info(f"Modelo cargado: {mensaje}")
        print(f"[OK] {mensaje} | {len(modelo.reactions)} reacciones")
    except Exception as e:
        logging.error(f"Fallo al cargar modelo: {e}")
        print(f"[ERROR] {e}", file=sys.stderr)
        return 2

    # --- Inyectar ruta si es necesario ---
    if args.quimico not in modelo.reactions:
        if args.inyectar and args.quimico == "EX_caft_e":
            print(f"[INFO] Inyectando ruta heterologa para {args.quimico}...")
            modelo, args.quimico = inyectar_ruta_heterologa(modelo, args.quimico)
            if modelo is None:
                print(f"[ERROR] No se pudo inyectar la ruta.", file=sys.stderr)
                return 3
        else:
            print(f"[ERROR] '{args.quimico}' no existe en el modelo. Usa --inyectar si es una ruta heterologa.",
                  file=sys.stderr)
            return 3

    # --- Ejecutar metodo ---
    print(f"[INFO] Ejecutando metodo '{args.metodo}' (seed={args.seed}, gen={args.generaciones}, "
          f"pop={args.poblacion}, mut={args.mutaciones})...")

    try:
        if args.metodo == "ia":
            mejor = metodo_ia(modelo, args.quimico, args)
        elif args.metodo == "optknock":
            mejor = metodo_optknock(modelo, args.quimico, args)
        elif args.metodo == "integral":
            # Para integral en CLI, hacemos IA + termo (sin mapa visual que requiere interfaz)
            mejor = metodo_ia(modelo, args.quimico, args)
            logging.info("Modo integral: se ejecutara termo despues de IA")
        else:
            print(f"[ERROR] Metodo desconocido: {args.metodo}", file=sys.stderr)
            return 4
    except Exception as e:
        logging.error(f"Fallo en el metodo {args.metodo}: {e}")
        logging.error(traceback.format_exc())
        print(f"[ERROR] {e}", file=sys.stderr)
        return 4

    # --- Resumen ---
    print(f"\n[OK] {args.metodo.upper()} completado")
    print(f"     Fitness: {mejor['fitness']:.4f} mmol/gDW/h")
    print(f"     Genes ({len(mejor['genes'])}): {', '.join(mejor['genes'])}")

    # --- Traducir nombres (para salida) ---
    id_uniprot = get_organism_id(args.modelo)
    print(f"\n[INFO] Traduciendo genes (UniProt taxon={id_uniprot})...")
    for g in mejor["genes"]:
        nombre, funcion = traducir_gen_uniprot(g, id_uniprot)
        print(f"       {g:>15} | {nombre:<15} | {funcion[:50]}")

    # --- Analisis profundo + termo ---
    if not args.no_termino:
        gen_biomasa = [r.id for r in modelo.reactions if "biomass" in r.id.lower()][0]
        try:
            print(f"\n[INFO] Analisis profundo + termodinamica...")
            analisis_profundo(modelo, gen_biomasa, mejor["genes"], args.quimico)
        except Exception as e:
            logging.warning(f"Analisis profundo fallo: {e}")
            print(f"[WARN] Analisis profundo fallo: {e}")

    # --- Reporte .md ---
    if not args.no_reporte and REPORTE_DISPONIBLE:
        try:
            print(f"\n[INFO] Generando reporte cientifico...")
            generar_documento_markdown(modelo, args.modelo, args.quimico, mejor["genes"], mejor["fitness"])
        except Exception as e:
            logging.warning(f"Reporte fallo: {e}")
            print(f"[WARN] Reporte fallo: {e}")

    # --- Exportar JSON/CSV ---
    try:
        parametros = {
            "generaciones": args.generaciones,
            "poblacion": args.poblacion,
            "mutaciones": args.mutaciones,
            "cli": True,
        }
        _exportar_resultado_json(
            modelo_id=args.modelo,
            quimico=args.quimico,
            metodo=args.metodo,
            genes=mejor["genes"],
            fitness=mejor["fitness"],
            seed=args.seed,
            parametros_extra=parametros,
        )
        print(f"\n[OK] Resultados exportados a '{args.output}/'")
    except Exception as e:
        logging.warning(f"Exportacion fallo: {e}")
        print(f"[WARN] Exportacion fallo: {e}")

    # --- Exit code ---
    if mejor["fitness"] > 0:
        logging.info("Session ended OK")
        print(f"\n[OK] Sesion completada. Log: {log_file}")
        return 0
    else:
        logging.warning("Fitness = 0 (no se encontro solucion)")
        print(f"\n[WARN] Fitness = 0. No se encontro cepa viable. Log: {log_file}")
        return 5


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("\n[INFO] Interrumpido por el usuario.", file=sys.stderr)
        sys.exit(130)
    except Exception as e:
        print(f"\n[FATAL] {e}", file=sys.stderr)
        logging.exception("Fatal error")
        sys.exit(1)