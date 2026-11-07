#!/usr/bin/env python
# -*- coding: utf-8 -*-
import warnings
warnings.filterwarnings("ignore")

import argparse
import logging
import os
import sys
import random
import traceback
from datetime import datetime

from metabo_engineer import (
    _detectar_gen_biomasa,
    _evolucionar,
    cargar_modelo_universal,
    evaluar_acoplamiento_fva,
    evaluar_bacteria,
    genes_candidatos_activos,
    configurar_ambiente,
    obtener_genes_no_esenciales,
    get_organism_id,
    traducir_gen_uniprot,
    analisis_profundo,
    inyectar_ruta_heterologa
)

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="cli_experto.py",
        description="MetaboEngineer Pro - CLI no interactivo."
    )
    parser.add_argument("-m", "--modelo", required=True)
    parser.add_argument("-q", "--quimico", required=True)
    parser.add_argument("--metodo", choices=["ia", "optknock", "integral"], default="ia")
    parser.add_argument("--ambiente", choices=["aerobio", "anaerobio"], default="aerobio")
    parser.add_argument("--generaciones", type=int, default=10)
    parser.add_argument("--poblacion", type=int, default=20)
    parser.add_argument("--mutaciones", type=int, default=4)
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument("--inyectar", action="store_true")
    return parser

def setup_logging():
    os.makedirs("logs", exist_ok=True)
    log_file = f"logs/cli_{datetime.now():%Y%m%d_%H%M%S}.log"
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s", handlers=[logging.FileHandler(log_file, encoding="utf-8")])
    return log_file

def metodo_ia(modelo, quimico, args):
    gen_biomasa = _detectar_gen_biomasa(modelo)
    modelo.reactions.get_by_id(quimico).bounds = (0, 1000.0)

    todos_los_genes = obtener_genes_no_esenciales(modelo, gen_biomasa)
    cache = {}
    return _evolucionar(
        todos_los_genes, args.mutaciones, args.generaciones, args.poblacion,
        lambda g: evaluar_bacteria(modelo, g, quimico, gen_biomasa, _cache=cache),
        crossover_split=3,
    )

def metodo_optknock(modelo, quimico, args):
    gen_biomasa = _detectar_gen_biomasa(modelo)
    modelo.reactions.get_by_id(quimico).bounds = (0, 1000.0)

    genes_candidatos = genes_candidatos_activos(modelo)
    if len(genes_candidatos) < 2:
        raise RuntimeError("Pocos genes activos.")

    cache = {}
    return _evolucionar(
        genes_candidatos, args.mutaciones, args.generaciones, args.poblacion,
        lambda g: evaluar_acoplamiento_fva(modelo, g, quimico, gen_biomasa, _cache=cache),
        crossover_split=2,
    )

FITNESS_FALLO_TOTAL = -1e5  # evaluar_bacteria devuelve -1e6 cuando una evaluacion falla


def verificar_resultado(mejor, metodo):
    """Imprime el estado del resultado y devuelve el codigo de salida (0 = hay acoplamiento)."""
    if mejor["fitness"] <= FITNESS_FALLO_TOTAL:
        print("[ERROR] TODAS las evaluaciones fallaron (fitness = -1e6). "
              "Revisa el modelo, el ambiente y el quimico objetivo.", file=sys.stderr)
        return 5
    if mejor["fitness"] <= 0:
        print(f"\n[WARN] {metodo.upper()} termino sin encontrar acoplamiento "
              "(fitness 0). Prueba con mas generaciones, poblacion o genes a apagar.")
    else:
        print(f"\n[OK] {metodo.upper()} completado")
    return 0


def main():
    parser = build_parser()
    args = parser.parse_args()
    log_file = setup_logging()
    
    if args.seed is None: args.seed = random.randint(0, 999999)
    random.seed(args.seed)
    print(f"[INFO] Seed: {args.seed}")

    print(f"[INFO] Cargando modelo '{args.modelo}'...")
    try:
        modelo, mensaje = cargar_modelo_universal(args.modelo)
        print(f"[OK] {mensaje} | {len(modelo.reactions)} reacciones")
    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr); return 2

    print(f"[INFO] Configurando ambiente: {args.ambiente.upper()}...")
    try:
        modelo = configurar_ambiente(modelo, args.ambiente)
    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr); return 2

    if args.quimico not in modelo.reactions:
        if args.inyectar and args.quimico == "EX_caft_e":
            print(f"[INFO] Inyectando ruta heterologa para {args.quimico}...")
            modelo, args.quimico = inyectar_ruta_heterologa(modelo, args.quimico)
            if modelo is None: return 3
        else:
            print(f"[ERROR] '{args.quimico}' no existe. Usa --inyectar.", file=sys.stderr); return 3

    print(f"[INFO] Ejecutando metodo '{args.metodo}'...")

    try:
        if args.metodo == "ia":
            mejor = metodo_ia(modelo, args.quimico, args)
        elif args.metodo == "optknock":
            mejor = metodo_optknock(modelo, args.quimico, args)
        elif args.metodo == "integral":
            mejor = metodo_ia(modelo, args.quimico, args)
    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr); logging.error(traceback.format_exc()); return 4

    codigo = verificar_resultado(mejor, args.metodo)
    if codigo != 0:
        return codigo

    print(f"     Fitness: {mejor['fitness']:.4f} mmol/gDW/h")
    print(f"     Genes ({len(mejor['genes'])}): {', '.join(mejor['genes'])}")

    id_uniprot = get_organism_id(args.modelo)
    print(f"\n[INFO] Traduciendo genes (taxon={id_uniprot})...")
    for g in mejor["genes"]:
        nombre, funcion = traducir_gen_uniprot(g, id_uniprot)
        print(f"       {g:>15} | {nombre:<15} | {funcion[:50]}")

    gen_biomasa = _detectar_gen_biomasa(modelo)
    try:
        print("\n[INFO] Analisis profundo + termodinamica...")
        analisis_profundo(modelo, gen_biomasa, mejor["genes"], args.quimico)
    except Exception as e:
        print(f"[WARN] Analisis profundo fallo: {e}")

    if args.metodo == "integral":
        try:
            from reporte_investigador import generar_documento_markdown
            print("\n[INFO] Generando reporte cientifico...")
            generar_documento_markdown(modelo, args.modelo, args.quimico, mejor["genes"], mejor["fitness"])
        except Exception as e:
            print(f"[WARN] Reporte fallo: {e}")

    print(f"\n[OK] Sesion completada. Log: {log_file}")
    return 0

if __name__ == "__main__":
    try: sys.exit(main())
    except KeyboardInterrupt: sys.exit(130)