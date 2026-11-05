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
    cargar_modelo_universal,
    inyectar_ruta_heterologa,
    evaluar_bacteria,
    get_organism_id,
    traducir_gen_uniprot,
    analisis_profundo,
    _exportar_resultado_json,
    REPORTE_DISPONIBLE,
    configurar_ambiente,
)

if REPORTE_DISPONIBLE:
    from reporte_investigador import generar_documento_markdown

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="cli_experto.py",
        description="MetaboEngineer Pro - CLI no interactivo para automatizacion."
    )
    parser.add_argument("-m", "--modelo", required=True, help="ID del modelo (ej. textbook, iML1515)")
    parser.add_argument("-q", "--quimico", required=True, help="Reaccion objetivo (ej. EX_succ_e)")
    parser.add_argument("--metodo", choices=["ia", "optknock", "integral"], default="ia", help="Algoritmo")
    parser.add_argument("--ambiente", choices=["aerobio", "anaerobio"], default="aerobio", help="Condición de oxígeno (default: aerobio)")
    parser.add_argument("--generaciones", type=int, default=10)
    parser.add_argument("--poblacion", type=int, default=20)
    parser.add_argument("--mutaciones", type=int, default=4)
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument("--inyectar", action="store_true")
    parser.add_argument("-o", "--output", default="reportes")
    parser.add_argument("--no-reporte", action="store_true")
    parser.add_argument("--no-termino", action="store_true")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"], default="INFO")
    parser.add_argument("--version", action="version", version="MetaboEngineer Pro CLI v5.9")
    return parser

def setup_logging(log_level_str: str) -> str:
    os.makedirs("logs", exist_ok=True)
    log_file = f"logs/cli_{datetime.now():%Y%m%d_%H%M%S}.log"
    logging.basicConfig(level=getattr(logging, log_level_str), format="%(asctime)s [%(levelname)s] %(message)s", handlers=[logging.FileHandler(log_file, encoding="utf-8")])
    return log_file

def main():
    parser = build_parser()
    args = parser.parse_args()
    log_file = setup_logging(args.log_level)
    logging.info(f"Args: {vars(args)}")
    os.makedirs(args.output, exist_ok=True)

    if args.seed is None:
        args.seed = random.randint(0, 999999)
        print(f"[INFO] Seed generada: {args.seed}")
    random.seed(args.seed)

    print(f"[INFO] Cargando modelo '{args.modelo}'...")
    try:
        modelo, mensaje = cargar_modelo_universal(args.modelo)
        print(f"[OK] {mensaje} | {len(modelo.reactions)} reacciones")
    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr); return 2

    # APLICAR AMBIENTE
    print(f"[INFO] Configurando ambiente: {args.ambiente.upper()}...")
    modelo = configurar_ambiente(modelo, args.ambiente)

    if args.quimico not in modelo.reactions:
        if args.inyectar and args.quimico == "EX_caft_e":
            print(f"[INFO] Inyectando ruta heterologa para {args.quimico}...")
            modelo, args.quimico = inyectar_ruta_heterologa(modelo, args.quimico)
            if modelo is None: return 3
        else:
            print(f"[ERROR] '{args.quimico}' no existe. Usa --inyectar.", file=sys.stderr); return 3

    print(f"[INFO] Ejecutando '{args.metodo}' (seed={args.seed}, gen={args.generaciones}, pop={args.poblacion}, mut={args.mutaciones})...")

    try:
        if args.metodo == "ia":
            mejor = metodo_ia(modelo, args.quimico, args)
        elif args.metodo == "optknock":
            mejor = metodo_optknock(modelo, args.quimico, args)
        elif args.metodo == "integral":
            mejor = metodo_ia(modelo, args.quimico, args)
    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr); logging.error(traceback.format_exc()); return 4

    print(f"\n[OK] {args.metodo.upper()} completado")
    print(f"     Fitness: {mejor['fitness']:.4f} mmol/gDW/h")
    print(f"     Genes ({len(mejor['genes'])}): {', '.join(mejor['genes'])}")

    id_uniprot = get_organism_id(args.modelo)
    print(f"\n[INFO] Traduciendo genes (taxon={id_uniprot})...")
    for g in mejor["genes"]:
        nombre, funcion = traducir_gen_uniprot(g, id_uniprot)
        print(f"       {g:>15} | {nombre:<15} | {funcion[:50]}")

    if not args.no_termino:
        gen_biomasa = [r.id for r in modelo.reactions if "biomass" in r.id.lower()][0]
        try:
            print(f"\n[INFO] Analisis profundo + termodinamica...")
            analisis_profundo(modelo, gen_biomasa, mejor["genes"], args.quimico)
        except Exception as e:
            print(f"[WARN] Analisis profundo fallo: {e}")

    if not args.no_reporte and REPORTE_DISPONIBLE:
        try:
            print(f"\n[INFO] Generando reporte cientifico...")
            generar_documento_markdown(modelo, args.modelo, args.quimico, mejor["genes"], mejor["fitness"])
        except Exception as e:
            print(f"[WARN] Reporte fallo: {e}")

    try:
        parametros = {"generaciones": args.generaciones, "poblacion": args.poblacion, "mutaciones": args.mutaciones, "ambiente": args.ambiente, "cli": True}
        _exportar_resultado_json(args.modelo, args.quimico, args.metodo, mejor["genes"], mejor["fitness"], args.seed, parametros)
        print(f"\n[OK] Resultados exportados a '{args.output}/'")
    except Exception as e:
        print(f"[WARN] Exportacion fallo: {e}")

    print(f"\n[OK] Sesion completada. Log: {log_file}")
    return 0

# ==========================================
# METODOS DUPLICADOS PARA EL CLI
# ==========================================
from rich.progress import track
def metodo_ia(modelo, quimico, args):
    gen_biomasa = [r.id for r in modelo.reactions if "biomass" in r.id.lower()][0]
    modelo.reactions.get_by_id(quimico).bounds = (0, 1000.0)
    
    # SOLUCIÓN CIENTÍFICA: Solo buscar genes del metabolismo central (Glucólisis + TCA)
    # Esto evita que la IA apague transportes inútiles
    genes_metabolismo_central = []
    reacciones_centrales = ['PYK', 'PFK', 'FBA', 'GAPD', 'PGK', 'PGM', 'ENO', 'PDH', 'CS', 'ACONTa', 'ACONTb', 'ICDHyr', 'AKGDH', 'SUCOAS', 'SUCDi', 'FUM', 'MDH', 'ME1', 'ME2', 'PPC', 'PPS']
    
    for rxn_id in reacciones_centrales:
        if rxn_id in modelo.reactions:
            for gene in modelo.reactions.get_by_id(rxn_id).genes:
                if gene.id not in genes_metabolismo_central:
                    genes_metabolismo_central.append(gene.id)
                    
    if len(genes_metabolismo_central) < 5:
        # Fallback por si el modelo no tiene esos IDs exactos
        genes_metabolismo_central = [g.id for g in modelo.genes][:100] 
        
    print(f"[INFO] Restringiendo busqueda a {len(genes_metabolismo_central)} genes del metabolismo central...")
    todos_los_genes = genes_metabolismo_central

    poblacion = []
    for _ in range(args.poblacion):
        genes = random.sample(todos_los_genes, min(args.mutaciones, len(todos_los_genes)))
        poblacion.append({"genes": genes, "fitness": evaluar_bacteria(modelo, genes, quimico, gen_biomasa)})
        
    for gen in track(range(args.generaciones), description="Evolucionando..."):
        poblacion.sort(key=lambda x: x["fitness"], reverse=True)
        mejores = poblacion[:5]
        nueva = mejores.copy()
        while len(nueva) < args.poblacion:
            p1, p2 = random.choice(mejores)["genes"], random.choice(mejores)["genes"]
            hijo = list(set(p1[:3] + p2[3:]))
            if len(hijo) < args.mutaciones: hijo.append(random.choice(todos_los_genes))
            nueva.append({"genes": hijo, "fitness": evaluar_bacteria(modelo, hijo, quimico, gen_biomasa)})
        poblacion = nueva
    return sorted(poblacion, key=lambda x: x["fitness"], reverse=True)[0]

def metodo_optknock(modelo, quimico, args):
    gen_biomasa = [r.id for r in modelo.reactions if "biomass" in r.id.lower()][0]
    modelo.reactions.get_by_id(quimico).bounds = (0, 1000.0)
    sol_wt = modelo.optimize()
    activas = {r for r, v in sol_wt.fluxes.items() if abs(v) > 0.01 and not r.startswith("EX_")}
    candidatos = list({g.id for r in activas if r in modelo.reactions for g in modelo.reactions.get_by_id(r).genes})
    if len(candidatos) < 2: raise RuntimeError("Pocos genes activos.")
    def evaluar(genes):
        m = modelo.copy()
        for g in genes:
            if g in m.genes: m.genes.get_by_id(g).knock_out()
        m.objective = gen_biomasa
        s = m.optimize()
        return (s.fluxes.get(quimico, 0.0) / s.objective_value) if (s.status == "optimal" and s.objective_value > 0.01) else 0.0
    poblacion = []
    for _ in range(args.poblacion):
        genes = random.sample(candidatos, min(args.mutaciones, len(candidatos)))
        poblacion.append({"genes": genes, "fitness": evaluar(genes)})
    for gen in track(range(args.generaciones), description="Evolucionando..."):
        poblacion.sort(key=lambda x: x["fitness"], reverse=True)
        mejores = poblacion[:5]
        nueva = mejores.copy()
        while len(nueva) < args.poblacion:
            p1, p2 = random.choice(mejores)["genes"], random.choice(mejores)["genes"]
            hijo = list(set(p1[:2] + p2[2:]))
            if len(hijo) < args.mutaciones: hijo.append(random.choice(candidatos))
            nueva.append({"genes": hijo, "fitness": evaluar(hijo)})
        poblacion = nueva
    return sorted(poblacion, key=lambda x: x["fitness"], reverse=True)[0]

if __name__ == "__main__":
    try: sys.exit(main())
    except KeyboardInterrupt: sys.exit(130)