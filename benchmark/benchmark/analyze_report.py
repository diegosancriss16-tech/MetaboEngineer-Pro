# benchmark/analyze_report.py
import sys
import argparse
import json
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import metabo_engineer as me
from rich.console import Console

console = Console()

def analizar_reporte(ruta_json):
    with open(ruta_json, "r", encoding="utf-8") as f:
        data = json.load(f)
        
    modelo_id = data.get("modelo")
    quimico = data.get("quimico")
    genes = data.get("genes_knockout", [])
    
    console.print(f"[cyan]Analizando reporte:[/cyan] {ruta_json.name}")
    console.print(f"  Modelo: {modelo_id}")
    console.print(f"  Quimico: {quimico}")
    console.print(f"  Genes ({len(genes)}): {', '.join(genes)}")
    
    if not genes:
        console.print("[red]ERROR[/red] El reporte no contiene genes.")
        return

    # Cargar modelo
    try:
        modelo, msg = me.cargar_modelo_universal(modelo_id)
        console.print(f"[green]OK[/green] {msg}")
    except Exception as e:
        console.print(f"[red]ERROR[/red] No se pudo cargar el modelo: {e}")
        return

    # Detectar ambiente desde el reporte (si existe) o asumir anaerobio para succinato
    ambiente = data.get("ambiente", "anaerobio")
    modelo = me.configurar_ambiente(modelo, ambiente)
    
    gen_biomasa = me._detectar_gen_biomasa(modelo)
    
    # Asegurar que el químico esté abierto
    if quimico in modelo.reactions:
        modelo.reactions.get_by_id(quimico).bounds = (0, 1000.0)
    else:
        console.print(f"[red]ERROR[/red] El químico '{quimico}' no existe en el modelo.")
        return

    console.print(f"\n[cyan]Calculando retención @90% en {ambiente}...[/cyan]")
    me.analisis_profundo(modelo, gen_biomasa, genes, quimico)

def main():
    parser = argparse.ArgumentParser(description="Analiza la retención de un reporte JSON.")
    parser.add_argument("--reporte", required=True, help="Ruta al archivo JSON del reporte.")
    args = parser.parse_args()
    
    ruta = Path(args.reporte)
    if not ruta.exists():
        console.print(f"[red]ERROR[/red] El archivo '{ruta}' no existe.")
        return
        
    analizar_reporte(ruta)

if __name__ == "__main__":
    main()