# benchmark/distill_design.py
import sys
import json
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import metabo_engineer as me
from rich.console import Console
from rich.table import Table

console = Console()

def main():
    modelo_id = "iML1515"
    quimico = "EX_succ_e"
    ambiente = "anaerobio"
    
    # Los genes exactos que reportó tu GA
    genes_ga = ["b1380", "b1241", "b1478", "b4287", "b3738", "b3583", "b0905", "b1622", "b4086", "b0355", "b4195"]
    
    modelo, _ = me.cargar_modelo_universal(modelo_id)
    modelo = me.configurar_ambiente(modelo, ambiente)
    gen_biomasa = me._detectar_gen_biomasa(modelo)
    modelo.reactions.get_by_id(quimico).bounds = (0, 1000.0)

    console.print(f"[cyan]Evaluando diseño original ({len(genes_ga)} genes)...[/cyan]")
    fit_orig = me.evaluar_bacteria(modelo, genes_ga, quimico, gen_biomasa)
    
    tabla = Table(title="Destilación Leave-One-Out del Diseño GA")
    tabla.add_column("Gen Eliminado", style="cyan")
    tabla.add_column("Fitness Resultante", style="green")
    tabla.add_column("Delta", style="magenta")
    
    mejor_fit = fit_orig
    mejor_genes = genes_ga.copy()
    
    tabla.add_row("[bold]Ninguno (Original)[/bold]", f"{fit_orig:.4f}", "-")
    
    for gen in genes_ga:
        subset = [g for g in genes_ga if g != gen]
        fit = me.evaluar_bacteria(modelo, subset, quimico, gen_biomasa)
        delta = fit - fit_orig
        
        color = "green" if delta > 0 else "red"
        tabla.add_row(gen, f"{fit:.4f}", f"[{color}]{delta:+.4f}[/{color}]")
        
        if fit > mejor_fit:
            mejor_fit = fit
            mejor_genes = subset
            
    console.print(tabla)
    console.print(f"\n[bold green]¡El mejor diseño destilado tiene {len(mejor_genes)} genes y un Fitness de {mejor_fit:.4f}![/bold green]")
    console.print(f"Genes: {', '.join(mejor_genes)}")
    
    # Guardar en JSON
    out_data = {
        "fitness": mejor_fit,
        "genes": mejor_genes,
        "original_fitness": fit_orig,
        "original_genes": genes_ga
    }
    out_path = PROJECT_ROOT / "benchmark" / "results" / "distilled_best_design.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(out_data, f, indent=2)
    console.print(f"[green]OK[/green] Guardado en {out_path}")

if __name__ == "__main__":
    main()