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
    
    # Los 11 genes que el GA encontró originalmente
    genes_ga = ["b1380", "b1241", "b1478", "b4287", "b3738", "b3583", "b0905", "b1622", "b4086", "b0355", "b4195"]
    
    modelo, _ = me.cargar_modelo_universal(modelo_id)
    modelo = me.configurar_ambiente(modelo, ambiente)
    gen_biomasa = me._detectar_gen_biomasa(modelo)
    modelo.reactions.get_by_id(quimico).bounds = (0, 1000.0)

    console.print(f"[cyan]1. Evaluando diseño original del GA ({len(genes_ga)} genes)...[/cyan]")
    fit_orig = me.evaluar_bacteria(modelo, genes_ga, quimico, gen_biomasa)
    console.print(f"   Fitness original: {fit_orig:.4f}")
    
    # 2. Destilación Iterativa
    current_genes = list(genes_ga)
    current_fit = fit_orig
    improved = True
    
    console.print("\n[cyan]2. Aplicando Destilación Iterativa (Leave-One-Out)...[/cyan]")
    
    while improved:
        improved = False
        best_drop = None
        best_drop_fit = current_fit
        
        for gen in current_genes:
            subset = [g for g in current_genes if g != gen]
            fit = me.evaluar_bacteria(modelo, subset, quimico, gen_biomasa)
            # Si el fitness no baja, nos podemos deshacer del gen
            if fit >= current_fit - 1e-4:
                if fit > best_drop_fit:
                    best_drop = gen
                    best_drop_fit = fit
                elif fit == current_fit:
                    best_drop = gen # Eliminar pasajero
                    
        if best_drop:
            console.print(f"   -> Eliminado gen: {best_drop} (Fitness: {best_drop_fit:.4f})")
            current_genes.remove(best_drop)
            current_fit = best_drop_fit
            improved = True
            
    console.print(f"\n[bold green]¡El diseño destilado (núcleo) tiene {len(current_genes)} genes y un Fitness de {current_fit:.4f}![/bold green]")
    console.print(f"Genes: {', '.join(current_genes)}")
    
    # 3. Guardar JSONs (en la raíz para evitar .gitignore)
    out_dir = PROJECT_ROOT / "reportes"
    out_dir.mkdir(parents=True, exist_ok=True)
    
    with open(out_dir / "design_ga_11genes.json", "w") as f:
        json.dump({"fitness": fit_orig, "genes": genes_ga}, f, indent=2)
        
    with open(out_dir / "design_distilled_core.json", "w") as f:
        json.dump({"fitness": current_fit, "genes": current_genes}, f, indent=2)

    console.print(f"[green]OK[/green] JSONs guardados en reportes/")

if __name__ == "__main__":
    main()