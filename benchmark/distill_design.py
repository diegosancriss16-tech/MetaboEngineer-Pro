# benchmark/distill_design.py
import sys
import json
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import metabo_engineer as me
from rich.console import Console

console = Console()

def main():
    modelo_id = "iML1515"
    quimico = "EX_succ_e"
    ambiente = "anaerobio"
    
    # Los genes que el GA puro (semilla 2024) encontró
    genes_ga = ["b1090", "b2415", "b3732", "b3708"]
    
    modelo, _ = me.cargar_modelo_universal(modelo_id)
    modelo = me.configurar_ambiente(modelo, ambiente)
    gen_biomasa = me._detectar_gen_biomasa(modelo)
    modelo.reactions.get_by_id(quimico).bounds = (0, 1000.0)

    console.print(f"[cyan]1. Evaluando diseño original del GA Puro ({len(genes_ga)} genes)...[/cyan]")
    fit_orig = me.evaluar_bacteria(modelo, genes_ga, quimico, gen_biomasa)
    console.print(f"   Fitness original: {fit_orig:.4f}")
    
    current_genes = list(genes_ga)
    current_fit = fit_orig
    
    console.print("\n[cyan]2. Aplicando Destilación Iterativa (Leave-One-Out)...[/cyan]")
    
    # CORREGIDO: Lógica de destilación segura (sin sobrescribir mejoras con neutros)
    while len(current_genes) > 1:
        cand = [(me.evaluar_bacteria(modelo, [g for g in current_genes if g != x], quimico, gen_biomasa), x) for x in current_genes]
        fit, gen = max(cand)
        if fit < current_fit - 1e-4:
            break
        console.print(f"   -> Eliminado gen: {gen} (Fitness: {fit:.4f})")
        current_genes.remove(gen)
        current_fit = fit
            
    console.print(f"\n[bold green]¡El diseño destilado (núcleo) tiene {len(current_genes)} genes y un Fitness de {current_fit:.4f}![/bold green]")
    console.print(f"Genes: {', '.join(current_genes)}")
    
    # Guardar JSON en reportes/
    out_dir = PROJECT_ROOT / "reportes"
    out_dir.mkdir(parents=True, exist_ok=True)
    
    with open(out_dir / "design_ga_pure.json", "w") as f:
        json.dump({"fitness": fit_orig, "genes": genes_ga}, f, indent=2)
        
    with open(out_dir / "design_distilled_core.json", "w") as f:
        json.dump({"fitness": current_fit, "genes": current_genes}, f, indent=2)

    console.print(f"[green]OK[/green] JSONs guardados en reportes/")

if __name__ == "__main__":
    main()
