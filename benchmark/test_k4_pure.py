# benchmark/test_k4_pure.py
import sys
from pathlib import Path
PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import random
import metabo_engineer as me
from rich.console import Console
from rich.panel import Panel

console = Console()

def main():
    modelo_id = "iML1515"
    quimico = "EX_succ_e"
    ambiente = "anaerobio"
    
    console.print(Panel.fit("[bold green]Experimento: GA Puro con k=4 (sin seeding)[/bold green]"))
    
    modelo, _ = me.cargar_modelo_universal(modelo_id)
    modelo = me.configurar_ambiente(modelo, ambiente)
    gen_biomasa = me._detectar_gen_biomasa(modelo)
    modelo.reactions.get_by_id(quimico).bounds = (0, 1000.0)
    
    todos_los_genes = me.obtener_genes_no_esenciales(modelo, gen_biomasa)
    
    # Semillas aleatorias para probar si converge al núcleo
    mejor_global = None
    for seed in [42, 100, 2024, 777, 1234]:
        console.print(f"\n[cyan]Corriendo GA Puro (seed={seed}, k=4)...[/cyan]")
        random.seed(seed)
        _cache = {}
        mejor = me._evolucionar(
            todos_los_genes, k=4, generaciones=15, pop_size=50,
            fitness_fn=lambda g: me.evaluar_bacteria(modelo, g, quimico, gen_biomasa, _cache=_cache),
            crossover_split=2
        )
        console.print(f"  Fitness: {mejor['fitness']:.4f} | Genes: {', '.join(mejor['genes'])}")
        
        if mejor_global is None or mejor["fitness"] > mejor_global["fitness"]:
            mejor_global = mejor

    console.print(f"\n[bold green]Mejor resultado del GA Puro k=4:[/bold green]")
    console.print(f"  Fitness: {mejor_global['fitness']:.4f}")
    console.print(f"  Genes: {', '.join(mejor_global['genes'])}")
    
    # Comparamos con el nucleo destilado
    nucleo = ["b1380", "b1241", "b1478", "b3738"]
    fit_nucleo = me.evaluar_bacteria(modelo, nucleo, quimico, gen_biomasa)
    console.print(f"\n[bold yellow]Núcleo teórico (4 genes):[/bold yellow]")
    console.print(f"  Fitness: {fit_nucleo:.4f}")
    console.print(f"  Genes: {', '.join(nucleo)}")

if __name__ == "__main__":
    main()
