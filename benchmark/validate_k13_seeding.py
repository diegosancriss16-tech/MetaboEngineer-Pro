# benchmark/validate_k13_seeding.py
import sys
import argparse
import random
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import metabo_engineer as me
from rich.console import Console
from rich.panel import Panel

console = Console()

def resolver_genes(modelo, nombres):
    """Resuelve nombres comunes (ldhA) o IDs de BiGG (b1380) a IDs reales del modelo."""
    resueltos = []
    for nombre in nombres:
        # 1. ¿Existe como gene.id directo?
        if nombre in modelo.genes:
            resueltos.append(nombre)
            continue
        # 2. ¿Existe como gene.name? (búsqueda por nombre)
        encontrado = False
        for g in modelo.genes:
            if g.name == nombre or g.name.lower() == nombre.lower():
                resueltos.append(g.id)
                encontrado = True
                break
        if not encontrado:
            console.print(f"[yellow]WARN[/yellow] No se encontró '{nombre}' ni como ID ni como nombre.")
    return resueltos

def main():
    parser = argparse.ArgumentParser(description="Valida si el GA cruza el valle con seeding.")
    parser.add_argument("--modelo", default="iML1515")
    parser.add_argument("--quimico", default="EX_succ_e")
    parser.add_argument("--ambiente", default="anaerobio", choices=["aerobio", "anaerobio"])
    parser.add_argument("--generaciones", type=int, default=20)
    parser.add_argument("--poblacion", type=int, default=50)
    parser.add_argument("--seed", type=int, default=100)
    args = parser.parse_args()

    console.print(Panel.fit("[bold green]Experimento de Seeding (KJ060 extendido)[/bold green]"))
    
    try:
        modelo, msg = me.cargar_modelo_universal(args.modelo)
        console.print(f"[green]OK[/green] {msg}")
    except Exception as e:
        console.print(f"[red]ERROR[/red] {e}")
        return

    modelo = me.configurar_ambiente(modelo, args.ambiente)
    gen_biomasa = me._detectar_gen_biomasa(modelo)
    
    if args.quimico not in modelo.reactions:
        console.print(f"[red]ERROR[/red] El químico '{args.quimico}' no existe.")
        return
    modelo.reactions.get_by_id(args.quimico).bounds = (0, 1000.0)

    # Semilla con nombres comunes (el script los traducirá a IDs de BiGG)
    semilla_nombres = [
        "ldhA", "adhE", "adhP", "mhpF", "ackA", "tdcD", "purT", 
        "focA", "focB", "pflB", "yfiD", "alt", "atpE"
    ]
    
    console.print(f"[cyan]Resolviendo {len(semilla_nombres)} nombres a IDs de BiGG...[/cyan]")
    semilla_real = resolver_genes(modelo, semilla_nombres)
    
    if not semilla_real:
        console.print("[red]ERROR[/red] Ningún gen de la semilla se pudo resolver. Abortando.")
        return
    
    console.print(f"[green]OK[/green] Semilla resuelta ({len(semilla_real)} genes): {', '.join(semilla_real)}")
    
    todos_los_genes = me.obtener_genes_no_esenciales(modelo, gen_biomasa)
    _cache_eval = {}
    
    def fitness_fn(g):
        return me.evaluar_bacteria(modelo, g, args.quimico, gen_biomasa, _cache=_cache_eval)

    random.seed(args.seed)
    mejor = me._evolucionar(
        pool_genes=todos_los_genes, 
        k=len(semilla_real), 
        generaciones=args.generaciones, 
        pop_size=args.poblacion, 
        fitness_fn=fitness_fn, 
        crossover_split=3, 
        descripcion="Evolucionando con semilla...",
        semillas=[semilla_real]
    )
    
    console.print(f"\n[green]OK[/green] Evolución completada.")
    console.print(f"  Fitness final: {mejor['fitness']:.4f}")
    console.print(f"  Genes finales ({len(mejor['genes'])}): {', '.join(mejor['genes'])}")
    
    me._exportar_resultado_json(
        args.modelo, args.quimico, "ia_seeded", mejor["genes"], mejor["fitness"], args.seed, 
        {"gens": args.generaciones, "pop": args.poblacion, "seeded": True}
    )
    
    console.print("\n[cyan]Calculando retención @90% del ganador...[/cyan]")
    me.analisis_profundo(modelo, gen_biomasa, mejor["genes"], args.quimico)

if __name__ == "__main__":
    main()
