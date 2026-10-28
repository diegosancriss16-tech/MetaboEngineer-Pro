import warnings
warnings.filterwarnings("ignore")

import cobra
from cobra.io import load_model
from rich.console import Console
from rich.table import Table
from rich.progress import track

console = Console()

def buscar_id_gen(modelo, nombre_comun):
    """Busca el ID real de un gen basándose en su nombre común o parte de su ID."""
    for gen in modelo.genes:
        # Busca en el ID (ej. b3919) o en el nombre completo (ej. tpiA)
        if nombre_comun.lower() in gen.id.lower() or (gen.name and nombre_comun.lower() in gen.name.lower()):
            return gen.id
    return None

def probar_letalidad_sintetica():
    console.print("🧬 Iniciando Test de Letalidad Sintética (Doble Knockout)...\n", style="bold cyan")
    
    # 1. Cargar el modelo de E. coli
    modelo = load_model("textbook")
    console.print(f"✅ Modelo cargado: {modelo.id} ({len(modelo.genes)} genes)")
    
    # 2. Calcular el crecimiento normal (Wild-Type)
    solucion_wt = modelo.optimize()
    crecimiento_wt = solucion_wt.objective_value
    console.print(f"📈 Crecimiento Wild-Type: {crecimiento_wt:.4f} h⁻¹\n")
    
    # 3. Pares de genes a probar (usando nombres comunes)
    pares_a_probar = [
        ('tpi', 'fba'),   # Letal conocido
        ('tpi', 'gap'),   # Letal conocido
        ('pgi', 'zwf'),   # Bloqueo de pentosas
        ('ppc', 'pck'),   # Ruta anaplerótica
        ('sdh', 'frd'),   # Respiración vs Fermentación
        ('fba', 'gap')    # Control
    ]
    
    resultados = []
    
    console.print("🔬 Buscando IDs y probando dobles mutaciones...\n")
    
    for par in track(pares_a_probar, description="Evaluando pares..."):
        nombre1, nombre2 = par
        
        # Buscar los IDs reales en el modelo
        gen1 = buscar_id_gen(modelo, nombre1)
        gen2 = buscar_id_gen(modelo, nombre2)
        
        if gen1 and gen2:
            # Clonar el modelo
            modelo_mutante = modelo.copy()
            
            # Apagar ambos genes
            modelo_mutante.genes.get_by_id(gen1).knock_out()
            modelo_mutante.genes.get_by_id(gen2).knock_out()
            
            # Calcular nuevo crecimiento
            solucion_mut = modelo_mutante.optimize()
            crecimiento_mut = solucion_mut.objective_value if solucion_mut.status == 'optimal' else 0.0
                
            # Calcular la reducción porcentual
            reduccion = ((crecimiento_wt - crecimiento_mut) / crecimiento_wt) * 100 if crecimiento_wt > 0 else 0
            
            # Determinar si es letal sintético (crecimiento casi 0)
            es_letal = "☠️ SÍ" if crecimiento_mut < 0.01 else "✅ NO"
            
            resultados.append({
                "Gen 1": f"{nombre1} ({gen1})",
                "Gen 2": f"{nombre2} ({gen2})",
                "Crecimiento": round(crecimiento_mut, 4),
                "Reducción (%)": round(reduccion, 1),
                "¿Letal?": es_letal
            })
        else:
            console.print(f"⚠️ No se encontró uno de los genes: {nombre1} o {nombre2}", style="yellow")

    # 4. Mostrar resultados en tabla
    if resultados:
        tabla = Table(title="Resultados del Test de Letalidad Sintética")
        tabla.add_column("Gen 1", style="cyan")
        tabla.add_column("Gen 2", style="cyan")
        tabla.add_column("Crecimiento", style="magenta")
        tabla.add_column("Reducción", style="yellow")
        tabla.add_column("¿Letal?", style="green")
        
        for r in resultados:
            tabla.add_row(r["Gen 1"], r["Gen 2"], str(r["Crecimiento"]), f"{r['Reducción (%)']}%", r["¿Letal?"])
            
        console.print(tabla)
        console.print("\n💡 Nota: 'Letal' significa que la célula no puede crecer si le faltan ambos genes al mismo tiempo.")

if __name__ == "__main__":
    probar_letalidad_sintetica()