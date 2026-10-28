import warnings
warnings.filterwarnings("ignore")

import re
from cobra.io import read_sbml_model, load_model
from equilibrator_api import ComponentContribution, Q_
from rich.console import Console
from rich.table import Table

console = Console()

def construir_reaccion(modelo, reaccion):
    """Construye la fórmula usando anotaciones de KEGG o IDs de BiGG."""
    lado_izquierdo = []
    lado_derecho = []
    
    for metabolito, coeficiente in reaccion.metabolites.items():
        # 1. Obtener el ID de KEGG de forma segura (verificando si es lista o string)
        kegg_annot = metabolito.annotation.get('kegg.compound')
        
        if isinstance(kegg_annot, list):
            kegg_id = kegg_annot[0] if kegg_annot else None
        elif isinstance(kegg_annot, str):
            kegg_id = kegg_annot
        else:
            kegg_id = None
            
        # 2. Formatear el término. Si hay KEGG usarlo, si no, usar BiGG directo.
        if kegg_id and kegg_id.startswith('C'):
            termino = f"kegg:{kegg_id}"
        else:
            # eQuilibrator sabe leer IDs de BiGG directamente si están en su caché
            termino = f"bigg:{metabolito.id}"
            
        if abs(coeficiente) > 1:
            termino = f"{int(abs(coeficiente))} {termino}"
            
        if coeficiente < 0:
            lado_izquierdo.append(termino)
        else:
            lado_derecho.append(termino)
            
    izq_str = " + ".join(lado_izquierdo) if lado_izquierdo else "bigg:h2o_c"
    der_str = " + ".join(lado_derecho) if lado_derecho else "bigg:h2o_c"
    
    return f"{izq_str} = {der_str}"

def analizar_termodinamica_universal(modelo, reacciones_ids):
    console.print(f"\n[bold cyan]🌡️ Iniciando Termodinámica Universal (V2.2 - Híbrida)...[/bold cyan]")
    console.print("Usando base de datos local (BiGG + KEGG). Esto será instantáneo.\n")
    
    cc = ComponentContribution()
    cc.p_h = Q_(7.0)
    cc.p_mg = Q_(3.0)
    cc.ionic_strength = Q_("0.1 M")
    
    tabla = Table(title="Termodinámica Automática (Cualquier Reacción)")
    tabla.add_column("Reacción", style="cyan")
    tabla.add_column("Fórmula de Entrada", style="dim")
    tabla.add_column("ΔG' (kJ/mol)", style="magenta")
    tabla.add_column("Viabilidad", style="green")
    
    for r_id in reacciones_ids:
        if r_id not in modelo.reactions:
            continue
            
        reaccion = modelo.reactions.get_by_id(r_id)
        formula = construir_reaccion(modelo, reaccion)
        
        try:
            eq = cc.parse_reaction_formula(formula)
            dg = cc.physiological_dg_prime(eq).value.magnitude
            viable = "✅ SÍ" if dg < 0 else "❌ NO"
            
            color = "green" if dg < 0 else "red"
            # Acortamos la fórmula para que se vea bonita en la tabla
            formula_mostrar = formula[:40] + "..." if len(formula) > 40 else formula
            
            tabla.add_row(
                r_id, 
                formula_mostrar, 
                f"[{color}]{round(dg, 2)}[/{color}]", 
                viable
            )
        except Exception as e:
            # Si falla, mostramos el error específico
            error_msg = str(e).split('\n')[0][:30]
            tabla.add_row(r_id, formula[:20], "Error", error_msg)
            
    console.print(tabla)
    console.print("\n💡 [bold]Ventaja V2.2:[/bold] Ahora usa los IDs de BiGG directamente con eQuilibrator, evitando errores de traducción.")

if __name__ == "__main__":
    console.print("Cargando modelo iMM904...")
    try:
        # Cargar el modelo de levadura
        modelo = read_sbml_model("modelos/iMM904.xml")
        
        # Probemos reacciones clave
        reacciones_a_probar = ['PGI', 'PFK', 'FBA', 'PGK', 'ENO', 'PYK', 'PDC1', 'ADH1', 'MET1']
        
        analizar_termodinamica_universal(modelo, reacciones_a_probar)
        
    except Exception as e:
        console.print(f"[red]Error cargando el modelo: {e}[/red]")