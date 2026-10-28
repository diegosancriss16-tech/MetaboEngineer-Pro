import warnings
warnings.filterwarnings("ignore")

import re
from cobra.io import read_sbml_model, load_model
from equilibrator_api import ComponentContribution, Q_
from rich.console import Console
from rich.table import Table

console = Console()


def construir_reaccion(modelo, reaccion):
    """Construye la formula usando anotaciones de KEGG o IDs de BiGG."""
    lado_izquierdo = []
    lado_derecho = []

    for metabolito, coeficiente in reaccion.metabolites.items():
        kegg_annot = metabolito.annotation.get('kegg.compound')
        if isinstance(kegg_annot, list):
            kegg_id = kegg_annot[0] if kegg_annot else None
        elif isinstance(kegg_annot, str):
            kegg_id = kegg_annot
        else:
            kegg_id = None

        if kegg_id and kegg_id.startswith('C'):
            termino = f"kegg:{kegg_id}"
        else:
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
    console.print(f"\n[cyan]INFO[/cyan] Termodinamica Universal (BiGG + KEGG + eQuilibrator)\n")

    cc = ComponentContribution()
    cc.p_h = Q_(7.0)
    cc.p_mg = Q_(3.0)
    cc.ionic_strength = Q_("0.1 M")

    tabla = Table(title="Termodinamica Automatica")
    tabla.add_column("Reaccion", style="cyan")
    tabla.add_column("Formula de Entrada", style="dim")
    tabla.add_column("DG' (kJ/mol)", style="magenta")
    tabla.add_column("Viable", style="green")

    for r_id in reacciones_ids:
        if r_id not in modelo.reactions:
            continue

        reaccion = modelo.reactions.get_by_id(r_id)
        formula = construir_reaccion(modelo, reaccion)

        try:
            eq = cc.parse_reaction_formula(formula)
            dg = cc.physiological_dg_prime(eq).value.magnitude
            viable = "SI" if dg < 0 else "NO"
            color = "green" if dg < 0 else "red"
            formula_mostrar = formula[:40] + "..." if len(formula) > 40 else formula
            tabla.add_row(r_id, formula_mostrar, f"[{color}]{round(dg, 2)}[/{color}]", viable)
        except Exception as e:
            error_msg = str(e).split('\n')[0][:30]
            tabla.add_row(r_id, formula[:20], "Error", error_msg)

    console.print(tabla)


if __name__ == "__main__":
    console.print("[cyan]INFO[/cyan] Probando con iMM904...")
    try:
        modelo = read_sbml_model("modelos/iMM904.xml")
        analizar_termodinamica_universal(modelo, ['PGI', 'PFK', 'FBA', 'PGK', 'ENO', 'PYK'])
    except Exception as e:
        console.print(f"[red]ERROR[/red] {e}")