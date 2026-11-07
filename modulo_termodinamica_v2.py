import warnings
warnings.filterwarnings("ignore")

from equilibrator_api import ComponentContribution, Q_
from rich.console import Console
from rich.table import Table

console = Console()

def construir_reaccion(modelo, reaccion):
    lado_izquierdo = []
    lado_derecho = []

    for metabolito, coeficiente in reaccion.metabolites.items():
        kegg_annot = metabolito.annotation.get('kegg.compound')
        if isinstance(kegg_annot, list): kegg_id = kegg_annot[0] if kegg_annot else None
        elif isinstance(kegg_annot, str): kegg_id = kegg_annot
        else: kegg_id = None

        if kegg_id and kegg_id.startswith('C'): termino = f"kegg:{kegg_id}"
        else: termino = f"bigg:{metabolito.id}"

        # CORREGIDO Bug 6: No truncar con int(). Usar :g para decimales (ej. 0.5)
        if abs(coeficiente) != 1: 
            termino = f"{abs(coeficiente):g} {termino}"
            
        if coeficiente < 0: lado_izquierdo.append(termino)
        else: lado_derecho.append(termino)

    izq_str = " + ".join(lado_izquierdo) if lado_izquierdo else "bigg:h2o_c"
    der_str = " + ".join(lado_derecho) if lado_derecho else "bigg:h2o_c"
    return f"{izq_str} = {der_str}"

def analizar_termodinamica_universal(modelo, reacciones_ids):
    console.print("\n[cyan]INFO[/cyan] Termodinamica Universal\n")
    cc = ComponentContribution()
    cc.p_h, cc.p_mg, cc.ionic_strength = Q_(7.0), Q_(3.0), Q_("0.1 M")

    tabla = Table(title="Termodinamica Automatica")
    tabla.add_column("Reaccion", style="cyan"); tabla.add_column("Formula", style="dim"); tabla.add_column("DG' (kJ/mol)", style="magenta"); tabla.add_column("Viable", style="green")

    for r_id in reacciones_ids:
        if r_id not in modelo.reactions: continue
        
        if r_id.startswith("EX_"):
            tabla.add_row(r_id, "Intercambio", "N/A", "N/A")
            continue
            
        reaccion = modelo.reactions.get_by_id(r_id)
        if r_id.startswith("CAFT") or r_id.startswith("XAN"):
            tabla.add_row(r_id, "Metabolito Sintetico", "N/A", "N/A")
            continue

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
            tabla.add_row(r_id, formula[:20], "Error API", error_msg)
    console.print(tabla)