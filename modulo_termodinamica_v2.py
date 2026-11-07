import re
import warnings
warnings.filterwarnings("ignore")

from equilibrator_api import ComponentContribution, Q_
from rich.console import Console
from rich.table import Table

console = Console()

_SUFIJO_COMPARTIMENTO = re.compile(r"_[a-z]{1,2}$")


def _id_bigg_sin_compartimento(met_id):
    """'atp_c' -> 'atp'. eQuilibrator indexa BiGG sin el sufijo de compartimento."""
    return _SUFIJO_COMPARTIMENTO.sub("", met_id)


def construir_reaccion(modelo, reaccion):
    """Formula de eQuilibrator para una reaccion interna.

    Lanza ValueError si la reaccion es de frontera o tiene un solo lado: una reaccion asi
    no tiene un DG' bien definido y completarla con un sustrato inventado daria un numero falso.
    """
    if reaccion.boundary:
        raise ValueError(f"{reaccion.id} es una reaccion de frontera (intercambio/demanda/sumidero)")

    lado_izquierdo = []
    lado_derecho = []

    for metabolito, coeficiente in reaccion.metabolites.items():
        kegg_annot = metabolito.annotation.get('kegg.compound')
        if isinstance(kegg_annot, list): kegg_id = kegg_annot[0] if kegg_annot else None
        elif isinstance(kegg_annot, str): kegg_id = kegg_annot
        else: kegg_id = None

        if kegg_id and kegg_id.startswith('C'): termino = f"kegg:{kegg_id}"
        else: termino = f"bigg.metabolite:{_id_bigg_sin_compartimento(metabolito.id)}"

        if abs(coeficiente) != 1: 
            termino = f"{abs(coeficiente):g} {termino}"
            
        if coeficiente < 0: lado_izquierdo.append(termino)
        else: lado_derecho.append(termino)

    if not lado_izquierdo or not lado_derecho:
        raise ValueError(f"{reaccion.id} tiene un solo lado; no se puede calcular DG'")
    return f"{' + '.join(lado_izquierdo)} = {' + '.join(lado_derecho)}"


def analizar_termodinamica_universal(modelo, reacciones_ids, fluxes=None):
    """Analiza el DG' de las reacciones. Si se pasa `fluxes` (diccionario de flujos),
    tiene en cuenta el sentido de la reaccion para decidir si es viable."""
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

        formula = ""
        try:
            formula = construir_reaccion(modelo, reaccion)
            eq = cc.parse_reaction_formula(formula)
            dg = cc.physiological_dg_prime(eq).value.magnitude
            
            # FIX: Evaluar el sentido del flujo
            if fluxes is not None and r_id in fluxes:
                flujo = fluxes[r_id]
                if abs(flujo) < 1e-6:
                    viable = "N/A (Flujo 0)"
                    color = "dim"
                elif flujo > 0:
                    viable = "SI" if dg < 0 else "NO"
                    color = "green" if dg < 0 else "red"
                else: # flujo < 0 (reaccion inversa)
                    viable = "SI" if dg > 0 else "NO"
                    color = "green" if dg > 0 else "red"
            else:
                # Fallback si no se pasan flujos
                viable = "SI" if dg < 0 else "NO"
                color = "green" if dg < 0 else "red"

            formula_mostrar = formula[:40] + "..." if len(formula) > 40 else formula
            tabla.add_row(r_id, formula_mostrar, f"[{color}]{round(dg, 2)}[/{color}]", viable)
        except Exception as e:
            error_msg = str(e).split('\n')[0][:30]
            tabla.add_row(r_id, formula[:20], "Error API", error_msg)
    console.print(tabla)