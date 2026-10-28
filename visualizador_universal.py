import warnings
warnings.filterwarnings("ignore")

import os
import networkx as nx
from pyvis.network import Network
import cobra
from cobra.io import load_model
from rich.console import Console

console = Console()


def generar_mapa_universal(modelo_id="iYO844", quimico="EX_lys__L_e"):
    console.print(f"\n[cyan]INFO[/cyan] Generando Mapa Metabolico para {modelo_id}...\n")

    ruta_local = f"modelos/{modelo_id}.xml"
    try:
        if os.path.exists(ruta_local):
            modelo = cobra.io.read_sbml_model(ruta_local)
        else:
            modelo = load_model(modelo_id)
    except Exception as e:
        console.print(f"[red]ERROR[/red] No se pudo cargar el modelo: {e}")
        return

    if quimico not in modelo.reactions:
        console.print(f"[red]ERROR[/red] Quimico {quimico} no valido en este modelo.")
        return

    modelo.reactions.get_by_id(quimico).bounds = (0, 1000.0)
    modelo.objective = quimico
    solucion = modelo.optimize()

    G = nx.DiGraph()
    metabolitos_ignorados = {
        'h_c', 'h_e', 'h2o_c', 'h2o_e', 'atp_c', 'adp_c', 'amp_c',
        'pi_c', 'ppi_c', 'nad_c', 'nadh_c', 'nadp_c', 'nadph_c', 'co2_c', 'co2_e'
    }

    for reaccion in modelo.reactions:
        flujo = solucion.fluxes[reaccion.id]
        if abs(flujo) > 1.0 and not reaccion.id.startswith("EX_"):
            reactivos = [m.id for m, coef in reaccion.metabolites.items()
                         if coef < 0 and m.id not in metabolitos_ignorados]
            productos = [m.id for m, coef in reaccion.metabolites.items()
                         if coef > 0 and m.id not in metabolitos_ignorados]

            for r in reactivos:
                for p in productos:
                    color = "blue" if flujo > 0 else "red"
                    G.add_edge(r, p, title=f"{reaccion.id} (Flux: {flujo:.1f})", color=color)

    if len(G.nodes) > 0:
        largest_cc = max(nx.weakly_connected_components(G), key=len)
        G = G.subgraph(largest_cc)

    net = Network(height="800px", width="100%", bgcolor="#1e1e1e", font_color="white", directed=True)
    net.from_nx(G)
    net.repulsion(node_distance=150, spring_length=200)

    if not os.path.exists("reportes"):
        os.makedirs("reportes")

    ruta_html = f"reportes/mapa_global_{modelo_id}_{quimico}.html"
    net.save_graph(ruta_html)

    console.print(f"[green]OK[/green] Mapa guardado en: [cyan]{ruta_html}[/cyan]")