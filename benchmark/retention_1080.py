# benchmark/retention_1080.py
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import json
import glob
import metabo_engineer as me
from rich.console import Console

console = Console()

# 1. Cargar iML1515
try:
    modelo, msg = me.cargar_modelo_universal("iML1515")
    console.print(f"[green]OK[/green] {msg}")
except Exception as e:
    console.print(f"[red]ERROR[/red] {e}")
    raise SystemExit(1)

# 2. Configurar anaerobio
modelo = me.configurar_ambiente(modelo, "anaerobio")
gen_biomasa = me._detectar_gen_biomasa(modelo)
quimico = "EX_succ_e"
modelo.reactions.get_by_id(quimico).bounds = (0, 1000.0)

# 3. Buscar SOLO reportes de iML1515 + EX_succ_e (cualquier método)
patrones = [
    "ia_iML1515_EX_succ_e_*.json",
    "integral_iML1515_EX_succ_e_*.json",
    "optknock_iML1515_EX_succ_e_*.json",
    "optknock_exacto_iML1515_EX_succ_e_*.json",
]

candidatos = []
for patron in patrones:
    for ruta in glob.glob(str(PROJECT_ROOT / "reportes" / patron)):
        try:
            with open(ruta, encoding="utf-8") as f:
                d = json.load(f)
            fit = d.get("fitness", 0)
            # Rango amplio: 9-12 captura tanto 10.1242 como 10.80
            if 9.0 <= fit <= 12.0:
                candidatos.append((fit, d.get("fecha", ""), ruta, d.get("genes_knockout", [])))
        except Exception:
            pass

if not candidatos:
    console.print("[yellow]No encontré reporte de iML1515/EX_succ_e con fitness entre 9 y 12.[/yellow]")
    console.print()
    console.print("[cyan]Reportes disponibles en total:[/cyan]")
    for ruta in sorted(glob.glob(str(PROJECT_ROOT / "reportes" / "*.json"))):
        try:
            with open(ruta, encoding="utf-8") as f:
                d = json.load(f)
            console.print(f"  {Path(ruta).name}:")
            console.print(f"    fit={d.get('fitness')}, modelo={d.get('modelo')}, quimico={d.get('quimico')}, metodo={d.get('metodo')}")
        except Exception:
            pass
    console.print()
    console.print("[yellow]Para generar el reporte de 10.80, corre el experimento de seeding:[/yellow]")
    console.print("  python benchmark\\validate_k13_seeding.py")
    console.print()
    console.print("[cyan]Alternativa manual:[/cyan] Edita este script y reemplaza la búsqueda por una lista de genes hardcoded.")
    raise SystemExit(0)

candidatos.sort(reverse=True)
fit, fecha, ruta, genes = candidatos[0]
console.print(f"[cyan]Usando el mejor reporte encontrado:[/cyan] {ruta}")
console.print(f"  Fitness: {fit}")
console.print(f"  Fecha:   {fecha}")
console.print(f"  Genes ({len(genes)}): {', '.join(genes)}")

# 4. Validar genes en el modelo
genes_validos = [g for g in genes if g in modelo.genes]
if len(genes_validos) != len(genes):
    invalidos = set(genes) - set(genes_validos)
    console.print(f"[yellow]WARN[/yellow] Genes no en modelo: {invalidos}")
    genes = genes_validos

if not genes:
    console.print("[red]ERROR[/red] Ningún gen válido. Abortando.")
    raise SystemExit(1)

# 5. Calcular retención @90%
console.print(f"\n[cyan]Calculando retención @90% para {len(genes)} genes...[/cyan]")
me.analisis_profundo(modelo, gen_biomasa, genes, quimico)