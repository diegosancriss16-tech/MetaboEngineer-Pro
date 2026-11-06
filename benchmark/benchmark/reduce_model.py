import sys
from pathlib import Path
PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
from cobra.io import read_sbml_model, write_sbml_model
from cobra.flux_analysis import fastcc
import metabo_engineer as me

print("Cargando iML1515...")
model = read_sbml_model("iML1515.xml")
model.solver = 'glpk'

# Aplicar el mismo medio que usas en tus benchmarks
model = me.configurar_ambiente(model, "anaerobio")
model.solver = 'glpk'

print(f"Reacciones originales: {len(model.reactions)}")
print("Ejecutando fastCC (esto puede tardar unos minutos)...")

# fastCC devuelve un modelo consistente
reduced_model = fastcc(model)

print(f"Reacciones después de fastCC: {len(reduced_model.reactions)}")
write_sbml_model(reduced_model, "iML1515_reduced.xml")
print("✅ Modelo reducido guardado como iML1515_reduced.xml")