from cobra.io import load_model, write_sbml_model

print("Descargando iJN1463 desde la API de BiGG...")
try:
    # COBRApy maneja las redirecciones y la API por nosotros
    modelo = load_model("iJN1463")
    write_sbml_model(modelo, "iJN1463.xml")
    print("✅ Modelo guardado exitosamente como iJN1463.xml")
except Exception as e:
    print(f"❌ Error: {e}")
