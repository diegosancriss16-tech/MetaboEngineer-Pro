"""
Tests de regresión basados en el feedback del revisor.
Comprueban: cálculo FVA-min (garantizado), bloqueo aerobio y reproducibilidad real.
"""
import pytest
import os
import sys
import re
import subprocess
from cobra.io import load_model
import metabo_engineer as me

@pytest.fixture(scope="module")
def modelo_textbook():
    return load_model("textbook")

def strip_ansi(text):
    """Elimina los códigos de color ANSI (Rich) del texto capturado."""
    return re.sub(r'\x1b\[[0-9;]*m', '', text)

def test_ambiente_aerobio_bloquea_fermentacion(modelo_textbook):
    """Verifica que el ambiente aerobio bloquee la fermentación en (0.0, 0.0)."""
    modelo = me.configurar_ambiente(modelo_textbook, "aerobio")
    for rxn_id in ["EX_etoh_e", "EX_lac__D_e", "EX_for_e"]:
        if rxn_id in modelo.reactions:
            bounds = modelo.reactions.get_by_id(rxn_id).bounds
            assert bounds == (0.0, 0.0), f"{rxn_id} no fue bloqueada correctamente."

def test_analisis_profundo_calcula_minimo_garantizado(modelo_textbook, capsys):
    """Verifica numéricamente que el valor impreso como 'Garantizado' sea igual a evaluar_bacteria y <= Cota Superior."""
    gen_biomasa = [r.id for r in modelo_textbook.reactions if "biomass" in r.id.lower()][0]
    quimico = "EX_succ_e"
    modelo_textbook.reactions.get_by_id(quimico).bounds = (0, 1000.0)
    
    # Ejecutamos la impresión de consola
    me.analisis_profundo(modelo_textbook, gen_biomasa, [], quimico)
    captured = capsys.readouterr()
    clean_out = strip_ansi(captured.out)
    
    # 1. Comprobar etiquetas
    assert "Garantizado @100%" in clean_out
    assert "Cota Superior" in clean_out
    
    # 2. Extraer valores numéricos mediante Regex
    match_garantizado = re.search(r"Garantizado @100%\):\s*([\d.]+)", clean_out)
    match_cota = re.search(r"Cota Superior \(@90% crecimiento\):\s*([\d.]+)", clean_out)
    
    assert match_garantizado is not None, "No se encontró el valor garantizado en la salida."
    assert match_cota is not None, "No se encontró el valor de la cota superior."
    
    val_garantizado = float(match_garantizado.group(1))
    val_cota = float(match_cota.group(1))
    
    # 3. Comparar con evaluar_bacteria (el estándar de oro del GA)
    fitness_ga = me.evaluar_bacteria(modelo_textbook, [], quimico, gen_biomasa)
    
    assert abs(val_garantizado - max(0.0, fitness_ga)) < 1e-4, \
        f"El valor impreso ({val_garantizado}) no coincide con evaluar_bacteria ({fitness_ga})."
        
    # 4. El mínimo garantizado debe ser <= a la cota superior (sentido físico)
    assert val_garantizado <= val_cota + 1e-4, \
        f"El mínimo garantizado ({val_garantizado}) es mayor que la cota superior ({val_cota})."

def test_reproducibilidad_semilla_evita_hash_random(tmp_path):
    """
    Ejecuta el GA real en subprocesos con diferentes PYTHONHASHSEED.
    Si el operador de cruce usa set(), los resultados variarán.
    """
    script_content = """
import sys
import io
import metabo_engineer as me
from cobra.io import load_model
from rich.prompt import Prompt, IntPrompt

# Mock inputs para que no pida interactividad
Prompt.ask = lambda *args, **kwargs: "EX_succ_e"
IntPrompt.ask = lambda *args, **kwargs: 1 # Presupuesto pequeño

# Mock genes para evitar cálculo pesado
me.obtener_genes_no_esenciales = lambda modelo, gen_biomasa: ["genA", "genB", "genC", "genD", "genE", "genF", "genG", "genH"]

# Mock fitness que depende fuertemente del orden de los genes (usando ord)
def mock_eval(modelo, genes, quimico, biomasa):
    return float(sum(ord(c) for g in genes for c in g) % 1000)

me.evaluar_bacteria = mock_eval
me.analisis_profundo = lambda *args, **kwargs: None
me._exportar_resultado_json = lambda *args, **kwargs: None

modelo = load_model("textbook")

# Capturamos la salida de consola del GA
old_stdout = sys.stdout
sys.stdout = io.StringIO()
me.disenar_cepa_ia(modelo, "textbook", "textbook", seed=42)
output = sys.stdout.getvalue()
sys.stdout = old_stdout

# Imprimimos solo la línea del fitness final para comparar
for line in output.split("\\n"):
    if "Diseno completado" in line:
        print(line.strip())
        break
"""
    script_path = tmp_path / "repro_test.py"
    script_path.write_text(script_content)
    
    # Ejecutamos con PYTHONHASHSEED=1
    env1 = {**os.environ, "PYTHONHASHSEED": "1", "PYTHONPATH": os.getcwd()}
    res1 = subprocess.run([sys.executable, str(script_path)], capture_output=True, text=True, env=env1)
    
    # Ejecutamos con PYTHONHASHSEED=2
    env2 = {**os.environ, "PYTHONHASHSEED": "2", "PYTHONPATH": os.getcwd()}
    res2 = subprocess.run([sys.executable, str(script_path)], capture_output=True, text=True, env=env2)
    
    out1 = [line for line in res1.stdout.split("\n") if "Diseno completado" in line]
    out2 = [line for line in res2.stdout.split("\n") if "Diseno completado" in line]
    
    assert len(out1) > 0 and len(out2) > 0, f"El GA no produjo salida esperada. Error: {res1.stderr}"
    assert out1[0] == out2[0], f"La reproducibilidad falló entre distintos PYTHONHASHSEED: {out1[0]} vs {out2[0]}"