"""
Tests de regresión basados en el feedback del revisor.
Comprueban: cálculo FVA-min, bloqueo aerobio, reproducibilidad real y relleno seguro.
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
    """Verifica numéricamente que el 'Garantizado' sea igual a evaluar_bacteria y <= Cota Superior."""
    gen_biomasa = [r.id for r in modelo_textbook.reactions if "biomass" in r.id.lower()][0]
    quimico = "EX_succ_e"
    modelo_textbook.reactions.get_by_id(quimico).bounds = (0, 1000.0)
    
    me.analisis_profundo(modelo_textbook, gen_biomasa, [], quimico)
    captured = capsys.readouterr()
    clean_out = strip_ansi(captured.out)
    
    assert "Garantizado @100%" in clean_out
    assert "Cota Superior" in clean_out
    
    match_garantizado = re.search(r"Garantizado @100%\):\s*([\d.]+)", clean_out)
    match_cota = re.search(r"Cota Superior \(@90% crecimiento\):\s*([\d.]+)", clean_out)
    
    assert match_garantizado is not None, "No se encontró el valor garantizado."
    assert match_cota is not None, "No se encontró la cota superior."
    
    val_garantizado = float(match_garantizado.group(1))
    val_cota = float(match_cota.group(1))
    fitness_ga = me.evaluar_bacteria(modelo_textbook, [], quimico, gen_biomasa)
    
    assert abs(val_garantizado - max(0.0, fitness_ga)) < 1e-4, "El impreso no coincide con evaluar_bacteria."
    assert val_garantizado <= val_cota + 1e-4, "El mínimo garantizado > cota superior."

def test_reproducibilidad_semilla_evita_hash_random(tmp_path):
    """
    Ejecuta el GA real con genes reales y diferentes PYTHONHASHSEED.
    Compara la lista completa de genes del mejor individuo.
    """
    script_content = """
import sys
import io
import random
import metabo_engineer as me
from cobra.io import load_model
from rich.prompt import Prompt, IntPrompt

# 1. Cargar modelo y usar genes REALES (primeros 20)
modelo = load_model("textbook")
todos_los_genes = [g.id for g in modelo.genes][:20]
me.obtener_genes_no_esenciales = lambda modelo, gen_biomasa: todos_los_genes

# 2. Mock de fitness que depende del orden de los genes
def mock_eval(modelo, genes, quimico, biomasa):
    return float(sum(ord(c) for g in genes for c in g) % 1000)
me.evaluar_bacteria = mock_eval

# 3. Forzar parámetros del GA: 5 generaciones, 12 población, 4 genes
respuestas_int = iter([5, 12, 4])
def fake_int_ask(*args, **kwargs):
    return next(respuestas_int)
IntPrompt.ask = fake_int_ask
Prompt.ask = lambda *args, **kwargs: "EX_succ_e"

# 4. Evitar efectos secundarios
me.analisis_profundo = lambda *args, **kwargs: None
me._exportar_resultado_json = lambda *args, **kwargs: None

# 5. Ejecutar GA
random.seed(42)
mejor = me.disenar_cepa_ia(modelo, "textbook", "textbook", seed=42)
print("MEJOR:" + ",".join(mejor["genes"]))
"""
    script_path = tmp_path / "repro_test.py"
    script_path.write_text(script_content)
    
    env1 = {**os.environ, "PYTHONHASHSEED": "1", "PYTHONPATH": os.getcwd()}
    res1 = subprocess.run([sys.executable, str(script_path)], capture_output=True, text=True, env=env1, timeout=30)
    
    env2 = {**os.environ, "PYTHONHASHSEED": "2", "PYTHONPATH": os.getcwd()}
    res2 = subprocess.run([sys.executable, str(script_path)], capture_output=True, text=True, env=env2, timeout=30)
    
    out1 = [line for line in res1.stdout.split("\n") if line.startswith("MEJOR:")]
    out2 = [line for line in res2.stdout.split("\n") if line.startswith("MEJOR:")]
    
    assert len(out1) > 0, f"Error en ejecución 1: {res1.stderr}"
    assert len(out2) > 0, f"Error en ejecución 2: {res2.stderr}"
    assert out1[0] == out2[0], f"Reproducibilidad falló: {out1[0]} vs {out2[0]}"

def test_relleno_no_se_cuelga_con_pocos_genes(tmp_path):
    """
    Fuerza el caso donde hay menos genes disponibles que k.
    Si el código no usa min(), el subprocess se colgará y el timeout lo matará.
    """
    script_content = """
import random
import metabo_engineer as me
from cobra.io import load_model
from rich.prompt import Prompt, IntPrompt

modelo = load_model("textbook")
# Solo 2 genes disponibles
me.obtener_genes_no_esenciales = lambda modelo, gen_biomasa: ["b0001", "b0002"]
me.evaluar_bacteria = lambda *args: 1.0  # Fitness trivial para rapidez
me.analisis_profundo = lambda *args: None
me._exportar_resultado_json = lambda *args: None

respuestas_int = iter([2, 5, 4])  # 2 gens, 5 pop, 4 mutaciones (k > genes disponibles)
IntPrompt.ask = lambda *args, **kwargs: next(respuestas_int)
Prompt.ask = lambda *args, **kwargs: "EX_succ_e"

random.seed(42)
me.disenar_cepa_ia(modelo, "textbook", "textbook", seed=42)
print("OK")
"""
    script_path = tmp_path / "stress_test.py"
    script_path.write_text(script_content)
    
    env = {**os.environ, "PYTHONPATH": os.getcwd()}
    try:
        res = subprocess.run([sys.executable, str(script_path)], capture_output=True, text=True, env=env, timeout=10)
        assert "OK" in res.stdout
    except subprocess.TimeoutExpired:
        pytest.fail("El GA se colgó. El bucle while no tiene límite cuando k > genes disponibles.")