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
    return re.sub(r'\x1b\[[0-9;]*m', '', text)

def test_aerobio_reabre_oxigeno():
    """Verifica que pasar de anaerobio a aerobio restaure los límites del oxígeno."""
    modelo = load_model("textbook")
    modelo_anaerobio = me.configurar_ambiente(modelo, "anaerobio")
    modelo_aerobio = me.configurar_ambiente(modelo_anaerobio, "aerobio")
    
    assert modelo_aerobio.reactions.EX_o2_e.bounds == (-1000.0, 1000.0)
    assert modelo_aerobio.slim_optimize() == pytest.approx(0.8739, abs=1e-3)

def test_cache_cambia_con_el_medio(tmp_path, monkeypatch):
    """Verifica que el hash del medio invalide el caché si cambia la glucosa."""
    monkeypatch.chdir(tmp_path)
    m = load_model("textbook")
    m2 = m.copy()
    m2.reactions.EX_glc__D_e.lower_bound = -5
    
    me.obtener_genes_no_esenciales(m, "Biomass_Ecoli_core")
    me.obtener_genes_no_esenciales(m2, "Biomass_Ecoli_core")
    
    cache_files = list((tmp_path / "benchmark/results").glob("essential_cache_*.json"))
    assert len(cache_files) == 2

def test_analisis_profundo_calcula_minimo_garantizado(modelo_textbook, capsys):
    """Verifica el cálculo numérico y el regex de las etiquetas actuales."""
    gen_biomasa = [r.id for r in modelo_textbook.reactions if "biomass" in r.id.lower()][0]
    quimico = "EX_succ_e"
    modelo_textbook.reactions.get_by_id(quimico).bounds = (0, 1000.0)
    
    # Corremos analisis_profundo con lista vacía (comportamiento WT)
    me.analisis_profundo(modelo_textbook, gen_biomasa, [], quimico)
    captured = capsys.readouterr()
    clean_out = strip_ansi(captured.out)
    
    assert "Garantizado @100%" in clean_out
    assert "Cota Superior" in clean_out
    
    match_garantizado = re.search(r"Garantizado @100%\):\s*([\d.]+)", clean_out)
    match_cota = re.search(r"Cota Superior \(@90% max\):\s*([\d.]+)", clean_out)
    
    assert match_garantizado is not None, "No se encontró el valor garantizado."
    assert match_cota is not None, "No se encontró la cota superior."
    
    val_garantizado = float(match_garantizado.group(1))
    val_cota = float(match_cota.group(1))
    
    # Calcular fitness del WT directamente (sin evaluar_bacteria, que ahora exige >=1 gen)
    modelo_wt = modelo_textbook.copy()
    modelo_wt.objective = gen_biomasa
    sol_bio = modelo_wt.optimize()
    max_bio = sol_bio.objective_value
    modelo_wt.reactions.get_by_id(gen_biomasa).lower_bound = max_bio * 1.0
    modelo_wt.objective = quimico
    modelo_wt.objective.direction = 'min'
    sol_prod = modelo_wt.optimize()
    fitness_wt = max(0.0, sol_prod.objective_value)
    
    assert abs(val_garantizado - fitness_wt) < 1e-4
    assert val_garantizado <= val_cota + 1e-4

def test_reproducibilidad_semilla_evita_hash_random(tmp_path):
    """Ejecuta el GA en subprocesos con distinto PYTHONHASHSEED."""
    script_content = """
import sys
import io
import random
import metabo_engineer as me
from cobra.io import load_model
from rich.prompt import Prompt, IntPrompt

modelo = load_model("textbook")
todos_los_genes = [g.id for g in modelo.genes][:20]
me.obtener_genes_no_esenciales = lambda modelo, gen_biomasa: todos_los_genes

def mock_eval(*args, **kwargs):
    genes = args[1]
    return float(sum(ord(c) for g in genes for c in g) % 1000)

me.evaluar_bacteria = mock_eval

respuestas_int = iter([5, 12, 4])
def fake_int_ask(*args, **kwargs):
    return next(respuestas_int)
IntPrompt.ask = fake_int_ask
Prompt.ask = lambda *args, **kwargs: "EX_succ_e"

me.analisis_profundo = lambda *args, **kwargs: None
me._exportar_resultado_json = lambda *args, **kwargs: None

random.seed(42)
mejor = me.disenar_cepa_ia(modelo, "textbook", "textbook", seed=42)
print("MEJOR:" + ",".join(mejor["genes"]))
"""
    script_path = tmp_path / "repro_test.py"
    # FIX: Forzar UTF-8 para evitar el SyntaxError por caracteres no ASCII
    script_path.write_text(script_content, encoding='utf-8')
    
    env1 = {**os.environ, "PYTHONHASHSEED": "1", "PYTHONPATH": os.getcwd()}
    res1 = subprocess.run([sys.executable, str(script_path)], capture_output=True, text=True, env=env1, timeout=180)
    
    env2 = {**os.environ, "PYTHONHASHSEED": "2", "PYTHONPATH": os.getcwd()}
    res2 = subprocess.run([sys.executable, str(script_path)], capture_output=True, text=True, env=env2, timeout=180)
    
    out1 = [line for line in res1.stdout.split("\n") if line.startswith("MEJOR:")]
    out2 = [line for line in res2.stdout.split("\n") if line.startswith("MEJOR:")]
    
    assert len(out1) > 0, f"Error en ejecución 1: {res1.stderr}"
    assert len(out2) > 0, f"Error en ejecución 2: {res2.stderr}"
    assert out1[0] == out2[0], f"Reproducibilidad falló: {out1[0]} vs {out2[0]}"