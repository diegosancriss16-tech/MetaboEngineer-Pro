# tests/test_regression.py
import pytest
import os
import sys
import hashlib
from cobra.io import load_model
import metabo_engineer as me

def test_ambiente_aerobio_bloquea_fermentacion():
    modelo = load_model("textbook")
    modelo = me.configurar_ambiente(modelo, "aerobio")
    for rxn in ["EX_etoh_e", "EX_lac__D_e", "EX_for_e"]:
        if rxn in modelo.reactions:
            assert modelo.reactions.get_by_id(rxn).bounds == (0.0, 0.0)

def test_cargar_modelo_textbook_ok():
    modelo, msg = me.cargar_modelo_universal("textbook")
    assert modelo is not None
    assert len(modelo.reactions) > 0

def test_cargar_modelo_vacio_falla():
    with pytest.raises(Exception) as exc:
        me.cargar_modelo_universal("")
    assert "vacio" in str(exc.value).lower() or "invalido" in str(exc.value).lower()

def test_get_organism_id():
    assert me.get_organism_id("iYO844") == "224308"
    assert me.get_organism_id("iMM904") == "559292"
    assert me.get_organism_id("desconocido") == "83333"

def test_sha256_archivo(tmp_path):
    a1 = tmp_path / "a.txt"
    a1.write_text("hola mundo")
    a2 = tmp_path / "b.txt"
    a2.write_text("hola mundo")
    h1 = me._sha256_archivo(str(a1))
    h2 = me._sha256_archivo(str(a2))
    assert h1 == h2
    assert len(h1) == 64

def test_relleno_no_se_cuelga_con_pocos_genes(tmp_path, monkeypatch):
    # Este test protege el bucle while de relleno
    import subprocess
    script_content = """
import random
import metabo_engineer as me
from cobra.io import load_model

modelo = load_model("textbook")
me.obtener_genes_no_esenciales = lambda m, b: ["b0001", "b0002"] # Solo 2 genes
me.evaluar_bacteria = lambda *args, **kwargs: 1.0
me.analisis_profundo = lambda *args, **kwargs: None
me._exportar_resultado_json = lambda *args, **kwargs: None

random.seed(42)
# k=4, pero solo hay 2 genes. No debe colgarse.
mejor = me._evolucionar(
    pool_genes=["b0001", "b0002"], k=4, generaciones=1, pop_size=5, 
    fitness_fn=me.evaluar_bacteria
)
print("OK")
"""
    script_path = tmp_path / "stress_test.py"
    script_path.write_text(script_content, encoding='utf-8')
    env = {**os.environ, "PYTHONPATH": os.getcwd()}
    res = subprocess.run([sys.executable, str(script_path)], capture_output=True, text=True, env=env, timeout=10)
    assert "OK" in res.stdout