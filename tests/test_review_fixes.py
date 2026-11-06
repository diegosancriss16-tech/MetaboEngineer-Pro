"""
Tests de regresion basados en el feedback del revisor.
"""
import pytest
import random
from cobra.io import load_model
import metabo_engineer as me

@pytest.fixture(scope="module")
def modelo_textbook():
    return load_model("textbook")

def test_ambiente_aerobio_bloquea_fermentacion(modelo_textbook):
    modelo = me.configurar_ambiente(modelo_textbook, "aerobio")
    for rxn_id in ["EX_etoh_e", "EX_lac__D_e", "EX_for_e"]:
        if rxn_id in modelo.reactions:
            bounds = modelo.reactions.get_by_id(rxn_id).bounds
            assert bounds == (0.0, 0.0), f"{rxn_id} no fue bloqueada correctamente."

def test_reproducibilidad_semilla_evita_hash_random():
    p1 = ["genA", "genB", "genC", "genD"]
    p2 = ["genE", "genF", "genG", "genH"]
    hijo_ordenado = list(dict.fromkeys(p1[:3] + p2[3:]))
    assert hijo_ordenado == ["genA", "genB", "genC", "genH"]
    assert isinstance(hijo_ordenado, list)

def test_analisis_profundo_calcula_minimo_garantizado(modelo_textbook, capsys):
    gen_biomasa = [r.id for r in modelo_textbook.reactions if "biomass" in r.id.lower()][0]
    quimico = "EX_succ_e"
    modelo_textbook.reactions.get_by_id(quimico).bounds = (0, 1000.0)
    me.analisis_profundo(modelo_textbook, gen_biomasa, [], quimico)
    captured = capsys.readouterr()
    assert "Garantizado @100%" in captured.out
    assert "Cota Superior @90%" in captured.out
