"""Tests de las funciones compartidas extraidas en la refactorizacion
(evaluador FVA, genes candidatos y CLI sin GA duplicado)."""

import argparse
import random

import metabo_engineer as me
import cli_experto


def test_genes_candidatos_activos_ordenados_y_validos(modelo_textbook):
    genes = me.genes_candidatos_activos(modelo_textbook)
    assert genes, "El WT de textbook debe tener genes activos"
    assert genes == sorted(genes)
    ids_modelo = {g.id for g in modelo_textbook.genes}
    assert set(genes) <= ids_modelo


def test_evaluar_acoplamiento_fva_usa_cache(modelo_textbook, gen_biomasa_textbook):
    genes = me.genes_candidatos_activos(modelo_textbook)[:1]
    cache = {}
    m = modelo_textbook.copy()
    m.reactions.get_by_id("EX_succ_e").bounds = (0, 1000.0)

    a = me.evaluar_acoplamiento_fva(m, genes, "EX_succ_e", gen_biomasa_textbook, _cache=cache)
    assert isinstance(a, float) and a >= 0.0
    assert len(cache) == 1

    b = me.evaluar_acoplamiento_fva(m, list(reversed(genes)), "EX_succ_e", gen_biomasa_textbook, _cache=cache)
    assert a == b
    assert len(cache) == 1


def test_cli_metodos_usan_motor_compartido(modelo_textbook, monkeypatch, tmp_path):
    """El CLI debe devolver el mismo tipo de resultado que _evolucionar."""
    monkeypatch.chdir(tmp_path)  # el cache de esencialidad se escribe en el cwd
    random.seed(1)
    args = argparse.Namespace(generaciones=2, poblacion=6, mutaciones=2)

    for metodo in (cli_experto.metodo_ia, cli_experto.metodo_optknock):
        m = modelo_textbook.copy()
        mejor = metodo(m, "EX_succ_e", args)
        assert set(mejor) == {"genes", "fitness"}
        assert len(mejor["genes"]) == len(set(mejor["genes"]))


def test_cli_no_reimplementa_el_ga():
    import inspect
    src = inspect.getsource(cli_experto)
    assert "_evolucionar(" in src
    assert "random.sample(" not in src
