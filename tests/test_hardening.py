"""Tests de endurecimiento: solver sin fallo silencioso, avisos del CLI,
modulo compartido de organismos y compatibilidad del alias de OptKnock."""

from unittest.mock import MagicMock

import pytest

import cli_experto
import metabo_engineer as me
from organismos import get_organism_id


class TestSolverFuncional:
    def test_modelo_sano_se_devuelve_sin_cambios(self, modelo_textbook):
        m = modelo_textbook.copy()
        assert me._asegurar_solver_funcional(m) is m

    def test_lanza_si_ningun_solver_optimiza(self):
        m = MagicMock()
        m.slim_optimize.return_value = float("nan")
        with pytest.raises(RuntimeError, match="Ningun solver"):
            me._asegurar_solver_funcional(m)

    def test_cae_a_glpk_si_el_primero_falla(self):
        m = MagicMock()
        m.slim_optimize.side_effect = [float("nan"), 0.87]
        assert me._asegurar_solver_funcional(m) is m
        assert m.solver == "glpk"

    def test_ambiente_imposible_lanza_en_vez_de_continuar(self, modelo_textbook, monkeypatch):
        monkeypatch.setattr(me, "_asegurar_solver_funcional",
                            lambda m: (_ for _ in ()).throw(RuntimeError("sin solver")))
        with pytest.raises(RuntimeError):
            me.configurar_ambiente(modelo_textbook, "aerobio")


class TestVerificarResultadoCLI:
    def test_fallo_total_devuelve_codigo_5(self, capsys):
        assert cli_experto.verificar_resultado({"fitness": -1e6, "genes": ["g"]}, "ia") == 5
        assert "TODAS las evaluaciones fallaron" in capsys.readouterr().err

    def test_sin_acoplamiento_avisa_pero_no_es_error(self, capsys):
        assert cli_experto.verificar_resultado({"fitness": 0.0, "genes": ["g"]}, "ia") == 0
        salida = capsys.readouterr().out
        assert "WARN" in salida and "[OK]" not in salida

    def test_con_acoplamiento_es_ok(self, capsys):
        assert cli_experto.verificar_resultado({"fitness": 1.5, "genes": ["g"]}, "optknock") == 0
        assert "[OK] OPTKNOCK completado" in capsys.readouterr().out


class TestOrganismos:
    @pytest.mark.parametrize("modelo_id, taxon", [
        ("iYO844", "224308"), ("iCW773", "196627"), ("iMM904", "559292"),
        ("iEK1008", "1773"), ("iML1515", "83333"), ("desconocido", "83333"),
        ("", "83333"), (None, "83333"),
    ])
    def test_taxon(self, modelo_id, taxon):
        assert get_organism_id(modelo_id) == taxon

    def test_misma_funcion_en_todos_los_modulos(self):
        pytest.importorskip("equilibrator_api")  # reporte_investigador la importa al cargarse
        import reporte_investigador
        assert me.get_organism_id is get_organism_id
        assert reporte_investigador.get_organism_id is get_organism_id

    def test_metabo_engineer_reexporta_la_funcion_compartida(self):
        assert me.get_organism_id is get_organism_id


def test_alias_optknock_exacto_retrocompatible():
    assert me.optknock_exacto_milp is me.optknock_exacto
