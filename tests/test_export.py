"""
Tests de exportación a JSON/CSV.

Cubre:
- Exportar corrida genera JSON individual
- Exportar corrida hace append a historial.csv
- Múltiples exportaciones acumulan filas en CSV
"""

import json
import os
import pandas as pd
import pytest

import metabo_engineer as me


class TestExportarResultado:
    def _setup_entorno(self, tmp_path, monkeypatch):
        """Helper: redirige historial y CWD a tmp_path."""
        historial = tmp_path / "historial.csv"
        monkeypatch.setattr(me, "_RUTA_HISTORIAL_CSV", str(historial))
        monkeypatch.chdir(tmp_path)
        os.makedirs(tmp_path / "reportes", exist_ok=True)
        return historial

    def test_exportar_genera_json(self, tmp_path, monkeypatch):
        """Debe crear un archivo .json con los datos correctos."""
        self._setup_entorno(tmp_path, monkeypatch)

        me._exportar_resultado_json(
            modelo_id="textbook",
            quimico="EX_succ_e",
            metodo="ia",
            genes=["b0001", "b0002"],
            fitness=8.5,
            seed=42,
            parametros_extra={"generaciones": 10, "poblacion": 20},
        )

        jsons = list((tmp_path / "reportes").glob("*.json"))
        assert len(jsons) == 1

        with open(jsons[0]) as f:
            data = json.load(f)

        assert data["modelo"] == "textbook"
        assert data["quimico"] == "EX_succ_e"
        assert data["metodo"] == "ia"
        assert data["fitness"] == 8.5
        assert data["seed"] == 42
        assert data["genes_knockout"] == ["b0001", "b0002"]
        assert data["num_genes"] == 2

    def test_exportar_genera_csv_primera_vez(self, tmp_path, monkeypatch):
        """Primera exportación → CSV con header."""
        historial = self._setup_entorno(tmp_path, monkeypatch)

        me._exportar_resultado_json("textbook", "EX_succ_e", "ia", ["b1"], 8.0, 42)

        assert historial.exists()
        df = pd.read_csv(historial)
        assert len(df) == 1
        assert "fecha" in df.columns
        assert "fitness" in df.columns
        assert df.iloc[0]["fitness"] == 8.0

    def test_exportar_append_csv_segunda_vez(self, tmp_path, monkeypatch):
        """Segunda exportación → 2 filas en CSV (sin repetir header)."""
        historial = self._setup_entorno(tmp_path, monkeypatch)

        me._exportar_resultado_json("textbook", "EX_succ_e", "ia", ["b1"], 8.0, 42)
        me._exportar_resultado_json("textbook", "EX_succ_e", "ia", ["b2"], 9.0, 42)

        df = pd.read_csv(historial)
        assert len(df) == 2
        assert list(df["fitness"]) == [8.0, 9.0]

    def test_exportar_genes_como_string_en_csv(self, tmp_path, monkeypatch):
        """Los genes se guardan como string separado por ';' en CSV."""
        historial = self._setup_entorno(tmp_path, monkeypatch)

        me._exportar_resultado_json("textbook", "EX_succ_e", "ia", ["b1", "b2", "b3"], 8.0, 42)

        df = pd.read_csv(historial)
        assert df.iloc[0]["genes"] == "b1;b2;b3"