"""
Fixtures compartidos para todos los tests.

conftest.py es un archivo especial de pytest: los fixtures definidos aquí
están disponibles automáticamente en todos los archivos test_*.py.
"""

import pytest
from cobra.io import load_model

import metabo_engineer as me


@pytest.fixture(scope="session")
def modelo_textbook():
    """
    Carga el modelo textbook UNA vez para toda la sesión de tests.

    - 'textbook' es un alias de 'e_coli_core', incluido en COBRApy.
    - No requiere internet ni archivo externo → funciona en CI sin curl.
    """
    return load_model("textbook")


@pytest.fixture
def gen_biomasa_textbook(modelo_textbook):
    """ID de la reacción de biomasa del modelo textbook."""
    return [r.id for r in modelo_textbook.reactions if "biomass" in r.id.lower()][0]


@pytest.fixture(autouse=True)
def sin_termodinamica_real(monkeypatch):
    """Los tests no deben depender de eQuilibrator ni de la red.

    En CI con Python 3.10 pip instala equilibrator-api 0.6.0, cuya primera consulta descarga la
    base de compuestos (~1.3 GB): un test numerico de analisis_profundo llegaba a tardar casi
    8 minutos. La termodinamica se prueba aparte, en tests/test_termodinamica.py.
    """
    monkeypatch.setattr(me, "analizar_termodinamica_cepa", lambda *args, **kwargs: None)
