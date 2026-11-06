"""
Fixtures compartidos para todos los tests.

conftest.py es un archivo especial de pytest: los fixtures definidos aquí
están disponibles automáticamente en todos los archivos test_*.py.
"""

import pytest
from cobra.io import load_model


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