"""
Fixtures compartidos para todos los tests.

conftest.py es un archivo especial de pytest: los fixtures definidos aquí
están disponibles automáticamente en todos los archivos test_*.py.
"""

import pytest
import metabo_engineer as me


@pytest.fixture(scope="session")
def modelo_textbook():
    """
    Carga el modelo textbook UNA vez para toda la sesión de tests.
    Es local (no requiere internet) y pequeño (95 reacciones) → rápido.
    """
    modelo, _ = me.cargar_modelo_universal("textbook")
    return modelo


@pytest.fixture
def gen_biomasa_textbook(modelo_textbook):
    """ID de la reacción de biomasa del modelo textbook."""
    return [r.id for r in modelo_textbook.reactions if "biomass" in r.id.lower()][0]