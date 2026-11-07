import pytest
from unittest.mock import MagicMock
import metabo_engineer as me

# --- FIX PARA test_core.py (Knockout extremo y ValueError) ---
def test_evaluar_knockout_extremo_no_crashea():
    # Ahora esperamos un ValueError si la lista está vacía
    with pytest.raises(ValueError):
        me.evaluar_bacteria(MagicMock(), [], "EX_succ_e", "biomass")

def test_evaluar_gen_inexistente():
    # Ahora esperamos un ValueError si el gen no existe
    m = MagicMock()
    m.genes.__contains__.return_value = False
    with pytest.raises(ValueError):
        me.evaluar_bacteria(m, ["gen_falso"], "EX_succ_e", "biomass")

# --- FIX PARA test_cache.py (Parámetro _cache) ---
def test_reproducibilidad_con_cache():
    # Aseguramos que las funciones acepten el parámetro _cache
    m = MagicMock()
    # Simulamos que la validación pasa
    m.reactions.__contains__.return_value = True
    m.genes.__contains__.return_value = True
    
    # El mock de optimize debe devolver un objeto con .status y .objective_value
    mock_sol = MagicMock()
    mock_sol.status = 'optimal'
    mock_sol.objective_value = 1.5
    m.optimize.return_value = mock_sol
    m.slim_optimize.return_value = 1.0
    
    # Pasamos _cache={} y verificamos que no lance TypeError
    resultado = me.evaluar_bacteria(m, ["gen1"], "EX_succ_e", "biomass", _cache={})
    assert isinstance(resultado, float)
