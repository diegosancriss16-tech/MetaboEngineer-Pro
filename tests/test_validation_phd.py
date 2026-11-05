import pytest
import cobra
import metabo_engineer as me

BIO = "Biomass_Ecoli_core"

@pytest.fixture(scope="module")
def textbook():
    modelo, _ = me.cargar_modelo_universal("textbook")
    return modelo

class TestLineaBase:
    def test_crecimiento_WT_coincide_con_cobrapy(self, textbook):
        sol_orquestador = textbook.slim_optimize()
        m_puro = cobra.io.load_model("textbook")
        sol_puro = m_puro.slim_optimize()
        assert sol_orquestador == pytest.approx(sol_puro, abs=1e-6)

    def test_crecimiento_anaerobio_coincide_con_cobrapy(self, textbook):
        m_orq = me.configurar_ambiente(textbook, "anaerobio")
        sol_orq = m_orq.slim_optimize()
        m_puro = cobra.io.load_model("textbook")
        m_puro.reactions.EX_o2_e.bounds = (0, 1000)
        sol_puro = m_puro.slim_optimize()
        assert sol_orq == pytest.approx(sol_puro, abs=1e-6)

class TestMetamorficas:
    def test_escalado_de_glucosa(self, textbook):
        """Si duplicamos la glucosa disponible, el flujo de consumo se duplica."""
        m = textbook.copy()
        m.reactions.EX_glc__D_e.lower_bound = -10
        sol_10 = m.optimize()
        flujo_10 = sol_10.fluxes['EX_glc__D_e']
        
        m.reactions.EX_glc__D_e.lower_bound = -20
        sol_20 = m.optimize()
        flujo_20 = sol_20.fluxes['EX_glc__D_e']
        
        assert flujo_20 == pytest.approx(flujo_10 * 2, rel=1e-3)

    def test_renombrar_reaccion_no_cambia_resultado(self, textbook):
        """Cambiar el ID de una reacción no debe afectar el resultado del FBA."""
        m = textbook.copy()
        val_original = m.slim_optimize()
        
        # Renombramos la primera reacción
        old_id = m.reactions[0].id
        m.reactions[0].id = 'NUEVO_ID_TEST'
        
        assert m.slim_optimize() == pytest.approx(val_original, abs=1e-6)

class TestControlesNegativos:
    def test_metabolito_sin_ruta_da_cero(self, textbook):
        fit = me.evaluar_bacteria(textbook, [], "EX_metabolito_falso_e", BIO)
        assert fit == 0.0

    def test_knockout_letal_da_fitness_cero(self, textbook):
        fit = me.evaluar_bacteria(textbook, ["b2779"], "EX_succ_e", BIO)
        assert fit == 0.0

    def test_no_supera_el_maximo_estequiometrico(self, textbook):
        """El fitness nunca debe superar el máximo teórico del FBA para ese químico."""
        m = textbook.copy()
        m.reactions.EX_succ_e.bounds = (0, 1000)
        
        # ¡CORRECCIÓN! Cambiamos el objetivo a succinato para calcular su máximo real
        m.objective = 'EX_succ_e'
        max_succ = m.slim_optimize(error_value=0.0)
        
        # Evaluamos la cepa WT
        fit_wt = me.evaluar_bacteria(m, [], "EX_succ_e", BIO)
        
        # El fitness no puede ser mayor al máximo estequiométrico
        assert fit_wt <= max_succ + 1e-6
