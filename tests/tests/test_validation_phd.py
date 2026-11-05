"""
Pruebas de Validación Nivel PhD - Fase 1 (A2, C5, C7)
Basado en las recomendaciones del peer review.
"""
import pytest
import cobra
import metabo_engineer as me

BIO = "Biomass_Ecoli_core"

@pytest.fixture(scope="module")
def textbook():
    modelo, _ = me.cargar_modelo_universal("textbook")
    return modelo

# =====================================================================
# A2. LÍNEA BASE IDÉNTICA A COBRAPY PURO
# =====================================================================
class TestLineaBase:
    def test_crecimiento_WT_coincide_con_cobrapy(self, textbook):
        """El crecimiento del WT debe ser exactamente igual a COBRApy puro."""
        # Nuestro orquestador
        sol_orquestador = textbook.slim_optimize()
        
        # COBRApy puro (copiando el modelo de nuevo)
        m_puro = cobra.io.load_model("textbook")
        sol_puro = m_puro.slim_optimize()
        
        assert sol_orquestador == pytest.approx(sol_puro, abs=1e-6)

    def test_crecimiento_anaerobio_coincide_con_cobrapy(self, textbook):
        """El crecimiento anaerobio debe coincidir con COBRApy puro."""
        m_orq = me.configurar_ambiente(textbook, "anaerobio")
        sol_orq = m_orq.slim_optimize()
        
        m_puro = cobra.io.load_model("textbook")
        m_puro.reactions.EX_o2_e.bounds = (0, 1000)
        sol_puro = m_puro.slim_optimize()
        
        assert sol_orq == pytest.approx(sol_puro, abs=1e-6)


# =====================================================================
# C5. PRUEBAS METAMÓRFICAS
# =====================================================================
class TestMetamorficas:
    def test_escalado_de_glucosa(self, textbook):
        """Si duplicamos la glucosa, los flujos se duplican, pero el rendimiento no cambia."""
        m = textbook.copy()
        m.reactions.EX_glc__D_e.lower_bound = -10
        sol_10 = m.optimize()
        
        m.reactions.EX_glc__D_e.lower_bound = -20
        sol_20 = m.optimize()
        
        # El crecimiento debe duplicarse
        assert sol_20.objective_value == pytest.approx(sol_10.objective_value * 2, rel=1e-4)

    def test_reordenar_reacciones_no_cambia_resultado(self, textbook):
        """El orden de las reacciones en el diccionario no debe afectar el FBA."""
        m = textbook.copy()
        # Forzamos un reordenamiento interno (esto simula cambiar el orden en el XML)
        reacciones_originales = list(m.reactions)
        m.reactions._reaction_dict.clear()
        for r in reversed(reacciones_originales):
            m.reactions._reaction_dict[r.id] = r
            
        assert m.slim_optimize() == pytest.approx(textbook.slim_optimize(), abs=1e-6)


# =====================================================================
# C7. CONTROLES NEGATIVOS
# =====================================================================
class TestControlesNegativos:
    def test_metabolito_sin_ruta_da_cero(self, textbook):
        """Si pedimos un químico que no existe, el fitness debe ser 0, no un error."""
        fit = me.evaluar_bacteria(textbook, [], "EX_metabolito_falso_e", BIO)
        assert fit == 0.0

    def test_knockout_letal_da_fitness_cero(self, textbook):
        """Apagar un gen esencial debe dar fitness 0."""
        # eno (b2779) es esencial para la glucólisis
        fit = me.evaluar_bacteria(textbook, ["b2779"], "EX_succ_e", BIO)
        assert fit == 0.0

    def test_no_supera_el_maximo_estequiometrico(self, textbook):
        """El fitness nunca debe superar el máximo teórico del FBA (no hay energía gratis)."""
        m = textbook.copy()
        m.reactions.EX_succ_e.bounds = (0, 1000)
        max_succ = m.slim_optimize(error_value=0.0)
        
        # Evaluamos la cepa WT (sin knockouts)
        fit_wt = me.evaluar_bacteria(m, [], "EX_succ_e", BIO)
        
        # El fitness no puede ser mayor al máximo absoluto
        assert fit_wt <= max_succ + 1e-6