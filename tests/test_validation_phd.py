import pytest
import cobra
import random
import metabo_engineer as me

BIO_TEXTBOOK = "Biomass_Ecoli_core"

@pytest.fixture(scope="module")
def textbook():
    modelo, _ = me.cargar_modelo_universal("textbook")
    return modelo

class TestLineaBase:
    def test_crecimiento_WT_coincide_con_cobrapy(self, textbook):
        m = textbook.copy()
        sol_orquestador = m.slim_optimize()
        m_puro = cobra.io.load_model("textbook")
        sol_puro = m_puro.slim_optimize()
        assert sol_orquestador == pytest.approx(sol_puro, abs=1e-6)

    def test_crecimiento_anaerobio_coincide_con_cobrapy(self, textbook):
        m_orq = me.configurar_ambiente(textbook.copy(), "anaerobio")
        sol_orq = m_orq.slim_optimize()
        m_puro = cobra.io.load_model("textbook")
        m_puro.reactions.EX_o2_e.bounds = (0, 1000)
        sol_puro = m_puro.slim_optimize()
        assert sol_orq == pytest.approx(sol_puro, abs=1e-6)

class TestMetamorficas:
    def test_escalado_de_todas_las_cotas(self, textbook):
        m = textbook.copy()
        cotas = []
        for rxn in m.reactions:
            if rxn.lower_bound != 0 and rxn.lower_bound != -1000:
                cotas.append((rxn, 'lower', rxn.lower_bound))
            if rxn.upper_bound != 0 and rxn.upper_bound != 1000:
                cotas.append((rxn, 'upper', rxn.upper_bound))
        sol_1 = m.optimize()
        flujo_1 = sol_1.fluxes['EX_glc__D_e']
        for rxn, tipo, valor in cotas:
            if tipo == 'lower': rxn.lower_bound = valor * 2
            else: rxn.upper_bound = valor * 2
        sol_2 = m.optimize()
        flujo_2 = sol_2.fluxes['EX_glc__D_e']
        assert flujo_2 == pytest.approx(flujo_1 * 2, rel=1e-4)

    def test_renombrar_reaccion_no_cambia_resultado(self, textbook):
        m = textbook.copy()
        val_original = m.slim_optimize()
        m_nuevo = cobra.Model("test")
        for met in m.metabolites: m_nuevo.add_metabolites([met])
        reacciones = list(m.reactions)
        random.shuffle(reacciones)
        for rxn in reacciones: m_nuevo.add_reactions([rxn])
        m_nuevo.objective = BIO_TEXTBOOK
        assert m_nuevo.slim_optimize() == pytest.approx(val_original, abs=1e-6)

class TestControlesNegativos:
    def test_metabolito_sin_ruta_da_cero(self, textbook):
        fit = me.evaluar_bacteria(textbook.copy(), [], "EX_metabolito_falso_e", BIO_TEXTBOOK)
        assert fit == 0.0

    def test_knockout_letal_da_fitness_cero(self, textbook):
        fit = me.evaluar_bacteria(textbook.copy(), ["b2779"], "EX_succ_e", BIO_TEXTBOOK)
        assert fit == 0.0

    def test_no_supera_el_maximo_estequiometrico(self, textbook):
        m = textbook.copy()
        m.reactions.EX_succ_e.bounds = (0, 1000)
        m.objective = 'EX_succ_e'
        max_succ = m.slim_optimize(error_value=0.0)
        fit_wt = me.evaluar_bacteria(m, [], "EX_succ_e", BIO_TEXTBOOK)
        assert fit_wt <= max_succ + 1e-6

class TestE4Paper2022:
    def test_cepa_del_paper_en_sus_condiciones(self):
        m, _ = me.cargar_modelo_universal("iML1515")
        BIO_I1515 = "BIOMASS_Ec_iML1515_core_75p37M"
        genes_ssd1 = ['b1854', 'b1676', 'b1479', 'b2463', 'b1852']
        m.reactions.EX_o2_e.lower_bound = -2.0
        m.reactions.EX_glc__D_e.lower_bound = -10.0
        m.reactions.get_by_id(BIO_I1515).bounds = (0, 0)
        m.objective = 'EX_succ_e'
        y_max_wt = m.slim_optimize(error_value=float('nan')) / 10.0
        m_mut = m.copy()
        for g in genes_ssd1:
            if g in m_mut.genes: m_mut.genes.get_by_id(g).knock_out()
        y_max_mut = m_mut.slim_optimize(error_value=float('nan')) / 10.0
        print(f"\n--- E4: Condiciones del Paper (Toya 2022) ---")
        print(f"WT  (maximo FBA): {y_max_wt:.4f} mol/mol")
        print(f"SSD1 (maximo FBA): {y_max_mut:.4f} mol/mol")
        assert y_max_wt >= 0.20 - 1e-6
        assert y_max_mut >= 0.48 - 1e-6
        assert y_max_mut == pytest.approx(y_max_wt, abs=0.01)
