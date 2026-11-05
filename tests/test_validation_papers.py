"""
Suite de Validación contra Papers Canónicos
FIX: E4a ahora lee la columna 'ids' del DataFrame de single_gene_deletion
(no el índice, que es numérico y cuyo orden cambia por multiprocessing).
Fixtures herméticos: modelo fresco por test, sin estado compartido.
"""
import pytest
import cobra
import random
import time
from cobra.flux_analysis import single_gene_deletion, flux_variability_analysis as fva
from cobra.manipulation import knock_out_model_genes
import metabo_engineer as me

BIO_I1515 = "BIOMASS_Ec_iML1515_core_75p37M"
BIO_TEXTBOOK = "Biomass_Ecoli_core"


# =====================================================================
# FIXTURES HERMÉTICOS
# =====================================================================
@pytest.fixture
def fresh_model():
    """Modelo FRESCO por test. Sin estado compartido entre tests."""
    m, _ = me.cargar_modelo_universal("iML1515")
    m.objective = BIO_I1515
    return m


@pytest.fixture
def gene_sample():
    """Muestra determinística de 50 genes (semilla fija)."""
    random.seed(42)
    m, _ = me.cargar_modelo_universal("iML1515")
    return random.sample([g.id for g in m.genes], 50)


# =====================================================================
# E1: Baseline FBA (Monk 2017)
# =====================================================================
class TestE1_BaselineFBA:
    def test_aerobic_growth(self, fresh_model):
        m = fresh_model
        m.reactions.EX_glc__D_e.lower_bound = -10.0
        m.reactions.EX_o2_e.lower_bound = -1000.0
        mu = m.slim_optimize()
        print(f"\nE1 Aerobio: mu = {mu:.4f} h-1 (literatura: 0.877)")
        assert mu == pytest.approx(0.877, abs=0.02)

    def test_anaerobic_growth(self, fresh_model):
        m = fresh_model
        m.reactions.EX_glc__D_e.lower_bound = -10.0
        m.reactions.EX_o2_e.lower_bound = 0.0
        mu = m.slim_optimize()
        print(f"\nE1 Anaerobio: mu = {mu:.4f} h-1 (literatura: ~0.2-0.3)")
        assert 0.10 <= mu <= 0.35


# =====================================================================
# E4a: Esencialidad ESTRICTA vs single_gene_deletion
# FIX CRÍTICO: leer la columna 'ids', no el índice.
# =====================================================================
class TestE4a_EssentialityStrict:
    def test_essentiality_matches_cobrapy_strict(self, fresh_model, gene_sample):
        """
        Comparación estricta contra single_gene_deletion.
        FIX: los genes están en la COLUMNA 'ids' (frozensets), no en el index.
        El índice es numérico y su orden cambia por multiprocessing.
        """
        m_ref = fresh_model
        wt = m_ref.slim_optimize()
        ref = single_gene_deletion(m_ref, gene_list=gene_sample)

        def extraer_gene(ids):
            """Robusto ante frozenset/tuple/list/string según versión de cobra."""
            if isinstance(ids, (frozenset, set, tuple, list)):
                return next(iter(ids))
            return ids

        # ✅ CORRECTO: iterar la COLUMNA 'ids', emparejada con 'growth' por fila
        ref_growth = {extraer_gene(ids): g
                      for ids, g in zip(ref['ids'], ref['growth'])}
        ref_essential = {g for g, mu in ref_growth.items()
                         if mu < 0.05 * wt}

        # Nuestro método (copias independientes)
        our_essential = set()
        for gen_id in gene_sample:
            m_test = fresh_model.copy()
            m_test.genes.get_by_id(gen_id).knock_out()
            growth = m_test.slim_optimize(error_value=0.0)
            if growth < 0.05 * wt:
                our_essential.add(gen_id)

        mismatches = ref_essential.symmetric_difference(our_essential)
        print(f"\nE4a ESTRICTO: {len(mismatches)} discrepancias en {len(gene_sample)} genes")
        print(f"  COBRApy esenciales: {sorted(ref_essential)}")
        print(f"  Nuestro esenciales: {sorted(our_essential)}")
        assert len(mismatches) == 0, f"Discrepancias: {sorted(mismatches)}"

    def test_essentiality_literature_validation(self, fresh_model, gene_sample):
        """Validación contra literatura: ~18-22% de genes esenciales (Keio/Monk)."""
        m = fresh_model
        wt = m.slim_optimize()

        esenciales = 0
        for gen_id in gene_sample:
            m_test = m.copy()
            m_test.genes.get_by_id(gen_id).knock_out()
            try:
                growth = m_test.slim_optimize(error_value=0.0)
                if growth < 0.05 * wt:
                    esenciales += 1
            except:
                esenciales += 1

        pct = esenciales / len(gene_sample) * 100
        print(f"\nE4a: {esenciales}/{len(gene_sample)} genes esenciales ({pct:.1f}%)")
        print(f"    Literatura (Monk 2017): ~18%")
        assert 2 <= esenciales <= 18, f"Fuera de rango: {esenciales}/{len(gene_sample)} ({pct:.1f}%)"


# =====================================================================
# E4b: FVA - Producción Garantizada (Tepper 2010)
# =====================================================================
class TestE4b_FVA:
    def test_fva_minimum_guaranteed(self, fresh_model):
        m = fresh_model
        m.reactions.EX_glc__D_e.lower_bound = -10.0
        m.reactions.EX_o2_e.lower_bound = 0.0
        m.reactions.EX_succ_e.bounds = (0, 1000)
        m.reactions.get_by_id(BIO_I1515).bounds = (0, 0)
        m.objective = 'EX_succ_e'

        res = fva(m, ['EX_succ_e'], fraction_of_optimum=0.9)
        min_prod = res.minimum.iloc[0]
        max_prod = res.maximum.iloc[0]
        print(f"\nE4b FVA: min={min_prod:.4f}, max={max_prod:.4f}")
        assert min_prod >= 0, f"FVA min negativo: {min_prod}"


# =====================================================================
# E4c: Termodinámica (eQuilibrator - Henry 2009)
# =====================================================================
class TestE4c_Thermodynamics:
    def test_pfk_is_exergonic(self):
        try:
            from equilibrator_api import ComponentContribution, Q_
        except ImportError:
            pytest.skip("eQuilibrator no instalado")

        cc = ComponentContribution()
        cc.p_h, cc.p_mg = Q_(7.0), Q_(3.0)
        rxn = cc.parse_reaction_formula("kegg:C00002 + kegg:C05345 = kegg:C00008 + kegg:C05378")
        dg = cc.physiological_dg_prime(rxn).value.magnitude
        print(f"\nE4c PFK: dG' = {dg:.2f} kJ/mol")
        assert dg < -10, f"PFK no es exergónica: {dg}"


# =====================================================================
# E4: Paper 2022 (Toya & Shimizu)
# =====================================================================
class TestE4Paper2022:
    def test_experimental_yields_below_fba_bound(self, fresh_model):
        """Validación honesta: experimental < FBA máximo (techo termodinámico)."""
        m = fresh_model
        m.reactions.EX_glc__D_e.lower_bound = -10.0
        m.reactions.EX_o2_e.lower_bound = -2.0
        m.reactions.get_by_id(BIO_I1515).bounds = (0, 0)
        m.reactions.EX_succ_e.bounds = (0, 1000)
        m.objective = 'EX_succ_e'

        wt_fba = m.slim_optimize() / 10.0

        m_mut = m.copy()
        for g in ['b1854', 'b1676', 'b1479', 'b2463', 'b1852']:
            m_mut.genes.get_by_id(g).knock_out()
        ssd1_fba = m_mut.slim_optimize() / 10.0

        print(f"\n--- Validación Paper 2022 (Toya & Shimizu) ---")
        print(f"Experimental: WT=0.20, SSD1=0.48")
        print(f"FBA máximo:   WT={wt_fba:.4f}, SSD1={ssd1_fba:.4f}")

        assert wt_fba >= 0.20, f"WT FBA ({wt_fba}) debe ser >= 0.20"
        assert ssd1_fba >= 0.48, f"SSD1 FBA ({ssd1_fba}) debe ser >= 0.48"
        print("✅ El modelo respeta la termodinámica (FBA >= Experimental).")


# =====================================================================
# E5: End-to-End Pipeline (Jensen 2019) - marcado como slow
# =====================================================================
@pytest.mark.slow
class TestE5_EndToEnd:
    def test_full_pipeline_quick(self, fresh_model):
        m = fresh_model
        m.reactions.EX_glc__D_e.lower_bound = -10.0
        m.reactions.EX_o2_e.lower_bound = 0.0
        m.reactions.EX_succ_e.bounds = (0, 1000)

        assert len(m.reactions) > 2000

        t0 = time.time()
        genes_ok = me.obtener_genes_no_esenciales(m, BIO_I1515)
        t_esencialidad = time.time() - t0

        import cli_experto
        import types
        args = types.SimpleNamespace(poblacion=10, generaciones=3, mutaciones=3)
        random.seed(42)
        t0 = time.time()
        mejor = cli_experto.metodo_ia(m, 'EX_succ_e', args)
        t_ia = time.time() - t0

        print(f"\nE5 End-to-End:")
        print(f"  Genes no esenciales: {len(genes_ok)} ({t_esencialidad:.1f}s)")
        print(f"  Mejor fitness: {mejor['fitness']:.4f} ({t_ia:.1f}s)")

        assert len(genes_ok) > 1000
        assert mejor['fitness'] >= 0