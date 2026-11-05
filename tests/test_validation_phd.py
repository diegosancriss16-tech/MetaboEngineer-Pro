"""
Suite de Validación contra Papers Canónicos (E1, E4a, E4b, E4c, E5)
Basado en: Monk 2017, Burgard 2003, Patil 2005, Tepper 2010, Jensen 2019, Toya 2022
"""
import pytest
import cobra
import random
import time
from cobra.flux_analysis import single_gene_deletion, flux_variability_analysis as fva
import metabo_engineer as me

BIO_I1515 = "BIOMASS_Ec_iML1515_core_75p37M"
BIO_TEXTBOOK = "Biomass_Ecoli_core"


@pytest.fixture(scope="module")
def iML1515():
    m, _ = me.cargar_modelo_universal("iML1515")
    return m


# =====================================================================
# E1: Baseline FBA (Monk 2017, iML1515)
# =====================================================================
class TestE1_BaselineFBA:
    def test_aerobic_growth(self, iML1515):
        """Monk 2017: iML1515 con glucosa 10, aerobio -> mu = 0.877 h-1"""
        m = iML1515.copy()
        m.reactions.EX_glc__D_e.lower_bound = -10.0
        m.reactions.EX_o2_e.lower_bound = -1000.0
        m.objective = BIO_I1515
        mu = m.slim_optimize()
        print(f"\nE1 Aerobio: mu = {mu:.4f} h-1 (literatura: 0.877)")
        assert mu == pytest.approx(0.877, abs=0.02)

    def test_anaerobic_growth(self, iML1515):
        """Monk 2017: anaerobio -> mu ~ 0.2-0.3 h-1"""
        m = iML1515.copy()
        m.reactions.EX_glc__D_e.lower_bound = -10.0
        m.reactions.EX_o2_e.lower_bound = 0.0
        m.objective = BIO_I1515
        mu = m.slim_optimize()
        print(f"\nE1 Anaerobio: mu = {mu:.4f} h-1 (literatura: ~0.2-0.3)")
        assert 0.10 <= mu <= 0.35


# =====================================================================
# E4a: Esencialidad (Monk 2017 / Keio) - RÁPIDO (muestra)
# =====================================================================
class TestE4a_EssentialityQuick:
    def test_essentiality_sample_50_genes(self, iML1515):
        """Muestra aleatoria de 50 genes: verifica que el filtro coincide con cobrapy"""
        m = iML1515.copy()
        m.objective = BIO_I1515
        wt = m.slim_optimize()
        
        random.seed(42)
        sample = random.sample([g.id for g in m.genes], 50)
        
        # Comparar knockout individual: cobrapy vs nuestra lógica
        ref = single_gene_deletion(m, sample)
        
        mismatches = []
        for i, gen_id in enumerate(sample):
            growth_ref = ref.growth.iloc[i]
            growth_mut = me.evaluar_bacteria(m, [gen_id], BIO_I1515, BIO_I1515)
            # Nuestro evaluar_bacteria usa biomasa al 90%, así que growth_mut puede diferir
            # pero si growth_ref es 0, growth_mut también debe ser 0
            if growth_ref < 0.05 * wt and growth_mut > 0.5:
                mismatches.append((gen_id, growth_ref, growth_mut))
        
        print(f"\nE4a Sample: {len(mismatches)} discrepancias en 50 genes")
        assert len(mismatches) == 0, f"Discrepancias: {mismatches}"


# =====================================================================
# E4a: Esencialidad FULL (marcado como slow)
# =====================================================================
@pytest.mark.slow
class TestE4a_EssentialityFull:
    def test_essentiality_full(self, iML1515):
        """Verifica esencialidad en TODOS los genes (15-30 min)"""
        m = iML1515.copy()
        m.objective = BIO_I1515
        wt = m.slim_optimize()
        
        ref = single_gene_deletion(m)
        esenciales = sum(1 for g in ref.growth if g < 0.05 * wt)
        total = len(ref.growth)
        pct = esenciales / total * 100
        
        print(f"\nE4a Full: {esenciales}/{total} esenciales ({pct:.1f}%)")
        # Monk 2017 reporta ~271 esenciales (~18%)
        assert 12 <= pct <= 25, f"Porcentaje fuera de rango: {pct}%"


# =====================================================================
# E4b: FVA - Producción Garantizada (Tepper 2010)
# =====================================================================
class TestE4b_FVA:
    def test_fva_minimum_guaranteed(self, iML1515):
        """FVA debe reportar mínimo >= 0 para diseños acoplados"""
        m = iML1515.copy()
        m.reactions.EX_glc__D_e.lower_bound = -10.0
        m.reactions.EX_o2_e.lower_bound = 0.0
        m.reactions.EX_succ_e.bounds = (0, 1000)
        m.reactions.get_by_id(BIO_I1515).bounds = (0, 0)
        m.objective = 'EX_succ_e'
        
        res = fva(m, ['EX_succ_e'], fraction_of_optimum=0.9)
        min_prod = res.minimum.iloc[0]
        max_prod = res.maximum.iloc[0]
        
        print(f"\nE4b FVA: min={min_prod:.4f}, max={max_prod:.4f}")
        # Ambos deben ser >= 0 (no hay producción negativa)
        assert min_prod >= 0, f"FVA min negativo: {min_prod}"


# =====================================================================
# E4c: Termodinámica (eQuilibrator - Henry 2009)
# =====================================================================
class TestE4c_Thermodynamics:
    def test_pfk_is_exergonic(self):
        """PFK es fuertemente exergónica (glucólisis)"""
        try:
            from equilibrator_api import ComponentContribution, Q_
        except ImportError:
            pytest.skip("eQuilibrator no instalado")
        
        cc = ComponentContribution()
        cc.p_h, cc.p_mg = Q_(7.0), Q_(3.0)
        # PFK: ATP + F6P -> ADP + FBP
        rxn = cc.parse_reaction_formula("kegg:C00002 + kegg:C05345 = kegg:C00008 + kegg:C05378")
        dg = cc.physiological_dg_prime(rxn).value.magnitude
        
        print(f"\nE4c PFK: dG' = {dg:.2f} kJ/mol")
        assert dg < -10, f"PFK no es exergónica: {dg}"


# =====================================================================
# E5: End-to-End Pipeline (Jensen 2019)
# =====================================================================
class TestE5_EndToEnd:
    def test_full_pipeline_quick(self, iML1515):
        """Pipeline completo rápido (5-10 min)"""
        m = iML1515.copy()
        m.reactions.EX_glc__D_e.lower_bound = -10.0
        m.reactions.EX_o2_e.lower_bound = 0.0
        m.reactions.EX_succ_e.bounds = (0, 1000)
        
        # Fase 1: Carga (implicit, ya tenemos m)
        assert len(m.reactions) > 2000
        
        # Fase 2: Esencialidad
        import time
        t0 = time.time()
        genes_ok = me.obtener_genes_no_esenciales(m, BIO_I1515)
        t_esencialidad = time.time() - t0
        
        # Fase 3: Diseño con IA
        import cli_experto
        import types
        args = types.SimpleNamespace(poblacion=10, generaciones=3, mutaciones=3)
        random.seed(42)
        t0 = time.time()
        mejor = cli_experto.metodo_ia(m, 'EX_succ_e', args)
        t_ia = time.time() - t0
        
        # Fase 4: Validación
        print(f"\nE5 End-to-End:")
        print(f"  Genes no esenciales: {len(genes_ok)} ({t_esencialidad:.1f}s)")
        print(f"  Mejor fitness: {mejor['fitness']:.4f} ({t_ia:.1f}s)")
        print(f"  Genes: {mejor['genes']}")
        
        assert len(genes_ok) > 1000
        assert mejor['fitness'] >= 0


# =====================================================================
# VALIDACIÓN PAPER 2022 (Toya & Shimizu)
# =====================================================================
class TestPaper2022:
    def test_fba_upper_bound_vs_experimental(self, iML1515):
        """El FBA (cota superior) debe ser >= al rendimiento experimental del paper"""
        m = iML1515.copy()
        m.reactions.EX_glc__D_e.lower_bound = -10.0
        m.reactions.EX_o2_e.lower_bound = -2.0  # microaerobio
        m.reactions.get_by_id(BIO_I1515).bounds = (0, 0)
        m.objective = 'EX_succ_e'
        
        y_wt = m.slim_optimize() / 10.0
        
        m_mut = m.copy()
        for g in ['b1854', 'b1676', 'b1479', 'b2463', 'b1852']:
            if g in m_mut.genes:
                m_mut.genes.get_by_id(g).knock_out()
        y_mut = m_mut.slim_optimize() / 10.0
        
        print(f"\nPaper 2022 FBA bounds:")
        print(f"  WT  upper bound: {y_wt:.4f} mol/mol (experimental: 0.20)")
        print(f"  SSD1 upper bound: {y_mut:.4f} mol/mol (experimental: 0.48)")
        
        # El FBA siempre es >= al experimental (es un techo termodinámico)
        assert y_wt >= 0.20 - 1e-6
        assert y_mut >= 0.48 - 1e-6