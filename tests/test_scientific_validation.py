"""
Validación CIENTÍFICA de MetaboEngineer Pro (V6.3)
"""
import os
import re
import subprocess
import sys
import textwrap
import types

import pytest
from cobra import Metabolite, Reaction, Model
from cobra.flux_analysis import single_gene_deletion
from cobra.flux_analysis import flux_variability_analysis as fva

import metabo_engineer as me

BIOMASA = "Biomass_Ecoli_core"

@pytest.fixture(scope="module")
def textbook():
    modelo, _ = me.cargar_modelo_universal("textbook")
    return modelo

@pytest.fixture
def sin_equilibrator(monkeypatch):
    try:
        import equilibrator_api
    except ImportError:
        stub = types.ModuleType("equilibrator_api")
        stub.ComponentContribution = None
        stub.Q_ = None
        monkeypatch.setitem(sys.modules, "equilibrator_api", stub)
    yield

class TestAnclasExternas:
    def test_crecimiento_aerobio_referencia(self, textbook):
        m = me.configurar_ambiente(textbook, "aerobio")
        assert m.slim_optimize() == pytest.approx(0.8739, abs=1e-3)

    def test_crecimiento_anaerobio_referencia(self, textbook):
        m = me.configurar_ambiente(textbook, "anaerobio")
        assert m.slim_optimize() == pytest.approx(0.2117, abs=1e-3)

    def test_residuo_estequiometrico_y_cotas(self, textbook):
        sol = textbook.optimize()
        S = __import__("cobra").util.array.create_stoichiometric_matrix(textbook)
        v = sol.fluxes.loc[[r.id for r in textbook.reactions]].values
        assert abs(S @ v).max() < 1e-7
        for r in textbook.reactions:
            assert r.lower_bound - 1e-7 <= sol.fluxes[r.id] <= r.upper_bound + 1e-7

    def test_modelo_balanceado_en_masa_y_carga(self, textbook):
        malas = {r.id: r.check_mass_balance() for r in textbook.reactions
                 if not r.boundary and "biomass" not in r.id.lower()
                 and r.check_mass_balance()}
        assert not malas, malas

    def test_no_hay_energia_gratis(self, textbook):
        with textbook as m:
            for ex in m.exchanges: 
                ex.bounds = (0, 0)
            m.reactions.ATPM.bounds = (0, 1000)
            m.objective = "ATPM"
            assert m.slim_optimize(error_value=0.0) == pytest.approx(0, abs=1e-6)


class TestLogicaGPR:
    def test_esencialidad_coincide_con_cobra(self, textbook):
        m = textbook.copy()
        wt = m.slim_optimize()
        ref = single_gene_deletion(m)
        esperado = {next(iter(i)) for i, g in zip(ref.ids, ref.growth) if g > 0.05 * wt}
        
        # CORRECCIÓN: Excluimos el pseudogen s0001 de la comparación
        esperado.discard("s0001") 
        
        obtenido = set(me.obtener_genes_no_esenciales(m, BIOMASA))
        assert obtenido == esperado

    def test_isoenzima_acnA_acnB_viables_por_separado(self, textbook):
        for gen in ("b1276", "b0118"):
            with textbook as m:
                m.genes.get_by_id(gen).knock_out()
                assert m.slim_optimize(error_value=0.0) > 0.8

    def test_evaluar_bacteria_respeta_GPR(self, textbook):
        fit = me.evaluar_bacteria(textbook, ["b1276"], BIOMASA, BIOMASA)
        assert fit > 0.3

    def test_cohorte_knockouts_viables_no_dan_fitness_cero(self, textbook):
        fallos = []
        for g in textbook.genes:
            with textbook as mm:
                mm.genes.get_by_id(g.id).knock_out()
                real = mm.slim_optimize(error_value=0.0)
            if real >= 0.5 and len(g.reactions) <= 2:
                if me.evaluar_bacteria(textbook, [g.id], BIOMASA, BIOMASA) <= 0:
                    fallos.append(g.id)
        assert not fallos, f"Viables en realidad pero fitness=0: {fallos}"

    def test_knockout_fitness_positivo_o_cero(self, textbook):
        import random
        m = textbook.copy()
        m.reactions.EX_succ_e.bounds = (0, 1000)
        ids = [g.id for g in m.genes]
        rnd = random.Random(0)
        for _ in range(30):
            ko = rnd.sample(ids, rnd.randint(1, 3))
            fit = me.evaluar_bacteria(m, ko, "EX_succ_e", BIOMASA)
            assert fit >= 0.0

    def test_cepa_acoplada_conocida_supera_al_WT(self, textbook):
        m = textbook.copy()
        m.reactions.EX_succ_e.bounds = (0, 1000)
        cepa = ["b3733", "b3870", "s0001"]
        
        with m as mm:
            for g in cepa: 
                mm.genes.get_by_id(g).knock_out()
            mm.objective = BIOMASA
            garantizado = fva(mm, ["EX_succ_e"], fraction_of_optimum=1.0).minimum.iloc[0]
            
        with m as mm:
            mm.objective = BIOMASA
            wt_garantizado = fva(mm, ["EX_succ_e"], fraction_of_optimum=1.0).minimum.iloc[0]
            
        assert garantizado > wt_garantizado + 1.0
        assert (me.evaluar_bacteria(m, cepa, "EX_succ_e", BIOMASA) > me.evaluar_bacteria(m, [], "EX_succ_e", BIOMASA))


class TestOptimizacion:
    def test_fitness_sin_knockouts_igual_a_fba_al_90_por_ciento(self, textbook):
        m = textbook.copy()
        m.reactions.EX_succ_e.bounds = (0, 1000)
        fit = me.evaluar_bacteria(m, [], "EX_succ_e", BIOMASA)
        
        with m as mm:
            mm.objective = BIOMASA
            max_bio = mm.slim_optimize()
            mm.reactions.get_by_id(BIOMASA).lower_bound = max_bio * 0.9
            mm.objective = "EX_succ_e"
            ref = mm.slim_optimize()
            
        assert fit == pytest.approx(ref, rel=1e-6)

    def test_optknock_reproducible_con_misma_semilla(self):
        codigo = textwrap.dedent("""
            import warnings; warnings.filterwarnings("ignore")
            import sys, types, random, io, contextlib
            with contextlib.redirect_stdout(io.StringIO()):
                import metabo_engineer as me, cli_experto as cli
            m,_ = me.cargar_modelo_universal("textbook")
            a = types.SimpleNamespace(poblacion=15, generaciones=5, mutaciones=3)
            random.seed(7)
            with contextlib.redirect_stdout(io.StringIO()):
                b = cli.metodo_optknock(m, "EX_succ_e", a)
            print(sorted(b["genes"]), round(b["fitness"], 4))
        """)
        salidas = set()
        for hs in ("1", "2"):
            env = dict(os.environ, PYTHONHASHSEED=hs)
            r = subprocess.run([sys.executable, "-c", codigo], capture_output=True, text=True, env=env, timeout=300)
            salidas.add(r.stdout.strip().splitlines()[-1])
        assert len(salidas) == 1, salidas


def _parsear_formula(f):
    izq, der = f.split(" = ")
    out = {}
    for lado, signo in ((izq, -1), (der, 1)):
        for t in lado.split(" + "):
            mt = re.match(r"^(?:(\d+(?:\.\d+)?) )?(?:bigg|kegg):(\S+)$", t.strip())
            out[mt.group(2)] = signo * float(mt.group(1) or 1)
    return out


class TestQuimica:
    def test_construir_reaccion_conserva_estequiometria(self, sin_equilibrator):
        from modulo_termodinamica_v2 import construir_reaccion
        mod = Model("t")
        h2, o2, h2o = (Metabolite(i, formula=f, compartment="c") for i, f in (("h2_c", "H2"), ("o2_c", "O2"), ("h2o_c", "H2O")))
        r = Reaction("R")
        mod.add_reactions([r])
        r.add_metabolites({h2: -1.0, o2: -0.5, h2o: 1.0})
        got = _parsear_formula(construir_reaccion(mod, r))
        assert got == {"h2_c": -1.0, "o2_c": -0.5, "h2o_c": 1.0}


class TestRutaCafeina:
    def test_inyeccion_no_crashea(self, textbook):
        mm, q = me.inyectar_ruta_heterologa(textbook.copy(), "EX_caft_e")
        assert mm is not None and q == "EX_caft_e"

    def test_ruta_cafeina_balanceada_en_masa(self, textbook):
        mm, _ = me.inyectar_ruta_heterologa(textbook.copy(), "EX_caft_e")
        formulas = {
            "xan_c": "C5H4N4O2", "tbrm_c": "C7H8N4O2", "caft_c": "C8H10N4O2",
            "amet_c": "C15H23N6O5S", "ahcys_c": "C14H20N6O5S", "h_c": "H"
        }
        for mid, f in formulas.items():
            if mid in mm.metabolites:
                mm.metabolites.get_by_id(mid).formula = f
        
        for rid in ("CAFT_RXN1", "CAFT_RXN2"):
            assert not mm.reactions.get_by_id(rid).check_mass_balance(), rid


class TestReadmeVsCodigo:
    @pytest.mark.xfail(strict=True, reason="README promete Rhea, pero el código devuelve N/A.")
    def test_rhea_id_realmente_se_consulta(self, sin_equilibrator, monkeypatch, tmp_path):
        import reporte_investigador as ri
        monkeypatch.chdir(tmp_path)
        class R:
            status_code = 200
            def json(self):
                return {"results": [{"genes": [{"geneName": {"value": "adhE"}}], "proteinDescription": {"recommendedName": {"fullName": {"value": "Aldehyde-alcohol dehydrogenase"}}}}]}
        monkeypatch.setattr(ri.requests, "get", lambda *a, **k: R())
        monkeypatch.setattr(ri.time, "sleep", lambda *_: None)
        _, _, rhea = ri.buscar_info_uniprot_rhea("b1241")
        assert rhea.startswith("RHEA:")

    @pytest.mark.xfail(strict=True, reason="configurar_ambiente no avisa si el modelo no tiene EX_o2_e.")
    def test_anaerobio_en_modelo_sin_o2_debe_avisar(self):
        with pytest.raises(Exception):
            me.configurar_ambiente(Model("sin_o2"), "anaerobio")