"""Tests de construir_reaccion: sin sustratos inventados y con IDs BiGG sin compartimento."""

import pytest
from cobra import Metabolite, Model, Reaction

pytest.importorskip("equilibrator_api")  # modulo_termodinamica_v2 la importa al cargarse

from modulo_termodinamica_v2 import _id_bigg_sin_compartimento, construir_reaccion


@pytest.mark.parametrize("met_id, esperado", [
    ("atp_c", "atp"), ("h2o_e", "h2o"), ("pe160_p", "pe160"), ("4fe4s_c", "4fe4s"),
    ("glc__D_e", "glc__D"), ("nadhx__R_c", "nadhx__R"), ("cur", "cur"),
])
def test_id_bigg_sin_compartimento(met_id, esperado):
    assert _id_bigg_sin_compartimento(met_id) == esperado


def _modelo_con_metabolitos_sin_kegg():
    m = Model("mini")
    a = Metabolite("foo_c", compartment="c")
    b = Metabolite("bar_c", compartment="c")
    r = Reaction("R1")
    m.add_reactions([r])
    r.add_metabolites({a: -1, b: 1})
    return m, r


def test_metabolito_sin_kegg_usa_namespace_bigg_metabolite_sin_sufijo():
    m, r = _modelo_con_metabolitos_sin_kegg()
    formula = construir_reaccion(m, r)
    assert formula == "bigg.metabolite:foo = bigg.metabolite:bar"
    assert "bigg:" not in formula and "h2o" not in formula


def test_reaccion_de_intercambio_no_se_completa_con_agua(modelo_textbook):
    ex = modelo_textbook.reactions.get_by_id("EX_succ_e")
    with pytest.raises(ValueError, match="frontera"):
        construir_reaccion(modelo_textbook, ex)


def test_reaccion_de_un_solo_lado_lanza_error():
    m = Model("mini")
    a = Metabolite("foo_c", compartment="c")
    r = Reaction("SOLO_SUSTRATO")
    m.add_reactions([r])
    r.add_metabolites({a: -1})
    with pytest.raises(ValueError):
        construir_reaccion(m, r)


def test_reaccion_normal_de_textbook_usa_kegg(modelo_textbook):
    formula = construir_reaccion(modelo_textbook, modelo_textbook.reactions.get_by_id("PGI"))
    assert formula.count("kegg:") == 2 and "=" in formula


class _EquilibratorFalso:
    """Sustituto sin red: cualquier calculo de DG' falla, como con un compuesto no resuelto."""
    p_h = p_mg = ionic_strength = None

    def parse_reaction_formula(self, formula):
        raise RuntimeError("compuesto no resuelto")


def test_reporte_marca_intercambio_y_no_inventa_dg(modelo_textbook, monkeypatch, tmp_path):
    import reporte_investigador as rep
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(rep, "ComponentContribution", _EquilibratorFalso)
    monkeypatch.setattr(rep, "buscar_info_uniprot_rhea", lambda g, o="83333": (g, "n/d", "N/A"))
    rep.generar_documento_markdown(modelo_textbook, "textbook", "EX_succ_e", ["b1818"], 1.0)
    texto = next((tmp_path / "reportes").glob("REPORTE_*.md")).read_text(encoding="utf-8")
    fila = next(l for l in texto.splitlines() if l.startswith("| `EX_succ_e`"))
    assert "Intercambio" in fila
    # Las reacciones internas que no se pueden calcular aparecen como N/D, no se omiten
    assert "N/D" in texto
