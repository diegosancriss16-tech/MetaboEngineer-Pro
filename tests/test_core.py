"""
Tests de funcionalidad core de MetaboEngineer Pro.

Cubre:
- Carga de modelos (éxito y fallo)
- Evaluación de bacterias (fitness)
- Detección de organismo (UniProt taxon ID)
- Hashing SHA256
"""

import pytest
import metabo_engineer as me


# ==========================================
# CARGA DE MODELOS
# ==========================================
class TestModelLoading:
    def test_cargar_modelo_textbook_ok(self):
        """textbook es local y siempre debe cargar."""
        modelo, msg = me.cargar_modelo_universal("textbook")
        assert modelo is not None
        assert len(modelo.reactions) > 0
        assert len(modelo.genes) > 0

    def test_cargar_modelo_string_vacio_falla(self):
        """String vacío debe lanzar excepción clara."""
        with pytest.raises(Exception) as exc:
            me.cargar_modelo_universal("")
        assert "vacio" in str(exc.value).lower() or "invalido" in str(exc.value).lower()

    def test_cargar_modelo_none_falla(self):
        """None debe lanzar excepción (no romper con AttributeError)."""
        with pytest.raises(Exception):
            me.cargar_modelo_universal(None)

    def test_cargar_modelo_espacios_falla(self):
        """Solo espacios debe tratarse como vacío."""
        with pytest.raises(Exception):
            me.cargar_modelo_universal("   ")


# ==========================================
# EVALUACIÓN DE BACTERIAS (FITNESS)
# ==========================================
class TestEvaluarBacteria:
    def test_evaluar_sin_genes_devuelve_fitness(self, modelo_textbook, gen_biomasa_textbook):
        """Sin knockouts, debe devolver un fitness > 0."""
        fitness = me.evaluar_bacteria(modelo_textbook, [], "EX_succ_e", gen_biomasa_textbook)
        assert fitness >= 0.0

    def test_evaluar_knockout_extremo_no_crashea(self, modelo_textbook, gen_biomasa_textbook):
        """Bloquear muchos genes no debe crashear, debe devolver >= 0."""
        todos_los_genes = [g.id for g in modelo_textbook.genes[:20]]
        fitness = me.evaluar_bacteria(modelo_textbook, todos_los_genes, "EX_succ_e", gen_biomasa_textbook)
        assert fitness >= 0.0

    def test_evaluar_gen_inexistente_no_crashea(self, modelo_textbook, gen_biomasa_textbook):
        """Gen que no existe en el modelo debe ignorarse sin crashear."""
        fitness = me.evaluar_bacteria(modelo_textbook, ["gen_falso_123"], "EX_succ_e", gen_biomasa_textbook)
        assert fitness >= 0.0


# ==========================================
# DETECCIÓN DE ORGANISMO (taxon UniProt)
# ==========================================
class TestOrganismId:
    @pytest.mark.parametrize("modelo_id,expected", [
        ("iYO844", "224308"),      # B. subtilis
        ("iCW773", "196627"),      # C. glutamicum
        ("iMM904", "559292"),      # S. cerevisiae
        ("iEK1008", "1773"),       # M. tuberculosis
        ("textbook", "83333"),     # E. coli
        ("iML1515", "83333"),      # E. coli
        ("desconocido", "83333"),  # fallback E. coli
        ("", "83333"),             # vacío → fallback
        (None, "83333"),           # None → fallback
    ])
    def test_get_organism_id(self, modelo_id, expected):
        assert me.get_organism_id(modelo_id) == expected


# ==========================================
# HASHING SHA256
# ==========================================
class TestSha256:
    def test_sha256_archivo_mismo_contenido(self, tmp_path):
        """Dos archivos con mismo contenido → mismo hash."""
        a1 = tmp_path / "a.txt"
        a1.write_text("hola mundo")
        a2 = tmp_path / "b.txt"
        a2.write_text("hola mundo")

        h1 = me._sha256_archivo(str(a1))
        h2 = me._sha256_archivo(str(a2))

        assert h1 == h2
        assert len(h1) == 64  # SHA256 en hex = 64 chars

    def test_sha256_archivo_distinto_contenido(self, tmp_path):
        """Contenidos distintos → hashes distintos."""
        a1 = tmp_path / "a.txt"
        a1.write_text("hola")
        a2 = tmp_path / "b.txt"
        a2.write_text("adios")

        assert me._sha256_archivo(str(a1)) != me._sha256_archivo(str(a2))