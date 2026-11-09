import pytest
import metabo_engineer as me

class TestEvaluarBacteria:
    def test_evaluar_lista_vacia_lanza_error(self, modelo_textbook, gen_biomasa_textbook):
        """Una lista vacía debe lanzar ValueError explícito."""
        with pytest.raises(ValueError):
            me.evaluar_bacteria(modelo_textbook, [], "EX_succ_e", gen_biomasa_textbook)

    def test_evaluar_gen_inexistente_lanza_error(self, modelo_textbook, gen_biomasa_textbook):
        """Un gen que no existe en el modelo debe lanzar ValueError."""
        with pytest.raises(ValueError):
            me.evaluar_bacteria(modelo_textbook, ["gen_falso_123"], "EX_succ_e", gen_biomasa_textbook)

    def test_evaluar_knockout_extremo_no_crashea(self, modelo_textbook, gen_biomasa_textbook):
        """Bloquear 20 genes válidos no debe crashear; devuelve -1e6 (letal) o >= 0."""
        todos_los_genes = [g.id for g in modelo_textbook.genes[:20]]
        fitness = me.evaluar_bacteria(modelo_textbook, todos_los_genes, "EX_succ_e", gen_biomasa_textbook)
        assert fitness == -1e6 or fitness >= 0.0