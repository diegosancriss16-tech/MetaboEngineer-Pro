"""
Tests del cache de UniProt offline.

Cubre:
- Guardar/cargar cache
- Limpiar cache
- Traducción de genes comunes (sin API)
"""

import json
import pytest
import metabo_engineer as me


class TestUniProtCache:
    def test_cache_vacio_inicial(self, tmp_path, monkeypatch):
        """Si no existe el archivo, devuelve dict vacío."""
        cache_file = tmp_path / "cache_nuevo.json"
        monkeypatch.setattr(me, "_RUTA_CACHE_UNIPROT", str(cache_file))

        cache = me._cargar_cache_uniprot()
        assert cache == {}

    def test_cache_roundtrip(self, tmp_path, monkeypatch):
        """Guardar y recargar debe preservar los datos."""
        cache_file = tmp_path / "cache_test.json"
        monkeypatch.setattr(me, "_RUTA_CACHE_UNIPROT", str(cache_file))

        cache = {"83333:b0001": ["thrA", "Bifunctional aspartokinase"]}
        me._guardar_cache_uniprot(cache)

        # Verificar que el archivo existe y tiene JSON válido
        assert cache_file.exists()
        with open(cache_file) as f:
            data = json.load(f)
        assert data == cache

        # Recargar con la función
        cache2 = me._cargar_cache_uniprot()
        assert cache2 == cache

    def test_cache_corrupto_devuelve_vacio(self, tmp_path, monkeypatch):
        """Si el JSON está corrupto, no debe crashear → devuelve {}."""
        cache_file = tmp_path / "cache_corrupto.json"
        cache_file.write_text("{ ESTO NO ES JSON VALIDO")
        monkeypatch.setattr(me, "_RUTA_CACHE_UNIPROT", str(cache_file))

        cache = me._cargar_cache_uniprot()
        assert cache == {}

    def test_limpiar_cache_existente(self, tmp_path, monkeypatch):
        """limpiar_cache_uniprot debe borrar el archivo."""
        cache_file = tmp_path / "cache_borrar.json"
        cache_file.write_text("{}")
        monkeypatch.setattr(me, "_RUTA_CACHE_UNIPROT", str(cache_file))

        me.limpiar_cache_uniprot()
        assert not cache_file.exists()

    def test_limpiar_cache_inexistente_no_crashea(self, tmp_path, monkeypatch):
        """Si no hay cache, no debe crashear."""
        cache_file = tmp_path / "no_existe.json"
        monkeypatch.setattr(me, "_RUTA_CACHE_UNIPROT", str(cache_file))

        # No debe lanzar excepción
        me.limpiar_cache_uniprot()


class TestTraducirGen:
    def test_gen_sin_digitos_devuelve_nombre_comun(self):
        """Genes sin dígitos (pgk, adhE) se consideran nombres comunes."""
        nombre, funcion = me.traducir_gen_uniprot("pgk", "83333")
        assert nombre == "pgk"
        assert funcion == "Nombre comun"

    def test_gen_con_digitos_va_a_cache_o_api(self):
        """Genes con dígitos (b0001) pasan por cache/API (no son 'Nombre comun')."""
        # No verificamos el resultado exacto porque depende de UniProt.
        # Pero sí que NO devuelva "Nombre comun" para genes tipo bXXXX
        nombre, funcion = me.traducir_gen_uniprot("b0001", "83333")
        assert funcion != "Nombre comun"