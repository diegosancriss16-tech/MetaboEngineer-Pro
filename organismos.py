# -*- coding: utf-8 -*-
"""Mapeo de IDs de modelo a taxon de UniProt (modulo comun sin dependencias pesadas)."""

_TAXON_POR_MODELO = {
    "iYO844": "224308",   # Bacillus subtilis 168
    "iCW773": "196627",   # Corynebacterium glutamicum
    "iMM904": "559292",   # Saccharomyces cerevisiae S288C
    "iEK1008": "1773",    # Mycobacterium tuberculosis
}
_TAXON_DEFECTO = "83333"  # Escherichia coli K-12


def get_organism_id(modelo_id: str) -> str:
    if not modelo_id:
        return _TAXON_DEFECTO
    for clave, taxon in _TAXON_POR_MODELO.items():
        if clave in modelo_id:
            return taxon
    return _TAXON_DEFECTO
