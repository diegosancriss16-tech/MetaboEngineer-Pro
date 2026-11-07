"""El GA debe ser reproducible aunque otras librerias consuman el `random` global
(optlang/uuid.uuid1 lo hace en Python <= 3.12)."""

import random
import uuid

import metabo_engineer as me

POOL = [f"g{i}" for i in range(60)]


def _fitness_determinista(genes):
    return float(sum(int(g[1:]) for g in genes) % 17)


def _correr(semilla, contaminar):
    random.seed(semilla)

    def fitness(genes):
        if contaminar:
            random.getrandbits(14)      # lo que hace uuid.uuid1() en Python <= 3.12
            uuid.uuid1()
            random.random()
        return _fitness_determinista(genes)

    return me._evolucionar(POOL, 3, 4, 12, fitness, crossover_split=2)


def test_misma_semilla_mismo_resultado():
    assert _correr(42, False) == _correr(42, False)


def test_resultado_inmune_al_consumo_del_random_global():
    limpio, contaminado = _correr(42, False), _correr(42, True)
    assert limpio["genes"] == contaminado["genes"]
    assert limpio["fitness"] == contaminado["fitness"]


def test_semillas_distintas_dan_busquedas_distintas():
    visto = {tuple(_correr(s, False)["genes"]) for s in range(8)}
    assert len(visto) > 1


def test_rng_explicito_tiene_prioridad_sobre_el_global():
    random.seed(1)
    a = me._evolucionar(POOL, 3, 3, 10, _fitness_determinista, rng=random.Random(99))
    random.seed(2)
    b = me._evolucionar(POOL, 3, 3, 10, _fitness_determinista, rng=random.Random(99))
    assert a == b
