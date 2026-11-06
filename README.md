> **Required Notice:** Copyright Diego San Cristobal (https://github.com/diegosancriss16-tech/MetaboEngineer-Pro)

# 🧬 MetaboEngineer Pro (v0.1.0-beta)

> Suite integral de ingeniería metabólica y biología sintética para el diseño racional de cepas microbianas. 100% open-source y auditable.

[![Python](https://img.shields.io/badge/Python-3.10%2B-blue.svg)](https://www.python.org/)
[![License: PolyForm Noncommercial](https://img.shields.io/badge/License-PolyForm%20Noncommercial%201.0.0-blue.svg)](https://polyformproject.org/licenses/noncommercial/1.0.0)
[![Status](https://img.shields.io/badge/Status-Beta-orange.svg)]()
[![Version](https://img.shields.io/badge/Version-v0.1.0--beta-brightgreen.svg)]()
[![COBRApy](https://img.shields.io/badge/COBRApy-%E2%89%A50.29-blue.svg)](https://opencobra.github.io/cobrapy/)
[![eQuilibrator](https://img.shields.io/badge/eQuilibrator-%E2%89%A50.4-purple.svg)](https://equilibrator.weizmann.ac.il/)
[![Tests](https://github.com/diegosancriss16-tech/MetaboEngineer-Pro/actions/workflows/tests.yml/badge.svg)](https://github.com/diegosancriss16-tech/MetaboEngineer-Pro/actions/workflows/tests.yml)

**MetaboEngineer Pro** integra algoritmos genéticos, acoplamiento de crecimiento (OptKnock), validación termodinámica (eQuilibrator) y extracción de datos genómicos (UniProt/Rhea) en una sola interfaz. Implementa un **Filtro de Esencialidad validado gen-a-gen** y simulación de ambientes aerobios/anaerobios.

---

## 📋 Tabla de Contenidos

- [Características](#-características)
- [Metodología Híbrida y Acoplamiento (Knife-Edge)](#-metodología-híbrida-y-acoplamiento-knife-edge)
- [Validación Biológica (Keio / Monk)](#-validación-biológica-keio--monk)
- [Resultados del GA (El problema de la estocasticidad)](#-resultados-del-ga-el-problema-de-la-estocasticidad)
- [Instalación](#-instalación)
- [Uso rápido](#-uso-rápido)
- [Limitaciones conocidas](#-limitaciones-conocidas)
- [Roadmap](#-roadmap)
- [Citación y Predecesores](#-citación-y-predecesores)
- [Licencia](#-licencia)

---

## ✨ Características

| Módulo | Descripción |
|---|---|
| **Búsqueda heurística (GA)** | Algoritmo genético para explorar combinaciones de knockouts. |
| **OptKnock exacto (k≤2)** | Enumeración exhaustiva para obtener el óptimo demostrable. |
| **Filtro de Esencialidad** | Validado gen-a-gen contra datos experimentales (Keio). Evita soluciones letales. |
| **Ambientes Controlados** | Simulación de condiciones **aerobias** o **anaerobias**. |
| **Termodinámica** | ΔG' de reacciones vía eQuilibrator. |
| **Reporte científico** | Generación automática de reporte `.md` con UniProt + Rhea + cinética + termo. |
| **Mapa visual** | Red metabólica interactiva (PyVis/HTML). |
| **Universal** | E. coli, S. cerevisiae, B. subtilis, M. tuberculosis y cualquier modelo SBML. |
| **Caché inteligente** | Los modelos se guardan, indexan y verifican (SHA256) automáticamente. |

---

## 🔬 Metodología Híbrida y Acoplamiento (Knife-Edge)

La herramienta utiliza enumeración exhaustiva para \(k \le 2\) y un GA para \(k \ge 3\). Para el diseño de cepas acopladas, implementamos un fitness basado en la **producción mínima garantizada a 100% del crecimiento máximo (FVA-mínimo)**.

**⚠️ Caveat Biológico (El problema del *Knife-Edge*):**
El acoplamiento crecimiento-producción es extremadamente frágil. Evaluamos el knockout canónico `atpB` y descubrimos que su garantía colapsa si la célula no crece exactamente en su óptimo (1.5562 @ 100% vs 0.0337 @ 90%). Esto significa que la métrica mide el "peor caso teórico en el óptimo", no el "peor caso en un biorreactor real". Este es un caveat conocido de la literatura y la razón por la que RobustKnock existe.

---

## 📊 Validación Biológica (Keio / Monk)

Para validar el motor de esencialidad, comparamos `single_gene_deletion` contra el dataset de Keio (1,323 genes en común). Mostramos el trade-off completo sin *cherry-picking*:

| Condición experimental | Precisión | Recall | F1-Score |
|---|---|---|---|
| **M9 glucosa (mínimo)** | 43.48% | **72.73%** | 54.42% |
| **LB rico (igualado a Keio)** | **82.89%** | 57.27% | **67.74%** |

**Interpretación:** El modelo muestra el trade-off clásico precision/recall. En LB (condiciones reales de Keio) la precisión sube al 82.89% (no hay falsos positivos de biosíntesis). En M9 el recall sube al 72.73%, consistente con la literatura para modelos FBA no curados (Monk et al. 2017 reporta ~85% tras curación manual). Este es un resultado **real, auditable y consistente con la literatura**.

*(Nota forense: Durante el desarrollo, una aparente discrepancia contra COBRApy resultó ser un bug en el test de validación (iteración posicional del DataFrame). Tras corregirlo, ambos métodos coinciden en el 100% de los genes. Este ciclo de aseguramiento de calidad está documentado en `tests/`).*

---

## 🧬 Resultados del GA (El problema de la estocasticidad)

Evaluamos el GA con fitness FVA-min@100% en `iML1515` (anaerobio, succinato, \(k=3\)) con 3 semillas:

| Semilla | Genes KO | Fitness (FVA-min @ 100%) |
|---|---|---|
| 42 | nuoI, fau, yejB | 0.0526 (Óptimo local) |
| **100** | iscA, arnD, **atpA** | **1.5561 (Éxito)** |
| 2024 | tktB, nagK, fepD | 0.0526 (Óptimo local) |

**Análisis honesto:**
*   **Hit-rate:** 1/3 (33%). El GA encontró la solución canónica (`atpA`, subunidad de la ATP sintasa) en 1 de 3 intentos. Las otras 2 semillas convergieron a knockouts silenciosos (fitness = WT).
*   **Pasajeros silenciosos:** Validamos la solución ganadora por partes. El gen `atpA` solo otorga un fitness de **1.5562**, mientras que los genes `iscA` y `arnD` juntos solo dan **0.0526** (igual al WT). El GA redescubrió el óptimo de \(k=1\) (`atpB`/`atpA`), pero no demostró que \(k=3\) encuentre algo mejor que \(k=1\). Los otros dos genes son pasajeros silenciosos.
*   **Limitación:** El GA es estocástico. Se recomienda ejecución multi-semilla (Best-of-N) o siembra de población con targets conocidos (mejora prioritaria para v0.2).

---

## 🚀 Instalación

### Requisitos
- Python **3.10+**
- 4 GB de RAM mínimo (8 GB recomendado para modelos grandes)
- Conexión a internet para descargar modelos y consultar APIs (UniProt, Rhea)

### Instalación rápida

```bash
# 1. Clona el repositorio
git clone https://github.com/diegosancriss16-tech/MetaboEngineer-Pro.git
cd MetaboEngineer-Pro

# 2. Crea y activa entorno virtual
python -m venv venv
# Windows:
venv\Scripts\activate
# Linux/Mac:
source venv/bin/activate

# 3. Instala dependencias
pip install -r requirements.txt

# 4. Ejecuta
python metabo_engineer.py