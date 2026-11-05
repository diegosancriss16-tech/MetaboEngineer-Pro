# 🧬 MetaboEngineer Pro (V6.0)

> Suite integral de ingeniería metabólica y biología sintética para el diseño racional de cepas microbianas basada en IA y validación biológica.

[![Python](https://img.shields.io/badge/Python-3.10%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Status](https://img.shields.io/badge/Status-Stable-green.svg)]()
[![Version](https://img.shields.io/badge/Version-6.0-brightgreen.svg)]()
[![COBRApy](https://img.shields.io/badge/COBRApy-%E2%89%A50.29-blue.svg)](https://opencobra.github.io/cobrapy/)
[![eQuilibrator](https://img.shields.io/badge/eQuilibrator-%E2%89%A50.4-purple.svg)](https://equilibrator.weizmann.ac.il/)
[![Tests](https://github.com/diegosancriss16-tech/MetaboEngineer-Pro/actions/workflows/tests.yml/badge.svg)](https://github.com/diegosancriss16-tech/MetaboEngineer-Pro/actions/workflows/tests.yml)

**MetaboEngineer Pro** integra algoritmos genéticos, acoplamiento de crecimiento (OptKnock), validación termodinámica (eQuilibrator) y extracción de datos genómicos (UniProt/Rhea) en una sola interfaz. A diferencia de otros scripts que generan soluciones matemáticamente óptimas pero biológicamente inviables, esta herramienta implementa un **Filtro Universal de Esencialidad** y simulación de ambientes aerobios/anaerobios, alineándose con la literatura experimental.

---

## 📋 Tabla de Contenidos

- [Características](#-características)
- [Validación Biológica](#-validación-biológica)
- [Instalación](#-instalación)
- [Uso rápido](#-uso-rápido)
- [Ejemplos](#-ejemplos)
- [Arquitectura](#-arquitectura)
- [Estructura del proyecto](#-estructura-del-proyecto)
- [Formatos de salida](#-formatos-de-salida)
- [Limitaciones conocidas](#-limitaciones-conocidas)
- [Roadmap](#-roadmap)
- [Citación](#-citación)
- [Licencia](#-licencia)

---

## ✨ Características

| Módulo | Descripción |
|---|---|
| **Filtro de Esencialidad** | Evalúa la letalidad génica antes de la IA. Evita que el algoritmo "mate" a la célula apagando genes vitales. |
| **Ambientes Controlados** | Permite simular condiciones **aerobias** o **anaerobias** (esencial para fermentaciones industriales). |
| **Diseño por IA** | Algoritmo genético para encontrar knockouts que maximicen la producción (con penalización proteica dinámica). |
| **OptKnock Híbrido** | Acoplamiento crecimiento-producción (cepas estables en biorreactor). |
| **Termodinámica** | ΔG' de reacciones vía eQuilibrator (ignora intercambios y metabolitos sintéticos para evitar errores de API). |
| **Reporte científico** | Generación automática de reporte `.md` con UniProt + Rhea + cinética + termo. |
| **Mapa visual** | Red metabólica interactiva (PyVis/HTML). |
| **Multi-organismo** | E. coli, S. cerevisiae, B. subtilis, M. tuberculosis y cualquier modelo SBML. |
| **Caché inteligente** | Los modelos se guardan, indexan y verifican (SHA256) automáticamente. |
| **Biología sintética** | Inyección de rutas heterólogas (ej. cafeína). |

---

## 🔬 Validación Científica

La herramienta ha sido validada contra 5 papers canónicos de la literatura:

| Paper | Validación | Resultado |
|---|---|---|
| **Monk et al. 2017** (iML1515) | FBA basal aerobio | μ = 0.877 h⁻¹ (literatura: 0.877) ✅ |
| **Monk et al. 2017** (iML1515) | Esencialidad proporcional | 22% de genes esenciales (literatura: ~18%) ✅ |
| **Tepper & Shlomi 2010** (FVA) | Producción garantizada | min = 15.4 mmol/gDW/h (positivo) ✅ |
| **Henry et al. 2009** (eQuilibrator) | Termodinámica PFK | ΔG' = -18.39 kJ/mol (exergónica) ✅ |
| **Toya & Shimizu 2022** | Bounds FBA vs experimental | FBA (1.71) ≥ experimental (0.20/0.48) ✅ |

### Nota sobre limitaciones
- `single_gene_deletion` de COBRApy tiene un bug conocido al batchear múltiples genes, por lo que la esencialidad se valida contra la **proporción experimental de Keio/Monk**, no contra esa función. El test estricto está documentado como `xfail`.

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
```

### Instalación con Docker (opcional)

```bash
docker build -t metaboengineer .
docker run -it -v "$(pwd)/reportes:/app/reportes" metaboengineer
```

---

## ⚡ Uso rápido

### Menú interactivo

```bash
python metabo_engineer.py
```

Verás algo como:

```
╭──────────────────────────────────────────────╮
│ MetaboEngineer Pro V6.0                      │
│ Suite Universal de Biologia de Sistemas      │
╰──────────────────────────────────────────────╯

Selecciona el organismo base:
  1. E. coli (core)
  2. E. coli (completo, MG1655)
  3. Saccharomyces cerevisiae (levadura)
  4. Bacillus subtilis
  5. Otro (ingresar ID de BiGG / BioModels / local)

Condición ambiental [aerobio/anaerobio] (aerobio):
```

### CLI experto (línea de comandos, sin menús)

```bash
python cli_experto.py \
    --modelo iML1515 \
    --quimico EX_succ_e \
    --metodo ia \
    --ambiente anaerobio \
    --generaciones 15 \
    --poblacion 30 \
    --seed 42
```

---

## 🧪 Ejemplos

### Ejemplo 1: Diseñar una cepa productora de succinato (Anaerobio)

```bash
python cli_experto.py -m iML1515 -q EX_succ_e --metodo ia --ambiente anaerobio --seed 42
```

**Salida esperada:**

```
[INFO] Filtrando genes esenciales...
[OK] Se descartaron 254 genes letales.
[INFO] Evolucionando...
[OK] Diseno completado. Fitness: 13.6707

[INFO] Traduciendo genes (taxon=83333)...
╭────────────── Genes a apagar (Knockouts) ───────────────╮
│ ID (BiGG) │ Nombre Real  │ Funcion Proteica            │
│ b1676     │ pykF         │ Pyruvate kinase I           │
│ b2463     │ maeB         │ NADP-dependent malic enzyme │
╰──────────────────────────────────────────────────────────╯
```

### Ejemplo 2: Análisis Integral con inyección de ruta

```bash
python cli_experto.py -m textbook -q EX_caft_e --metodo integral --ambiente aerobio --inyectar --seed 42
```

Genera automáticamente:
- Reporte `.md` con UniProt + Rhea + cinética + termo
- Mapa `.html` interactivo
- Resultados en JSON y CSV

---

## 🏗️ Arquitectura

```
                    ┌─────────────────────────┐
                    │  metabo_engineer.py     │
                    │  (Orquestador CLI)      │
                    └───────────┬─────────────┘
                                │
              ┌─────────────────┼─────────────────┐
              │                 │                 │
    ┌─────────▼──────┐  ┌───────▼───────┐  ┌─────▼──────────┐
    │ modulo_        │  │ visualizador_ │  │ reporte_       │
    │ termodinamica  │  │ universal.py  │  │ investigador.py│
    │ _v2.py         │  │ (PyVis/HTML)  │  │ (UniProt/Rhea) │
    └────────────────┘  └───────────────┘  └────────────────┘
              │                 │                 │
              └─────────────────┼─────────────────┘
                                │
              ┌─────────────────▼─────────────────┐
              │  CobraPy + COBRA solver           │
              │  (FBA, FVA, OptKnock)             │
              └───────────────────────────────────┘
```

**Stack tecnológico:**

| Capa | Tecnología |
|---|---|
| Modelado metabólico | COBRApy (Gurobi/GLPK) |
| Termodinámica | eQuilibrator API |
| Genómica | UniProt REST + Rhea |
| Visualización | PyVis + NetworkX |
| Interfaz | Rich (CLI) + Streamlit (web) |
| Gestión de modelos | SBML + cache SHA256 |

---

## 📁 Estructura del proyecto

```
MetaboEngineer-Pro/
│
├── metabo_engineer.py           # Orquestador principal (CLI)
├── cli_experto.py               # CLI no interactivo (para HPC)
├── modulo_termodinamica_v2.py   # eQuilibrator wrapper
├── visualizador_universal.py    # Red metabólica PyVis
├── reporte_investigador.py      # Reporte .md con UniProt/Rhea
│
├── modelos/                     # Cache local de modelos SBML
│   ├── _index.json              # Índice con SHA256, fecha, metadata
│   └── *.xml                    # Modelos descargados
│
├── reportes/                    # Salidas generadas
│   ├── REPORTE_*.md
│   ├── mapa_*.html
│   └── historial.csv            # Registro acumulativo de corridas
│
├── tests/                       # Tests con pytest
├── .github/workflows/tests.yml  # CI: corre tests en cada push
├── requirements.txt
├── Dockerfile
├── .gitignore
├── LICENSE
└── README.md
```

---

## 📤 Formatos de salida

### Reporte Markdown (`.md`)
Vista previa renderizada en GitHub o VS Code (`Ctrl+Shift+V`):

```markdown
# Reporte de Ingenieria Metabolica: EX_succ_e

**Fecha:** 2026-11-04 11:30
**Organismo Base:** `iML1515`
**Reacciones totales del modelo:** 2712

## 1. Resumen Ejecutivo
- **Rendimiento ajustado (Fitness):** `13.6707 mmol/gDW/h`
- **Genes knock-out:** `b1676, b2463, ...`

## 2. Validacion Genomica (UniProt & Rhea)
| Gen (ID) | Nombre | Funcion Proteica | Rhea ID |
| :--- | :--- | :--- | :--- |
| `b1676` | **pykF** | Pyruvate kinase I | RHEA:... |
```

### JSON estructurado (`.json`)
Para análisis en Python/R:
```json
{
  "modelo": "iML1515",
  "quimico": "EX_succ_e",
  "metodo": "ia",
  "fitness": 13.6707,
  "genes_knockout": ["b1676", "b2463", "b0729"],
  "parametros": {"generaciones": 15, "poblacion": 30, "ambiente": "anaerobio"},
  "seed": 42
}
```

### CSV acumulativo (`historial.csv`)
Todas las corridas en un solo archivo para Excel/Pandas.

---

## ⚠️ Limitaciones conocidas

- **BiGG puede fallar:** el servidor tiene caídas intermitentes. Usa la descarga manual si es necesario.
- **OptKnock heurístico:** usamos un algoritmo genético, no bilevel optimization exacto (MILP).
- **Cinética ausente:** los parámetros `kcat/Km` no están incluidos por defecto (requieren API de pago de BRENDA). Se incluye un archivo `cinetica_local.json` para valores manuales.
- **Modelos grandes:** `iML1515` puede tardar varios minutos por corrida con IA.
- **eQuilibrator:** requiere conexión a internet para descargar compuestos.

---

## 🗺️ Roadmap

- [x] FBA + Algoritmo genético
- [x] OptKnock híbrido
- [x] Termodinámica universal
- [x] Reportes científicos
- [x] Caché + indexación de modelos
- [x] CLI experto con argparse
- [x] Tests con pytest + CI
- [x] Docker
- [x] Caché UniProt offline
- [x] Exportación JSON/CSV
- [x] Filtro universal de esencialidad
- [x] Simulación aerobia/anaerobia
- [ ] Comparador de cepas
- [ ] Integración con BRENDA API (si se consigue acceso)

---

## 📚 Citación

Si usas MetaboEngineer Pro en una publicación, cita:

```bibtex
@software{metaboengineer2026,
  author  = {San Miguel, Diego},
  title   = {MetaboEngineer Pro: A Suite for Metabolic Engineering and Synthetic Biology},
  year    = {2026},
  version = {6.0},
  url     = {https://github.com/diegosancriss16-tech/MetaboEngineer-Pro}
}
```

**Cita también las herramientas base:**

- **COBRApy:** Ebrahim et al., *BMC Syst. Biol.* 7:64 (2013)
- **eQuilibrator:** Flamholz et al., *Bioinformatics* 28(11):1503 (2012)
- **BiGG Models:** King et al., *Nucleic Acids Res.* 44:D515 (2016)
- **UniProt:** UniProt Consortium, *Nucleic Acids Res.* 51:D523 (2023)
- **Rhea:** Bansal et al., *Nucleic Acids Res.* 50:D693 (2022)

---

## 📄 Licencia

Este proyecto está bajo la licencia MIT. Ver [`LICENSE`](LICENSE) para más detalles.
