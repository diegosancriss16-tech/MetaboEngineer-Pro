# 🧬 MetaboEngineer Pro

> Suite integral de ingeniería metabólica y biología sintética para el diseño racional de cepas microbianas.

[![Python](https://img.shields.io/badge/Python-3.10%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Status](https://img.shields.io/badge/Status-Stable-green.svg)]()
[![Version](https://img.shields.io/badge/Version-5.2-brightgreen.svg)]()
[![COBRApy](https://img.shields.io/badge/COBRApy-%E2%89%A50.29-blue.svg)](https://opencobra.github.io/cobrapy/)
[![eQuilibrator](https://img.shields.io/badge/eQuilibrator-%E2%89%A50.4-purple.svg)](https://equilibrator.weizmann.ac.il/)

**MetaboEngineer Pro** integra algoritmos genéticos, acoplamiento de crecimiento (OptKnock), validación termodinámica (eQuilibrator), extracción de datos genómicos (UniProt/Rhea) y visualización interactiva de redes metabólicas en una sola interfaz de línea de comandos.

---

## 📋 Tabla de Contenidos

- [Características](#-características)
- [Instalación](#-instalación)
- [Uso rápido](#-uso-rápido)
- [Ejemplos](#-ejemplos)
- [Arquitectura](#-arquitectura)
- [Estructura del proyecto](#-estructura-del-proyecto)
- [Formatos de salida](#-formatos-de-salida)
- [Repositorios soportados](#-repositorios-soportados)
- [Limitaciones conocidas](#-limitaciones-conocidas)
- [Roadmap](#-roadmap)
- [Citación](#-citación)
- [Licencia](#-licencia)

---

## ✨ Características

| Módulo | Descripción |
|---|---|
| **FBA Básico** | Análisis de flujos (Flux Balance Analysis) con COBRApy |
| **Diseño por IA** | Algoritmo genético para encontrar knockouts que maximicen la producción |
| **OptKnock Híbrido** | Acoplamiento crecimiento-producción (cepas estables en biorreactor) |
| **Termodinámica** | ΔG' de reacciones vía eQuilibrator (viabilidad energética) |
| **Reporte científico** | Generación automática de reporte `.md` con UniProt + Rhea + termo |
| **Mapa visual** | Red metabólica interactiva (PyVis/HTML) |
| **Multi-organismo** | E. coli, S. cerevisiae, B. subtilis, M. tuberculosis y cualquier modelo SBML |
| **Caché inteligente** | Los modelos se guardan, indexan y verifican (SHA256) automáticamente |
| **Biología sintética** | Inyección de rutas heterólogas (ej. cafeína) |
| **Interfaz web** | Streamlit opcional para uso sin terminal |

---

## 🚀 Instalación

### Requisitos

- Python **3.10+**
- 4 GB de RAM mínimo (8 GB recomendado para modelos grandes)
- Conexión a internet para descargar modelos y consultar APIs (UniProt, Rhea)

### Instalación rápida

```bash
# 1. Clona el repositorio
git clone https://github.com/<Diegosacris>/MetaboEngineer-Pro.git
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
│ MetaboEngineer Pro V5.2                      │
│ Suite Universal de Biologia de Sistemas      │
╰──────────────────────────────────────────────╯

Selecciona el organismo base:
  1. E. coli (core)
  2. E. coli (completo, MG1655)
  3. Saccharomyces cerevisiae (levadura)
  4. Bacillus subtilis
  5. Otro (ingresar ID de BiGG / BioModels / local)

Modelos ya en cache local (listos sin internet):
  [c1] iEK1008 — 1226 reacciones, 890 genes, 2025-01-15
  [c0] Seleccion normal (descargar nuevo)
```

### CLI experto (línea de comandos, sin menús)

```bash
python cli_experto.py \
    --modelo iML1515 \
    --quimico EX_succ_e \
    --metodo ia \
    --generaciones 20 \
    --seed 42 \
    --output reportes/
```

---

## 🧪 Ejemplos

### Ejemplo 1: Diseñar una cepa productora de succinato en E. coli

```bash
python metabo_engineer.py
# → Opción 2 (E. coli completo)
# → Opción 2 (Diseño por IA)
# → Químico: EX_succ_e
# → 10 generaciones, 20 individuos, 4 knockouts
```

**Salida esperada:**

```
[INFO] Diseno de cepas por IA para E. coli (completo, MG1655)...
Quimico objetivo (ej. EX_succ_e, EX_caft_e): EX_succ_e
Generaciones [10]:
Poblacion [20]:
Genes a apagar [4]:
Evolucionando... ━━━━━━━━━━━━━━━━━━━━━━━━ 100%

[OK] Diseno completado.
Fitness: 12.3456

[INFO] Traduciendo nombres de genes con UniProt...
╭────────────── Genes a apagar (Knockouts) ───────────────╮
│ ID (BiGG) │ Nombre Real  │ Funcion Proteica            │
│ b0351     │ mhpF         │ Acetaldehyde dehydrogenase  │
│ b1241     │ adhE         │ Aldehyde-alcohol DH         │
│ b2296     │ ackA         │ Acetate kinase              │
│ b3115     │ tdcD         │ Propionate kinase           │
╰──────────────────────────────────────────────────────────╯

[OK] Reporte guardado: reportes/REPORTE_iML1515_EX_succ_e_2025-01-15_11-30.md
```

### Ejemplo 2: Acoplar crecimiento a producción (OptKnock)

```bash
python metabo_engineer.py
# → Opción 2 (E. coli completo)
# → Opción 3 (OptKnock)
# → Químico: EX_succ_e
# → Máximo 3 knockouts
```

**Interpretación:** un ratio alto (> 1) significa que la cepa produce succinato de forma natural al crecer.

### Ejemplo 3: Análisis integral con termodinámica + mapa + reporte

```bash
python metabo_engineer.py
# → Opción 2 (E. coli completo)
# → Opción 4 (Análisis Integral)
```

Genera automáticamente:
- Reporte `.md` con UniProt + Rhea + cinética + termo
- Mapa `.html` interactivo
- Todos los análisis de fases 1-4 en consola

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
├── cli_experto.py                # CLI no interactivo (para HPC)
├── app_web.py                    # Interfaz Streamlit
├── modulo_termodinamica_v2.py    # eQuilibrator wrapper
├── visualizador_universal.py     # Red metabólica PyVis
├── reporte_investigador.py       # Reporte .md con UniProt/Rhea
│
├── modelos/                      # Cache local de modelos SBML
│   ├── _index.json               # Índice con SHA256, fecha, metadata
│   ├── iML1515.xml
│   └── iEK1008.xml
│
├── reportes/                     # Salidas generadas
│   ├── REPORTE_*.md
│   ├── mapa_*.html
│   └── historial.csv             # Registro acumulativo de corridas
│
├── tests/                        # Tests con pytest
│   └── test_core.py
│
├── .github/
│   └── workflows/
│       └── tests.yml             # CI: corre tests en cada push
│
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

**Fecha:** 2025-01-15 11:30
**Organismo Base:** `iML1515`
**Reacciones totales del modelo:** 2712

## 1. Resumen Ejecutivo

- **Rendimiento ajustado (Fitness):** `12.3456 mmol/gDW/h`
- **Genes knock-out:** `b0351, b1241, b2296, b3115`

## 2. Validacion Genomica (UniProt & Rhea)

| Gen (ID) | Nombre | Funcion Proteica | Rhea ID |
| :--- | :--- | :--- | :--- |
| `b0351` | **mhpF** | Acetaldehyde dehydrogenase | RHEA:26869 |
| `b1241` | **adhE** | Aldehyde-alcohol dehydrogenase | RHEA:26333 |
...
```

### JSON estructurado (`.json`)

Para análisis en Python/R:

```json
{
  "modelo": "iML1515",
  "quimico": "EX_succ_e",
  "fitness": 12.3456,
  "genes_knockout": ["b0351", "b1241", "b2296", "b3115"],
  "parametros": {"generaciones": 10, "poblacion": 20, "mutaciones": 4, "seed": 42},
  "fecha": "2025-01-15T11:30:00"
}
```

### CSV acumulativo (`historial.csv`)

Todas las corridas en un solo archivo:

```csv
fecha,modelo,quimico,metodo,fitness,genes,seed
2025-01-15T11:30,iML1515,EX_succ_e,ia,12.3456,"b0351;b1241;b2296;b3115",42
2025-01-15T12:00,iML1515,EX_succ_e,optknock,13.2100,"b0351;b2296",42
```

---

## 🌐 Repositorios soportados

| Repositorio | Método de acceso | Estado |
|---|---|---|
| **BiGG Models** | `load_model()` + descarga directa | ✅ Automático |
| **BioModels** | `load_model()` | ✅ Automático |
| **KBase** | Manual (descarga SBML) | ⚠️ Manual |
| **ModelSEED** | Manual (descarga SBML) | ⚠️ Manual |
| **MetaNetX** | Manual (descarga SBML) | ⚠️ Manual |
| **Local** | Carpeta `modelos/` | ✅ Automático |

**Instrucciones manuales:**

Si BiGG está caído, descarga el XML desde [bigg.ucsd.edu](http://bigg.ucsd.edu) y guárdalo en `modelos/<id>.xml`. La próxima ejecución lo cargará desde disco sin tocar internet.

---

## ⚠️ Limitaciones conocidas

- **BiGG puede fallar:** el servidor tiene caídas intermitentes. Usa la descarga manual si es necesario.
- **OptKnock heurístico:** usamos un algoritmo genético, no bilevel optimization exacto (MILP). Resultado aproximado.
- **Cinética ausente:** los parámetros `kcat/Km` no están incluidos por defecto (requieren API autenticada de BRENDA).
- **Modelos grandes:** `iML1515` y similares pueden tardar 3-10 min por corrida con IA.
- **eQuilibrator:** requiere conexión a internet para descargar compuestos.

---

## 🗺️ Roadmap

- [x] FBA + Algoritmo genético
- [x] OptKnock híbrido
- [x] Termodinámica universal
- [x] Reportes científicos
- [x] Caché + indexación de modelos
- [x] Interfaz web (Streamlit)
- [ ] CLI experto con argparse
- [ ] Tests con pytest + CI
- [ ] Docker
- [ ] Caché UniProt offline
- [ ] Exportación JSON/CSV
- [ ] Comparador de cepas
- [ ] Integración con BRENDA API (si se consigue acceso)

---

## 📚 Citación

Si usas MetaboEngineer Pro en una publicación, cita:

```bibtex
@software{metaboengineer2025,
  author  = {San Miguel, Diego},
  title   = {MetaboEngineer Pro: A Suite for Metabolic Engineering and Synthetic Biology},
  year    = {2025},
  version = {5.2},
  url     = {https://github.com/<Diegosacris>/MetaboEngineer-Pro}
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