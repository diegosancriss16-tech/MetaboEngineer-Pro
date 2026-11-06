# 🧬 MetaboEngineer Pro (v0.1.0-beta)

> Suite integral de ingeniería metabólica y biología sintética para el diseño racional de cepas microbianas. 100% open-source y auditable.

[![Python](https://img.shields.io/badge/Python-3.10%2B-blue.svg)](https://www.python.org/)
[![License: PolyForm Noncommercial](https://img.shields.io/badge/License-PolyForm%20Noncommercial%201.0.0-blue.svg)](https://polyformproject.org/licenses/noncommercial/1.0.0/)
[![Status](https://img.shields.io/badge/Status-Beta-orange.svg)]()
[![Version](https://img.shields.io/badge/Version-v0.1.0--beta-brightgreen.svg)]()
[![COBRApy](https://img.shields.io/badge/COBRApy-%E2%89%A50.29-blue.svg)](https://opencobra.github.io/cobrapy/)
[![eQuilibrator](https://img.shields.io/badge/eQuilibrator-%E2%89%A50.4-purple.svg)](https://equilibrator.weizmann.ac.il/)
[![Tests](https://github.com/diegosancriss16-tech/MetaboEngineer-Pro/actions/workflows/tests.yml/badge.svg)](https://github.com/diegosancriss16-tech/MetaboEngineer-Pro/actions/workflows/tests.yml)

**MetaboEngineer Pro** integra algoritmos genéticos, acoplamiento de crecimiento (OptKnock), validación termodinámica (eQuilibrator) y extracción de datos genómicos (UniProt/Rhea) en una sola interfaz. Implementa un **Filtro de Esencialidad validado gen-a-gen** y simulación de ambientes aerobios/anaerobios.

---

## 📋 Tabla de Contenidos

- [Características](#-características)
- [Metodología Híbrida y Acoplamiento](#-metodología-híbrida-y-acoplamiento)
- [Validación Biológica (Keio)](#-validación-biológica-keio)
- [Resultados del GA (Epistasis y Estocasticidad)](#-resultados-del-ga-epistasis-y-estocasticidad)
- [Validación Experimental (Cepa KJ060 y el Valle de Aptitud)](#-validación-experimental-cepa-kj060-y-el-valle-de-aptitud)
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
| **Universal** | *E. coli*, *S. cerevisiae*, *B. subtilis*, *M. tuberculosis* y cualquier modelo SBML. |
| **Caché inteligente** | Los modelos se guardan, indexan y verifican (SHA256) automáticamente. |

---

## 🔬 Metodología Híbrida y Acoplamiento

La herramienta utiliza enumeración exhaustiva para k ≤ 2 y un GA para k ≥ 3. Para el diseño de cepas acopladas, implementamos un fitness basado en la **producción mínima garantizada a 100% del crecimiento máximo (FVA-mínimo)**.

**⚠️ Caveat Biológico (El problema del *Knife-Edge*):**
El acoplamiento basado en optimización energética (ATP sintasa) es extremadamente frágil: la garantía de `atpB` colapsa si la célula crece por debajo de su óptimo (1.5562 @ 100% → 0.0337 @ 90%, retención 2.2%). En contraste, el diseño experimental de Jantama (bloqueo de sumideros redox) retiene el 81.9% de su acoplamiento a 90% del crecimiento — **37× más robusto** que la estrategia computacional canónica. Ver sección [Validación Experimental](#-validación-experimental-cepa-kj060-y-el-valle-de-aptitud).

---

## 📊 Validación Biológica (Keio / Monk)

Para validar el motor de esencialidad, comparamos `single_gene_deletion` contra el dataset de Keio (1,323 genes en común). Mostramos el trade-off completo sin *cherry-picking*:

| Condición experimental | Precisión | Recall | F1-Score |
|---|---|---|---|
| **M9 glucosa (mínimo)** | 43.48% | **72.73%** | 54.42% |
| **LB rico (igualado a Keio)** | **82.89%** | 57.27% | **67.74%** |

**Interpretación:** El modelo muestra el trade-off clásico precision/recall. En LB (condiciones reales de Keio) la precisión sube al 82.89% (no hay falsos positivos de biosíntesis). En M9 el recall sube al 72.73%, consistente con la literatura para modelos FBA no curados (Monk et al. 2017 reporta ~85% tras curación manual).

*(Nota: Durante el desarrollo, una aparente discrepancia contra COBRApy resultó ser un bug en el test de validación. Tras corregirlo, ambos métodos coinciden en el 100% de los genes. Documentado en `tests/`.)*

---

## 🧬 Resultados del GA (Epistasis y Estocasticidad)

### k=3 (experimento inicial, 3 semillas)

| Semilla | Genes KO | Fitness |
|---|---|---|
| 42 | nuoI, fau, yejB | 0.0526 (WT) |
| **100** | iscA, arnD, **atpA** | **1.5561** |
| 2024 | tktB, nagK, fepD | 0.0526 (WT) |

Hit-rate: **1/3 (33%)**.

### k=5 (presupuesto ampliado, 2 semillas)

| Semilla | Genes KO | Fitness |
|---|---|---|
| 42 | atpD + 4 pasajeros | 1.5562 |
| **100** | **atpF + nrdD** + 3 pasajeros | **1.5612** |

Hit-rate: **2/2 (100%)** con presupuesto real (pop=50, gens=20).

**Análisis por Leave-One-Out del diseño ganador (semilla 100):**
- `atpF` (ATP sintasa, subunidad b): núcleo del diseño. Sin él, fitness = 0.0526 (WT).
- **`nrdD` (ribonucleótido reductasa anaerobia)**: contribuidor epistático real. Sin él, fitness = 1.5562 (= atpF solo). Aporta +0.0050 al diseño.
- `gatB`, `crr`, `mglA`: pasajeros silenciosos (sin efecto individual ni combinado).

**Conclusiones:**
- El GA puede superar el óptimo de k=1 mediante epistasis intergénera (ATP sintasa + metabolismo de nucleótidos), aunque el margen es pequeño (+0.32%).
- El GA converge consistentemente a la ATP sintasa como mecanismo de acoplamiento, redescubriendo el diseño canónico de la literatura computacional.
- La estocasticidad es real: en k=3 (presupuesto pequeño) el hit-rate baja al 33%.

---

## 🔬 Validación Experimental (Cepa KJ060 y el Valle de Aptitud)

### Hallazgo 1: La cepa experimental SÍ existe en el modelo — pero requiere k=13

La cepa KJ060 (Jantama et al. 2008; 5 KOs: ldhA, adhE, ackA, focA, pflB; 1.2–1.6 mol/mol succinato experimental) **no es reproducible en iML1515 como deleciones monogénicas**: 4 de sus 5 genes tienen isoenzimas en el modelo (adhE→adhP/mhpF, ackA→tdcD/purT, focA→focB, pflB→yfiD), y solo LDH_D se bloquea efectivamente.

Al extender los knockouts a las **familias completas de isoenzimas (13 KOs)**, el diseño experimental alcanza:

| Diseño | k | FVA-min @100% | FVA-min @90% | Retención |
|---|---|---|---|---|
| WT | 0 | 0.0526 | — | — |
| atpB (óptimo in silico k=1) | 1 | 1.5562 | 0.0337 | 2.2% ❌ |
| GA (mejor semilla, k=5) | 5 | 1.5612 | — | — |
| **KJ060 extendido (experimental)** | **13** | **7.7988** | **6.3910** | **81.9% ✅** |

**Interpretación:** La estrategia experimental (bloquear sumideros redox) domina a la computacional (bloquear ATP sintasa) en ambos ejes: **5× mayor acoplamiento garantizado y 37× más robusta a la relajación del crecimiento**. El mutante de 13 KOs retiene el 42.2% del crecimiento WT (μ=0.0665 h⁻¹), consistente con las cepas de biorreactor.

### Hallazgo 2: El "Valle de Aptitud" (Fitness Valley)

El camino hacia el óptimo global no es monótono:

| Acumulado | n_KO | Fitness | Crecimiento (μ) |
|---|---|---|---|
| WT | 0 | 0.0526 | 0.1575 |
| +lactato (ldhA) | 1 | 0.0526 | 0.1575 |
| +etanol (adhE+adhP+mhpF) | 4 | 0.0394 | 0.1182 |
| +acetato (ackA+tdcD+purT) | 7 | 0.0305 | 0.0915 |
| +formato (focA+focB) | 9 | **7.7872** | 0.0665 |
| +PFL (pflB+yfiD+alt) | 13 | **7.7988** | 0.0665 |

Los primeros 7 knockouts **empeoran el fitness respecto al WT**. El salto de 255× ocurre únicamente al bloquear el transporte de formato — pero solo funciona si los 7 anteriores ya están fuera. Esto es epistasis extrema: cada knockout individualmente inútil (o deletéreo) se vuelve necesario en combinación.

**Implicación para el diseño asistido:** todos los intermedios del valle son viables (μ ≥ 0.0665 > umbral 0.05), pero un GA con selección greedy purga los intermedios peores-que-WT. Cruzar el valle requiere estrategias de *seeding* o fitness multi-objetivo.

---

## 🚀 Instalación

### Requisitos
- Python **3.10+**
- 4 GB de RAM mínimo (8 GB recomendado para modelos grandes)
- Conexión a internet para descargar modelos y consultar APIs (UniProt, Rhea)

### Instalación rápida

1. Clona el repositorio
git clone https://github.com/diegosancriss16-tech/MetaboEngineer-Pro.git
cd MetaboEngineer-Pro

2. Crea y activa entorno virtual
python -m venv venv

Windows:
venv\Scripts\activate

Linux/Mac:
source venv/bin/activate

3. Instala dependencias
pip install -r requirements.txt

4. Ejecuta
python metabo_engineer.py


---

## ⚡ Uso rápido

### Menú interactivo
bash
python metabo_engineer.py


### CLI experto (línea de comandos)
bash
python cli_experto.py
–modelo iML1515
–quimico EX_succ_e
–metodo ia
–ambiente anaerobio
–generaciones 15
–poblacion 30
–seed 42


---

## ⚠️ Limitaciones conocidas

- **BiGG puede fallar:** el servidor tiene caídas intermitentes. Usa la descarga manual si es necesario.
- **OptKnock heurístico:** usamos un GA, no bilevel optimization exacto (MILP). Hit-rate estocástico: 33% en k=3 con presupuesto pequeño, 100% en k=5 con presupuesto real.
- **Atrapamiento en el valle de aptitud:** el GA converge consistentemente al pico de ATP sintasa (~1.56) y no cruza el valle epistático hacia el óptimo global (7.79). El problema es de exploración, no de viabilidad.
- **Cinética ausente:** los parámetros `kcat/Km` no están incluidos por defecto (requieren API de pago de BRENDA).
- **eQuilibrator:** requiere conexión a internet para descargar compuestos.
- **Margen de viabilidad fino en el óptimo:** el diseño de 13 KOs tiene μ=0.0665, solo 1.33× sobre el umbral de viabilidad de la herramienta (0.05).

---

## 🗺️ Roadmap (v0.2.0 y más allá)

- [x] FBA + Algoritmo genético
- [x] OptKnock híbrido
- [x] Termodinámica universal
- [x] Reportes científicos
- [x] Caché + indexación de modelos
- [x] CLI experto con argparse
- [x] Tests con pytest + CI
- [x] Docker
- [x] Filtro universal de esencialidad
- [x] Simulación aerobia/anaerobia
- [x] Validación gen-a-gen contra Keio
- [x] Benchmark de solvers (GLPK, SCIP, HiGHS)
- [x] Validación experimental KJ060 (familias de isozimas, k=13)
- [x] Caracterización del valle de aptitud (5 intermedios, todos viables)
- [x] Leave-one-out del diseño epistático (atpF + nrdD)
- [x] Curva de robustez (retención 81.9% vs 2.2%)
- [ ] **GA con k=13:** ¿Puede el GA cruzar el valle de aptitud y redescubrir el diseño experimental? Criterios pre-registrados: (a) éxito total: fitness ≥ 7.0; (b) parcial: fitness > 1.5612; (c) fracaso informativo: fitness ≈ 1.5562.
- [ ] **Siembra de población:** Incluir familias completas de KJ060 en la población inicial para facilitar el cruce del valle.
- [ ] **Leave-family-out del diseño minimal:** ¿Cuál es el subconjunto mínimo de las 13 que acopla? La familia PFL aporta solo +0.0116.
- [ ] **Restricciones enzimáticas:** Integración con ECMpy (licencia MIT).
- [ ] **Publicación:** Preprint en bioRxiv documentando el benchmark, la validación experimental y el valle de aptitud.

---

## 📚 Citación y Predecesores

## Citación

Si usas este software, cita los algoritmos base en los que se fundamenta:

- **MetaboEngineer Pro** (2025).
- **OptGene**: Patil, K. R., Rocha, I., Förster, J., & Nielsen, J. (2005). Evolutionary programming as a platform for in silico metabolic engineering. *BMC Bioinformatics*, 6(1), 308.
- **OptKnock**: Burgard, A. P., Pharkya, P., & Maranas, C. D. (2003). Optknock: a bilevel programming framework for identifying gene knockout strategies for microbial strain optimization. *Biotechnology and Bioengineering*, 84(6), 647-657.
- **COBRApy**: Ebrahim, A., Lerman, J. A., Palsson, B. O., & Hyduke, D. R. (2013). COBRApy: COnstraints-Based Reconstruction and Analysis for Python. *BMC Systems Biology*, 7(1), 74.
- **Keio Collection**: Baba, T., et al. (2006). Construction of Escherichia coli K-12 in-frame, single-gene knockout mutants: the Keio collection. *Molecular Systems Biology*, 2(1), 2006-0008.

---

## 📄 Licencia

Este proyecto está licenciado bajo la **PolyForm Noncommercial License 1.0.0**.

> **Required Notice:** Copyright Diego San Cristobal (https://github.com/diegosancriss16-tech/MetaboEngineer-Pro)

### ✅ Lo que puedes hacer GRATIS
- Usar la herramienta para investigación académica, docencia o uso personal.
- Modificar el código, hacer fork y crear nuevas versiones.
- Compartir el código con otras personas.
- Usarla en universidades, instituciones públicas u organizaciones sin fines de lucro.

### ❌ Lo que NO puedes hacer sin permiso
- Vender la herramienta o un servicio basado en ella.
- Usarla dentro de una empresa con fines comerciales.
- Ofrecer servicios de consultoría pagados que la utilicen.

### 💼 ¿Necesitas uso comercial?
Si tu empresa u organización quiere usar MetaboEngineer Pro con fines comerciales, contáctame para acordar una **licencia comercial**:
📧 diegosancriss16@gmail.com

Se ofrece **licenciamiento dual**: el código sigue siendo libre para uso no comercial, y se puede negociar una licencia comercial por separado.

### ⚠️ Nota sobre "Open Source"
Esta licencia es **Source-Available** (código disponible), no OSI-Open-Source. El código es público, auditable y modificable, pero el uso comercial requiere permiso. Ver [`LICENSE`](LICENSE) para el texto legal completo.
