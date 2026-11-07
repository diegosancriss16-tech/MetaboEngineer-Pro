> **Required Notice:** Copyright Diego San Miguel Crisantos (https://github.com/diegosancriss16-tech/MetaboEngineer-Pro)

# 🧬 MetaboEngineer Pro (v0.2.0)

> Suite integral de ingeniería metabólica y biología sintética para el diseño racional de cepas microbianas. 100% open-source y auditable.

[![Python](https://img.shields.io/badge/Python-3.10%2B-blue.svg)](https://www.python.org/)
[![License: PolyForm Noncommercial](https://img.shields.io/badge/License-PolyForm%20Noncommercial%201.0.0-blue.svg)](https://polyformproject.org/licenses/noncommercial/1.0.0/)
[![Status](https://img.shields.io/badge/Status-Beta-orange.svg)]()
[![Version](https://img.shields.io/badge/Version-v0.2.0-brightgreen.svg)]()
[![COBRApy](https://img.shields.io/badge/COBRApy-%E2%89%A50.29-blue.svg)](https://opencobra.github.io/cobrapy/)
[![eQuilibrator](https://img.shields.io/badge/eQuilibrator-%E2%89%A50.4-purple.svg)](https://equilibrator.weizmann.ac.il/)
[![Tests](https://github.com/diegosancriss16-tech/MetaboEngineer-Pro/actions/workflows/tests.yml/badge.svg)](https://github.com/diegosancriss16-tech/MetaboEngineer-Pro/actions/workflows/tests.yml)

**MetaboEngineer Pro** integra algoritmos genéticos (con mutación real y reproducibilidad determinista), acoplamiento de crecimiento (OptKnock), validación termodinámica (eQuilibrator) y extracción de datos genómicos (UniProt/Rhea) en una sola interfaz. Implementa un **Filtro de Esencialidad validado gen-a-gen** y simulación de ambientes aerobios/anaerobios con cálculo dinámico de rendimiento biológico.

---

## 📋 Tabla de Contenidos

- [Características](#-características)
- [Metodología Híbrida y Acoplamiento](#-metodología-híbrida-y-acoplamiento)
- [Validación Biológica y Universalidad](#-validación-biológica-y-universalidad)
- [Robustez Adversaria y Seguridad de Datos](#-robustez-adversaria-y-seguridad-de-datos)
- [Resultados del GA (Epistasis, Seeding y Destilación)](#-resultados-del-ga-epistasis-seeding-y-destilación)
- [Validación Experimental (Cepa KJ060)](#-validación-experimental-cepa-kj060)
- [Instalación](#-instalación)
- [Uso rápido y Pipeline de Diseño](#-uso-rápido-y-pipeline-de-diseño)
- [Limitaciones conocidas](#-limitaciones-conocidas)
- [Roadmap](#-roadmap)
- [Citación y Predecesores](#-citación-y-predecesores)
- [Licencia](#-licencia)

---

## ✨ Características

| Módulo | Descripción |
|---|---|
| **Búsqueda heurística (GA)** | Algoritmo genético con cruce y mutación (libre de duplicados) y orden determinista para reproducibilidad. |
| **OptKnock exacto (k≤2)** | Enumeración exhaustiva para obtener el óptimo demostrable. |
| **Filtro de Esencialidad** | Validado gen-a-gen contra datos experimentales (Keio). Implementa **Caché inteligente** (reduce el tiempo de cómputo de 6 horas a segundos). |
| **Ambientes Controlados** | Simulación estricta de condiciones **aerobias** (bloqueo de fermentación) o **anaerobias** (bloqueo de O2). |
| **Pipeline de Diseño** | Script automatizado de 5 etapas (Explorar -> Destilar -> Verificar -> Documentar) listo para laboratorio. |
| **Termodinámica** | ΔG' de reacciones vía eQuilibrator. |
| **Reporte científico** | Generación automática de reporte `.md` con UniProt + Rhea + cálculo dinámico mol/mol. |
| **Universal** | *E. coli*, *B. subtilis*, *P. putida*, *S. cerevisiae* y cualquier modelo SBML. |

---

## 🔬 Metodología Híbrida y Acoplamiento

La herramienta utiliza enumeración exhaustiva para k ≤ 2 y un GA para k ≥ 3. Para el diseño de cepas acopladas, implementamos un fitness basado en la **producción mínima garantizada a 100% del crecimiento máximo (FVA-mínimo)**.

**⚠️ Caveat Biológico (El problema del *Knife-Edge*):**
El acoplamiento basado en optimización energética (ATP sintasa) es extremadamente frágil: la garantía de `atpB` colapsa si la célula crece por debajo de su óptimo (1.5562 @ 100% → 0.0337 @ 90%, retención 2.2%). En contraste, el diseño experimental de Jantama (bloqueo de sumideros redox) retiene el 81.9% de su acoplamiento a 90% del crecimiento — **37× más robusto** que la estrategia computacional clásica.

---

## 📊 Validación Biológica y Universalidad

Para validar el motor de esencialidad, comparamos `single_gene_deletion` contra el dataset de Keio (1,323 genes en común). Mostramos el trade-off completo sin *cherry-picking*:

| Condición experimental | Precisión | Recall | F1-Score |
|---|---|---|---|
| **M9 glucosa (mínimo)** | 43.48% | **72.73%** | 54.42% |
| **LB rico (igualado a Keio)** | **82.89%** | 57.27% | **67.74%** |

### Validación Cruzada Multi-Organismo

Para garantizar que el motor de diseño es **operativo en múltiples organismos**, se ejecutó el pipeline completo en tres especies filogenéticamente distantes:

| Organismo | Modelo | Reacciones | Biomasa Detectada | IDs de Genes | Resultado |
|---|---|---|---|---|---|
| *Escherichia coli* | iML1515 | 2712 | `BIOMASS_Ec_iML1515...` | `bXXXX` | 🏆 Diseño exitoso (atpB) |
| *Bacillus subtilis* | iYO844 | 1250 | `BIOMASS_BS_10` | `BSUXXXX` | ⚠️ GA se quedó corto (honesto) |
| *Pseudomonas putida* | iJN1463 | 2927 | `BIOMASS_KT2440_Core2` | `PP_XXXX` | ⚠️ GA se quedó corto (honesto) |

---

## 🛡️ Robustez Adversaria y Seguridad de Datos

Como herramienta científica, el software debe fallar de forma segura ante entradas inválidas (evitando "silencios mortales" que generen datos falsos). Se implementó una suite de 10 pruebas adversariales en `benchmark/adversarial_test.py`.

| Categoría | Descripción | Resultado |
|---|---|---|
| **S3** | Rechazo limpio (valida antes de calcular) | **9/10** |
| **S2** | Degradación graciosa (atrapa el error y avisa) | **1/10** |
| **S1** | Silencio mortal (devuelve datos falsos) | **0** ✅ |
| **S0** | Crash (error no controlado) | **0** ✅ |

**Certificación:** El software rechaza limpiamente químicos inexistentes, genes con errores tipográficos, ambientes inválidos y modelos corruptos. Cero fallos silenciosos garantizan que los resultados generados son siempre matemáticamente válidos sobre las premisas solicitadas.

---

## 🧬 Resultados del GA (Epistasis, Seeding y Destilación)

### k=5 (presupuesto ampliado, 2 semillas)
| Semilla | Genes KO | Fitness |
|---|---|---|
| 42 | atpD + 4 pasajeros | 1.5562 |
| **100** | **atpF + nrdD** + 3 pasajeros | **1.5612** |

**Análisis por Leave-One-Out del diseño ganador (semilla 100):**
- `atpF` (ATP sintasa): núcleo del diseño.
- **`nrdD` (ribonucleótido reductasa anaerobia)**: contribuidor epistático real (+0.0050).
- Los otros 3 genes son pasajeros silenciosos.

### k=13 (GA puro vs. GA Sembrado + Destilado)
El GA puro (sin intervención) converge a la ATP sintasa, demostrando la barrera del **Valle de Aptitud**. Al inyectar los "escalones" del valle (familias de isoenzimas de Jantama) en la población inicial, el GA **rompió la barrera**.

| Semilla | Método | Genes KO | Fitness |
|---|---|---|---|
| 42 | GA Puro | atpB + ydiO + 11 pasajeros | 1.6444 |
| 100 | GA Puro | atpB + 12 pasajeros | 1.5562 |
| 42 | Seeding | KJ060 (4 genes) + 9 optimizados | 10.8056 |
| **100** | **Seeding + Destilación (LOO)** | **Núcleo de 4 genes (ldhA, adhE, adhP, atpB)** | **10.1242** |

**El Tercer Acto:** El GA sembrado descubrió un diseño que combinaba la estrategia redox con la estrategia energética, alcanzando un fitness **in silico** de **10.1242 mmol/gDW/h** (~1.01 mol/mol). Sin embargo, el verdadero hallazgo fue revelado al aplicar una **destilación iterativa (Leave-One-Out)**: el GA había arrastrado 7 pasajeros silenciosos. El núcleo real de este diseño es de tan solo **4 genes** (`b1380`, `b1241`, `b1478`, `b3738`), los cuales retienen un **79.0%** de su producción al 90% de crecimiento, y lo hacen con un crecimiento robusto de **\(\mu = 0.1071\)** (el doble del umbral de viabilidad).

*Nota de rigor:* El diseño destilado de 4 genes supera el techo computacional *in silico* de KJ060 (7.7988), pero **no supera las mediciones experimentales reales en biorreactor (1.2–1.6 mol/mol, equivalentes a ~12–16 mmol/gDW/h)**.

*(Genes del núcleo destilado: `b1380` (ldhA), `b1241` (adhE), `b1478` (adhP), `b3738` (atpB))*

---

## 🔬 Validación Experimental (Cepa KJ060)

La cepa KJ060 (Jantama et al. 2008; 5 KOs) **no es reproducible como deleciones monogénicas** debido a la redundancia de isoenzimas. Al extender los knockouts a las **familias completas (13 KOs)**, el diseño experimental alcanza:

| Diseño | k | FVA-min @100% | FVA-min @90% | \(\mu\) | Retención |
|---|---|---|---|---|---|
| WT | 0 | 0.0526 | — | 0.1575 | — |
| atpB (óptimo in silico k=1) | 1 | 1.5562 | 0.0337 | — | 2.2% ❌ |
| KJ060 extendido (in silico) | 13 | 7.7988 | 6.3910 | 0.0665 | 81.9% |
| **Núcleo Destilado (in silico)** | **4** | **10.1242** | **7.9986** | **0.1071** | **79.0%** ✅ |

**El Valle de Aptitud:** Los primeros 7 knockouts de KJ060 empeoran el fitness (0.0394 y 0.0305). El salto de 255× ocurre al bloquear el transporte de formato (7.7872). Esto demuestra que el camino hacia el óptimo requiere cruzar pasos transitoriamente deletéreos (epistasis extrema), lo cual requiere estrategias de *seeding*.

---

## 🚀 Instalación

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
⚡ Uso rápido y Pipeline de Diseño
Menú interactivo:

bash
python metabo_engineer.py
CLI experto:

bash
python cli_experto.py --modelo iML1515 --quimico EX_succ_e --ambiente anaerobio --metodo ia --generaciones 20 --poblacion 50 --seed 42
Pipeline de Diseño (Listo para el Laboratorio):
El script benchmark/pipeline_design.py automatiza el flujo de trabajo real de un investigador (5 etapas: Validación -> GA -> Destilación LOO -> Verificación Biológica -> Veredicto).

bash
# Uso estándar (descubre la biomasa automáticamente)
python benchmark/pipeline_design.py --modelo iML1515.xml --target EX_succ_e --ambiente anaerobio --seed 100

# Uso con Seeding manual (para forzar genes conocidos)
python benchmark/pipeline_design.py --modelo iML1515.xml --target EX_succ_e --ambiente anaerobio --seeding b3738,b0002
⚠️ Limitaciones conocidas
Alcance multiorgánico: La herramienta está validada en E. coli y es operativa en otros organismos (levadura, B. subtilis, P. putida). Sin embargo, topologías metabólicas distintas pueden requerir ajustes en los parámetros del GA para superar óptimos locales.

BiGG puede fallar: el servidor tiene caídas intermitentes. Usa la descarga manual si es necesario.

OptKnock heurístico: usamos un GA, no bilevel optimization exacto (MILP).

Dependencia del Seeding: El GA puro no cruza el valle de aptitud hacia el óptimo global (converge a la ATP sintasa). Para encontrar el núcleo de 4 genes se requiere la estrategia de siembra (seeding) + destilación LOO.

Cinética ausente: los parámetros kcat/Km no están incluidos por defecto (requieren API de pago de BRENDA) y han sido removidos de los reportes automáticos para evitar confusiones.

eQuilibrator: requiere conexión a internet para descargar compuestos.

Margen de viabilidad: El núcleo destilado tiene 
μ
=
0.1071
μ=0.1071, lo que le otorga un margen holgado (2.14× sobre el umbral de viabilidad de 0.05).

🗺️ Roadmap (v0.3.0 y más allá)
☑ FBA + Algoritmo genético
☑ OptKnock híbrido
☑ Termodinámica universal
☑ Reportes científicos
☑ Caché + indexación de modelos (reducción de tiempo 6h -> seg)
☑ CLI experto con argparse
☑ Tests con pytest + CI
☑ Docker
☑ Filtro universal de esencialidad
☑ Simulación aerobia/anaerobio
☑ Validación gen-a-gen contra Keio
☑ Benchmark de solvers (GLPK, SCIP, HiGHS)
☑ Validación experimental KJ060 (familias de isozimas, k=13)
☑ Caracterización del valle de aptitud (5 intermedios, todos viables)
☑ Leave-one-out del diseño epistático (atpF + nrdD)
☑ GA con mutación real y reproducibilidad determinista
☑ GA con k=13 (puro): No cruzó el valle (convergió a atpB).
☑ Seeding de población + Destilación LOO iterativa: Rompió el valle y encontró el núcleo de 4 genes (10.12).
☑ Suite de Pruebas Adversariales (0 fallos S0/S1).
☑ Pipeline de Diseño Asistido (pipeline_design.py).
□ Leave-family-out del diseño minimal.
□ Restricciones enzimáticas: Integración con ECMpy (licencia MIT).
□ Publicación: Preprint en bioRxiv documentando el benchmark, la validación experimental y el valle de aptitud.
📚 Citación y Predecesores
bibtex
@software{metaboengineer2025,
  author  = {San Miguel Crisantos, Diego},
  title   = {MetaboEngineer Pro: A Suite for Metabolic Engineering and Synthetic Biology},
  year    = {2025},
  version = {v0.2.0},
  url     = {https://github.com/diegosancriss16-tech/MetaboEngineer-Pro}
}
OptGene: Patil, K. R., et al. (2005). BMC Bioinformatics, 6(1), 308.

OptKnock: Burgard, A. P., et al. (2003). Biotechnology and Bioengineering, 84(6), 647-657.

COBRApy: Ebrahim, A., et al. (2013). BMC Systems Biology, 7(1), 74.

Keio Collection: Baba, T., et al. (2006). Molecular Systems Biology, 2(1), 2006-0008.

iML1515: Monk, J. M., et al. (2017). Nature Biotechnology, 35(10), 904-908.

Jantama et al.: Jantama, K., et al. (2008). Biotechnology and Bioengineering, 99(5), 1140-1153.

eQuilibrator: Flamholz, A., et al. (2012). Nucleic Acids Research, 40(D1), D770-D775.

📄 Licencia
Este proyecto está licenciado bajo la PolyForm Noncommercial License 1.0.0.

Required Notice: Copyright Diego San Miguel Crisantos (https://github.com/diegosancriss16-tech/MetaboEngineer-Pro)

✅ Lo que puedes hacer GRATIS
Usar la herramienta para investigación académica, docencia o uso personal.

Modificar el código, hacer fork y crear nuevas versiones.

Compartir el código con otras personas.

Usarla en universidades, instituciones públicas u organizaciones sin fines de lucro.

❌ Lo que NO puedes hacer sin permiso
Vender la herramienta o un servicio basado en ella.

Usarla dentro de una empresa con fines comerciales.

Ofrecer servicios de consultoría pagados que la utilicen.

💼 ¿Necesitas uso comercial?
Si tu empresa u organización quiere usar MetaboEngineer Pro con fines comerciales, contáctame para acordar una licencia comercial:
📧 diegosancriss16@gmail.com

⚠️ Nota sobre "Open Source": Esta licencia es Source-Available, no OSI-Open-Source. El código es público, auditable y modificable, pero el uso comercial requiere permiso. Ver LICENSE para el texto legal completo.