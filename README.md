> **Required Notice:** Copyright Diego San Miguel Crisantos (https://github.com/diegosancriss16-tech/MetaboEngineer-Pro)

# 🧬 MetaboEngineer Pro (v0.2.8)

> Suite integral de ingeniería metabólica y biología sintética para el diseño racional de cepas microbianas. 100% source-available (uso académico) y auditable.

[![Python](https://img.shields.io/badge/Python-3.10%2B-blue.svg)](https://www.python.org/)
[![License: PolyForm Noncommercial](https://img.shields.io/badge/License-PolyForm%20Noncommercial%201.0.0-blue.svg)](https://polyformproject.org/licenses/noncommercial/1.0.0/)
[![Status](https://img.shields.io/badge/Status-Beta-orange.svg)]()
[![Version](https://img.shields.io/badge/Version-v0.2.8-brightgreen.svg)]()
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
- [Resultados del GA (Epistasis y Óptimos Globales)](#-resultados-del-ga-epistasis-y-óptimos-globales)
- [Comparación con la cepa KJ060 (in silico)](#-comparación-con-la-cepa-kj060-in-silico)
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
| **OptKnock por enumeración (k≤2)** | Enumeración exhaustiva sobre candidatos acotados (genes activos en el WT). No es un MILP exacto ni garantiza el óptimo global. |
| **Filtro de Esencialidad** | Validado gen-a-gen contra datos experimentales (Keio). Implementa **Caché inteligente** (reduce el tiempo de cómputo de 6 horas a segundos). |
| **Ambientes Controlados** | Simulación estricta de condiciones **aerobias** (bloqueo de fermentación) o **anaerobias** (bloqueo de O2). |
| **Pipeline de Diseño** | Script automatizado de 5 etapas (Explorar -> Destilar -> Verificar -> Documentar) listo para laboratorio. |
| **Termodinámica** | ΔG' de reacciones vía eQuilibrator (evaluando el sentido del flujo). |
| **Reporte científico** | Generación automática de reporte `.md` con UniProt + Rhea + cálculo dinámico mol/mol. |
| **Universal** | *E. coli*, *B. subtilis*, *P. putida*, *S. cerevisiae* y cualquier modelo SBML. |

---

## 🔬 Metodología Híbrida y Acoplamiento

La herramienta utiliza enumeración exhaustiva sobre candidatos acotados para k ≤ 2 y un GA para k ≥ 3. Para el diseño de cepas acopladas, implementamos un fitness basado en la **producción mínima garantizada a 100% del crecimiento máximo (FVA-mínimo)**.

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
| *Escherichia coli* | iML1515 | 2712 | `BIOMASS_Ec_iML1515...` | `bXXXX` | ✅ Pipeline operativo. Diseño computacional frágil (atpB); diseño robusto: ptsH + atpD |
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

**Certificación:** El software rechaza limpiamente químicos inexistentes, genes con errores tipográficos, ambientes inválidos y modelos corruptos. En las 10 pruebas adversariales no se detectaron fallos silenciosos (S1) ni crashes (S0); esto no constituye una garantía general, solo cobertura de los casos probados.

---

## 🧬 Resultados del GA (Epistasis y Óptimos Globales)

### k=3 y k=5 (Presupuesto Pequeño)
Con un presupuesto pequeño (población 20, 10 generaciones), el GA converge a la ATP sintasa (`atpB`), demostrando que el paisaje de aptitud tiene un pico local dominante.

### k=4 (GA Puro, Sin Seeding)
El experimento clave: **¿Puede el GA puro encontrar un diseño robusto sin inyectar genes de Jantama?** La respuesta es **SÍ, de forma estocástica** (1 de 5 semillas). Con población 50 y 15 generaciones, el GA puro encontró un diseño de 10.8049 en la semilla 2024, mientras que las otras 4 convergieron al óptimo local de la ATP sintasa o al WT.

| Semilla | Método | Genes KO | Fitness |
|---|---|---|---|
| 100 | GA Puro (k=4) | `b3722, b2239, b2835, b2130` | 0.0526 (WT) |
| 42 | GA Puro (k=4) | `b3721, b4036, b3735, b0505` | 1.5562 (ATP sintasa) |
| 777 | GA Puro (k=4) | `b3566, b1010, b0973, b3738` | 1.5562 (ATP sintasa) |
| 1234 | GA Puro (k=4) | `b3786, b3734, b3785, b1448` | 1.5562 (ATP sintasa) |
| **2024** | **GA Puro (k=4)** | **`b1090, b2415, b3732, b3708`** | **10.8049** 🏆 |
| **100** | **Seeding + Destilación (LOO)** | **`b1380, b1241, b1478, b3738`** | **10.1242** 🏆 |

**El hallazgo inesperado (Núcleo de 2 Genes):** Al aplicar una **destilación iterativa (Leave-One-Out)** al diseño de 10.8049 del GA puro (semilla 2024), se descubrió que el verdadero núcleo es de tan solo **2 genes** que acoplan la entrada de glucosa a la producción de succinato.

**Nota de reproducibilidad entre versiones:** La semilla 2024 con k=4, pop=50 y gens=15 encuentra un **núcleo de 2 genes** en cualquier versión de Python, pero el "socio" de la ATP sintasa varía por pequeñas diferencias numéricas:

| Entorno | Núcleo | Fitness | Socio de `atpD` |
|---|---|---|---|
| **Python 3.12 + COBRApy 0.32.1** | `ptsH` (`b2415`) + `atpD` (`b3732`) | **10.8049** | Sistema PTS (captación de glucosa) |
| **Python 3.11 + COBRApy 0.32.1 (Docker)** | `pgi` (`b4025`) + `atpD` (`b3733`) | **10.3909** | Glucosa-6-fosfato isomerasa (glucólisis) |

Ambos son soluciones de 2 genes que acoplan el crecimiento a la producción de succinato manipulando la entrada de glucosa (PTS o PGI) y la ATP sintasa. La diferencia viene de variaciones numéricas mínimas en `numpy`/`scipy`/`swiglpk` que en la última generación del GA cambian qué individuo gana por un margen mínimo. **Ambos diseños son igual de válidos.** Para reproducibilidad byte a byte, fija las versiones exactas en `requirements.txt` (pendiente para v0.2.9).

**El Tercer Acto (El Núcleo de 4 Genes):** El GA sembrado descubrió un diseño inicial de 11 genes que combinaba la estrategia redox con la estrategia energética (10.1242). Al aplicar una **destilación iterativa (Leave-One-Out)** eliminando pasajeros silenciosos, se descubrió que el núcleo real es de tan solo **4 genes** (`b1380`, `b1241`, `b1478`, `b3738`). Este núcleo retiene un **79.0%** de su producción al 90% de crecimiento, con un crecimiento robusto de **μ = 0.1071** (el doble del umbral de viabilidad).

*Nota de rigor:* Ambos diseños (10.80 y 10.12) superan el techo computacional *in silico* de KJ060 (7.7872), pero **no superan las mediciones experimentales reales en biorreactor (1.2–1.6 mol/mol, equivalentes a ~12–16 mmol/gDW/h)**.

*(Genes del núcleo destilado: `b1380` (ldhA), `b1241` (adhE), `b1478` (adhP), `b3738` (atpB))*

---

## 🔬 Comparación con la cepa KJ060 (in silico)

La cepa KJ060 (Jantama et al. 2008; 5 KOs) **no es reproducible como deleciones monogénicas** debido a la redundancia de isoenzimas. Al extender los knockouts a las **familias completas de isoenzimas (10 genes sin atpE)**, el diseño experimental alcanza:

| Diseño | k | FVA-min @100% | FVA-min @90% | μ | Retención |
|---|---|---|---|---|---|
| WT | 0 | 0.0526 | — | 0.1575 | — |
| atpB (óptimo in silico k=1) | 1 | 1.5562 | 0.0337 | — | 2.2% ❌ |
| **KJ060 familias de isoenzimas (10 genes, sin atpE)** | **10** | **7.7872** | **6.3805** | **0.0665** | **81.9%** |
| **Las mismas + atpE (11 genes)** | **11** | **7.1434** | **5.6754** | **0.0596** | **79.5%** |
| **Núcleo Destilado (in silico)** | **4** | **10.1242** | **7.9986** | **0.1071** | **79.0%** ✅ |
| **Núcleo de 2 genes (ptsH + atpD)** | **2** | **10.8056** | **9.0970** | **0.0695** | **84.2%** ✅ |

*(Nota: El valor original de 7.7988 citado en versiones anteriores de este documento no pudo reproducirse y ha sido retirado. El diseño de Jantama original es ΔldhA ΔadhE ΔackA Δ(focA-pflB), por lo que el atpE no pertenece a KJ060).*

**El Valle de Aptitud y las Soluciones Alternativas:** Los primeros 7 knockouts de KJ060 empeoran el fitness (0.0394 y 0.0305), lo que constituye un valle epistático para llegar al diseño de 13 genes. El GA sembrado logró cruzarlo inyectando los escalones intermedios. Sin embargo, el GA puro (sin seeding) **no cruza el valle**, sino que encuentra **soluciones alternativas compactas** (como el núcleo de 2 genes `ptsH+atpD` o el de 4 genes) que sortean el valle por otra ruta metabólica. Esto sugiere que el paisaje de aptitud admite varias soluciones compactas, aunque con 1 de 5 semillas la evidencia es limitada.

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
Nota sobre el modelo: iML1515.xml ya no se distribuye con el repositorio. Descárgalo de BiGG (bigg.ucsd.edu → iML1515 → SBML) y colócalo en la raíz del proyecto, o deja que la herramienta lo descargue la primera vez.

⚡ Uso rápido y Pipeline de Diseño
Menú interactivo:

bash
python metabo_engineer.py
CLI experto:

bash
python cli_experto.py --modelo iML1515 --quimico EX_succ_e --ambiente anaerobio --metodo ia --generaciones 20 --poblacion 50 --seed 42
Pipeline de Diseño (listo para el laboratorio): el script benchmark/pipeline_design.py automatiza el flujo de trabajo de un investigador (5 etapas: Validación → GA → Destilación LOO → Verificación biológica → Veredicto).

bash
# Uso estándar (descubre la biomasa automáticamente)
python benchmark/pipeline_design.py --modelo iML1515.xml --target EX_succ_e --ambiente anaerobio --seed 100

# Uso con seeding manual (para forzar genes conocidos)
python benchmark/pipeline_design.py --modelo iML1515.xml --target EX_succ_e --ambiente anaerobio --seeding b3738,b0002
⚠️ Limitaciones conocidas
Resultados in silico: ningún diseño publicado en este repositorio ha sido validado experimentalmente. Los valores son predicciones del modelo iML1515 en condiciones anaerobias con glucosa a 10 mmol/gDW/h.

Alcance multiorganismo: la herramienta está validada en E. coli y es operativa en otros organismos (levadura, B. subtilis, P. putida). Topologías metabólicas distintas pueden requerir ajustes en los parámetros del GA para superar óptimos locales.

BiGG puede fallar: el servidor tiene caídas intermitentes. Usa la descarga manual si es necesario.

OptKnock heurístico: usamos un GA y una enumeración sobre candidatos acotados, no una optimización bilevel exacta (MILP).

Dependencia de parámetros y azar: el GA puro encontró una solución compacta (núcleo de 2 genes) en 1 de 5 semillas con k=4, población 50 y 15 generaciones. Con presupuestos pequeños (k=3, población 20) converge a óptimos locales como la ATP sintasa.

Cinética ausente: los parámetros kcat/Km no están incluidos por defecto (requieren API de pago de BRENDA) y han sido removidos de los reportes automáticos para evitar confusiones.

Artefacto termodinámico en ATPS4r: el cálculo de eQuilibrator para la ATP sintasa (ATPS4r) da un ΔG' positivo (~52 kJ/mol) porque no incluye el potencial electroquímico de membrana. No es un error del código, sino una limitación metodológica del cálculo puramente estequiométrico.

eQuilibrator: requiere conexión a internet para descargar compuestos.

Margen de viabilidad: el núcleo de 4 genes tiene μ = 0.1071 (2.14× sobre el umbral de viabilidad de 0.05) y el de 2 genes μ = 0.0695 (1.4×). Con restricciones de glucosa severas (glucosa 5), el μ del núcleo de 4 genes cae a 0.0256.

🗺️ Roadmap (v0.3.0 y más allá)
☑ FBA + algoritmo genético
☑ OptKnock híbrido
☑ Termodinámica universal
☑ Reportes científicos
☑ Caché + indexación de modelos (reducción de tiempo 6 h → segundos)
☑ CLI experto con argparse
☑ Tests con pytest + CI
☑ Docker
☑ Filtro universal de esencialidad
☑ Simulación aerobia/anaerobia
☑ Validación gen-a-gen contra Keio
☑ Benchmark de solvers (GLPK, SCIP, HiGHS)
☑ Validación KJ060 in silico (familias de isozimas, k=13)
☑ Caracterización del valle de aptitud (5 intermedios, todos viables)
☑ Leave-one-out del diseño epistático (atpF + nrdD)
☑ GA con mutación real y reproducibilidad determinista
☑ GA puro k=4: solución compacta en 1 de 5 semillas (10.8049), reducida por destilación LOO a un núcleo de 2 genes (ptsH + atpD, 10.8056)
☑ Seeding + destilación LOO iterativa: núcleo de 4 genes (10.12)
☑ Suite de pruebas adversariales (0 fallos S0/S1 en los casos probados)
☑ Pipeline de diseño asistido (pipeline_design.py)
□ Validación con ≥30 semillas: ampliar el experimento del GA puro (k=4) para confirmar si el éxito del 1/5 es estadísticamente significativo.
□ Contraste con Jantama (2022): comparar los resultados en fase estacionaria (biomasa = 0) y microaerobiosis con el paper (0.48-0.52 mol/mol) para cerrar la brecha in silico/in vivo.
□ Análisis de flujos (MOMA/ROOM): evaluar la robustez del núcleo de 2 genes (ptsH + atpD) bajo MOMA/ROOM y sin bucles en FBA.
□ Leave-family-out del diseño minimal
□ Restricciones enzimáticas: integración con ECMpy (licencia MIT)
□ Publicación: preprint en bioRxiv documentando el benchmark, la validación y el valle de aptitud
📚 Citación y Predecesores
bibtex
@software{metaboengineer2025,
  author  = {San Miguel Crisantos, Diego},
  title   = {MetaboEngineer Pro: A Suite for Metabolic Engineering and Synthetic Biology},
  year    = {2025},
  version = {v0.2.8},
  url     = {https://github.com/diegosancriss16-tech/MetaboEngineer-Pro}
}
OptGene: Patil, K. R., et al. (2005). BMC Bioinformatics, 6(1), 308.

OptKnock: Burgard, A. P., et al. (2003). Biotechnology and Bioengineering, 84(6), 647-657.

COBRApy: Ebrahim, A., et al. (2013). BMC Systems Biology, 7(1), 74.

Keio Collection: Baba, T., et al. (2006). Molecular Systems Biology, 2(1), 2006-0008.

iML1515: Monk, J. M., et al. (2017). Nature Biotechnology, 35(10), 904-908.

Jantama et al.: Jantama, K., et al. (2008). Biotechnology and Bioengineering, 99(5), 1140-1153.

eQuilibrator: Beber, M. E., et al. (2022). eQuilibrator 3.0: a database solution for thermodynamic constant estimation. Nucleic Acids Research, 50(D1), D603-D609.

📄 Licencia
Este proyecto está licenciado bajo la PolyForm Noncommercial License 1.0.0.

Required Notice: Copyright Diego San Miguel Crisantos (https://github.com/diegosancriss16-tech/MetaboEngineer-Pro)

✅ Lo que puedes hacer gratis
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

⚠️ Nota sobre "Open Source": esta licencia es source-available, no OSI-Open-Source. El código es público, auditable y modificable, pero el uso comercial requiere permiso. Ver LICENSE para el texto legal completo.