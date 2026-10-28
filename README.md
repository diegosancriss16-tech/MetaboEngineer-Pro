# MetaboEngineer Pro 🧬

Herramienta de biología de sistemas de nivel profesional para el diseño de cepas microbianas. Utiliza Inteligencia Artificial (Algoritmos Genéticos), Termodinámica y Biología Sintética para predecir qué genes modificar y producir químicos de interés industrial.

## 🚀 Características Principales
- **Diseño por IA:** Algoritmo genético para maximizar la producción de químicos.
- **Biología Sintética:** Inyección automática de rutas heterólogas (ej. producción de cafeína).
- **Validación Termodinámica:** Integración con eQuilibrator para verificar la viabilidad energética de las rutas.
- **Análisis de Carga Proteica:** Penalización realista por el costo metabólico de fabricar enzimas.
- **OptKnock:** Búsqueda de acoplamiento de crecimiento (Growth Coupling).
- **Multi-Organismo y Multi-Repositorio:** Soporte para BiGG Models, BioModels y archivos SBML locales.
- **Traducción UniProt:** Identificación de nombres de genes y funciones proteicas reales.
- **Reportes Automáticos:** Generación de reportes en formato `.txt` y `.csv` para el laboratorio.

## 🛠️ Instalación
1. Clona este repositorio.
2. Crea un entorno virtual de Python:
   ```bash
   python -m venv venv


   # MetaboEngineer Pro: Plataforma de Diseño de Cepas e Ingeniería Metabólica *In Silico*

**Resumen**
MetaboEngineer Pro es una herramienta de línea de comandos desarrollada en Python para el diseño racional de cepas microbianas. Integra análisis de flujos metabólicos (FBA) mediante COBRApy, validación termodinámica con eQuilibrator, y algoritmos genéticos para la optimización de vías metabólicas. Además, incorpora un módulo de biología sintética para la inyección de rutas heterólogas y un motor de aprendizaje profundo que traduce genes a funciones proteicas usando la API de UniProt.

**Características Principales**
- **Multi-Organismo:** Soporta modelos de BiGG Models, BioModels y archivos locales (probado exitosamente en *E. coli*, *B. subtilis* y *C. glutamicum*).
- **Biología Sintética:** Capaz de introducir rutas de producción de compuestos no nativos, como la cafeína (vía Coffea arabica).
- **Optimización por IA:** Implementación de Algoritmos Genéticos para la búsqueda de knockouts que maximicen la producción, penalizando la carga proteica.
- **OptKnock Híbrido:** Búsqueda de acoplamiento de crecimiento para cepas industrialmente robustas.
- **Generación de Reportes:** Automatización de reportes de laboratorio en formatos `.txt` y `.csv` listos para experimentación.

**Resultados Demostrativos**
Durante las pruebas, la herramienta logró diseñar una cepa de *E. coli* MG1655 con una producción de cafeína de 8.34 mmol/gDW/h (ajustada por carga proteica), identificando los knockouts de `hemH`, `treA`, `wecE` y `b4088` como estrategia para redirigir el flujo de carbono y energía.

**Tecnologías Utilizadas:**
`Python 3.12`, `COBRApy`, `Equilibrator-API`, `Pandas`, `Scikit-learn`, `Rich` (UI), `UniProt REST API`.