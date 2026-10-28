
# MetaboEngineer Pro 🧬

Herramienta de biología de sistemas de nivel profesional para el diseño de cepas microbianas. Utiliza Inteligencia Artificial (Algoritmos Genéticos), Termodinámica y Biología Sintética para predecir qué genes modificar y producir químicos de interés industrial.

## 🚀 Características Principales
- **Diseño por IA:** Algoritmo genético para maximizar la producción de químicos, considerando una penalización por carga proteica.
- **Biología Sintética:** Inyección automática de rutas heterólogas (ej. producción de cafeína en *E. coli* vía *Coffea arabica*).
- **Validación Termodinámica:** Integración con `eQuilibrator` para verificar la viabilidad energética de las rutas metabólicas.
- **OptKnock:** Búsqueda de acoplamiento de crecimiento (Growth Coupling) para diseñar cepas industrialmente robustas.
- **Multi-Organismo y Multi-Repositorio:** Soporte para BiGG Models, BioModels y archivos SBML locales (probado en *E. coli*, *B. subtilis* y *C. glutamicum*).
- **Traducción UniProt & Rhea:** Identificación de nombres reales de genes, funciones proteicas e IDs de reacciones bioquímicas.
- **Reportes Automáticos:** Generación de reportes científicos en formato Markdown (`.md`) y mapas metabólicos interactivos (`.html`).

---

## 🛠️ Requisitos Previos
Antes de instalar, asegúrate de tener:
1. **Python 3.9 o superior** instalado en tu computadora.
2. **Git** instalado (para clonar el repositorio).
3. Una terminal o línea de comandos (CMD, PowerShell o Terminal).

---

## 📥 Instalación Paso a Paso

**1. Clona este repositorio:**
```bash
git clone https://github.com/diegosacriss16-tech/MetaboEngineer-Pro.git
cd MetaboEngineer-Pro
```

**2. Crea un entorno virtual (recomendado):**
*En Windows:*
```bash
python -m venv venv
venv\Scripts\activate
```
*En Mac/Linux:*
```bash
python3 -m venv venv
source venv/bin/activate
```

**3. Instala las dependencias necesarias:**
```bash
pip install -r requirements.txt
```

---

## 💻 Uso

Para iniciar la herramienta, ejecuta el siguiente comando en tu terminal:
```bash
python metabo_engineer.py
```

Verás un menú interactivo en la consola:
1. **Análisis de Flujos Básico (FBA):** Calcula la tasa de crecimiento del organismo.
2. **Diseño de Cepas por IA:** Ejecuta el algoritmo genético para encontrar los mejores knockouts.
3. **Escaneo de Letalidad:** Busca combinaciones de genes que, al apagarse, matan al organismo.
4. **Análisis Integral:** Ejecuta el diseño por IA, la termodinámica, genera el mapa visual y el reporte automático.

### 📂 Estructura de Carpetas
La herramienta creará automáticamente dos carpetas cuando la ejecutes por primera vez:
- `modelos/`: Aquí se guardan los modelos metabólicos (`.xml`) descargados.
- `reportes/`: Aquí se guardan los reportes científicos (`.md`) y los mapas visuales (`.html`).

---

## ⚠️ Solución de Problemas Comunes

**1. Error "Servidor BiGG no responde"**
A veces, los servidores de BiGG Models se saturan. Si intentas cargar un modelo grande (como `iML1515`) y falla, descarga el archivo manualmente:
- Ve a [BiGG Models](http://bigg.ucsd.edu/models).
- Busca el ID del modelo (ej. `iML1515`).
- Descarga el archivo `.xml` y guárdalo dentro de la carpeta `modelos/` de tu proyecto.
- Vuelve a ejecutar el programa; lo leerá desde tu disco duro.

**2. Error con la Termodinámica (eQuilibrator)**
La primera vez que uses el módulo termodinámico, el programa descargará una base de datos de 1.3 GB en segundo plano (esto puede tardar unos minutos). Por favor, **no cierres la terminal** hasta que termine de descargar.

---

## 🧪 Resultados Demostrativos
Durante las pruebas, la herramienta logró diseñar una cepa de *E. coli* MG1655 con una producción teórica de cafeína de **8.34 mmol/gDW/h** (ajustada por carga proteica), identificando los knockouts de `hemH`, `treA`, `wecE` y `b4088` como estrategia para redirigir el flujo de carbono.

## 🛠️ Tecnologías Utilizadas
`Python 3.12`, `COBRApy` (Análisis de flujos), `Equilibrator-API` (Termodinámica), `Pandas` (Manejo de datos), `Scikit-learn` (Modelos predictivos), `Rich` (Interfaz de terminal), `NetworkX` + `PyVis` (Visualización de redes), `UniProt REST API` (Genómica).
```