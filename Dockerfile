FROM python:3.10-slim

# Instalar dependencias del sistema (GLPK, SWIG, libexpat1 y libxml2 para SBML)
RUN apt-get update && apt-get install -y \
    glpk-utils \
    libglpk-dev \
    swig \
    gcc \
    g++ \
    libexpat1 \
    libexpat1-dev \
    libxml2 \
    libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Copiar e instalar dependencias de Python con timeout alto para evitar cortes de red
COPY requirements.txt .
RUN pip install --no-cache-dir --timeout 1000 --retries 10 -r requirements.txt
RUN pip install --no-cache-dir --timeout 1000 --retries 10 equilibrator-api

# Copiar el resto del proyecto
COPY . .

# Comando por defecto
CMD ["python", "metabo_engineer.py"]
