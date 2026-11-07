FROM python:3.10-slim

<<<<<<< HEAD
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
=======
FROM python:3.11-slim

LABEL maintainer="Diego San Cristobal"
LABEL description="MetaboEngineer Pro - Suite de Ingenieria Metabolica"
LABEL version="0.2.7"

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        libglpk-dev \
        glpk-utils \
        git \
        curl \
>>>>>>> 28960a205200ccc85657bf1d171a0525c8aed74b
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
