# ============================================
# MetaboEngineer Pro - Imagen Docker
# ============================================

FROM python:3.11-slim

LABEL maintainer="Diego San Cristobal"
LABEL description="MetaboEngineer Pro - Suite de Ingenieria Metabolica"
LABEL version="5.7"

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
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY . .

RUN mkdir -p modelos reportes logs && \
    chmod -R 777 modelos reportes logs

CMD ["python", "metabo_engineer.py"]