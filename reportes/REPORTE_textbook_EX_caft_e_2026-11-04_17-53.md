# Reporte de Ingenieria Metabolica: EX_caft_e

**Fecha:** 2026-11-04 17:53
**Organismo Base:** `textbook`
**Reacciones totales del modelo:** 101

## 1. Resumen Ejecutivo

- **Rendimiento ajustado (Fitness):** `700.0000 mmol/gDW/h`
- **Genes knock-out:** `b2283, b4122, b3734, b0727`

## 2. Validacion Genomica (UniProt & Rhea)

| Gen (ID) | Nombre | Funcion Proteica | Rhea ID |
| :--- | :--- | :--- | :--- |
| `b2283` | **b2283** | No encontrado | `N/A` |
| `b4122` | **b4122** | No encontrado | `N/A` |
| `b3734` | **b3734** | No encontrado | `N/A` |
| `b0727` | **b0727** | No encontrado | `N/A` |

## 3. Nota sobre Cinetica Enzimatica

Los parametros cineticos especificos (kcat/Km) **no estan incluidos** en este reporte porque requieren una conexion autenticada a la API de BRENDA (acceso no publico). Esta herramienta prioriza la universalidad: todos los datos genomicos se extraen dinamicamente desde UniProt y Rhea, por lo que funciona para cualquier organismo.

Si necesitas cinetica, puedes aportar tus propios valores desde BRENDA/SABIO-RK para las reacciones relevantes.

## 4. Analisis Termodinamico (eQuilibrator)

| Reaccion | DG' (kJ/mol) | Viable |
| :--- | :--- | :--- |
| `EX_caft_e` | N/A (Sintetico) | N/A |
| `NADH16` | -2729.60 | SI |
| `FUM` | -3.43 | SI |
| `ATPS4r` | 52.11 | NO |
| `AKGDH` | -28.98 | SI |

---
*Generado automaticamente por MetaboEngineer Pro V5.8 - Universal*
