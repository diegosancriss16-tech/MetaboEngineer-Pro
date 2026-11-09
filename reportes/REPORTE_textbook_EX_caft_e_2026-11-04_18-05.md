# Reporte de Ingenieria Metabolica: EX_caft_e

**Fecha:** 2026-11-04 18:05
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

## 3. Cinetica Enzimatica (Base de Datos Local)

Los siguientes parametros fueron extraidos de la base de datos local (`cinetica_local.json`).

| Reaccion | Enzima | kcat (1/s) | Km (mM) |
| :--- | :--- | :--- | :--- |
| `EX_caft_e` | Transporte de Cafeina | 0.0 | 0.0 |
| `NADH16` | NADH deshidrogenasa I | 500.0 | 0.05 |
| `FUM` | Fumarasa | 800.0 | 0.5 |
| `ATPS4r` | ATP sintasa | 100.0 | 0.2 |
| `AKGDH` | Alfa-cetoglutarato deshidrogenasa | 120.0 | 0.1 |
> **Nota de diseno:** El software prioriza la universalidad sin requerir licencias de pago (ej. BRENDA). El investigador puede agregar valores manualmente al JSON.

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
