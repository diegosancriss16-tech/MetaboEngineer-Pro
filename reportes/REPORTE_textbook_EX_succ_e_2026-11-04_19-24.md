# Reporte de Ingenieria Metabolica: EX_succ_e

**Fecha:** 2026-11-04 19:24
**Organismo Base:** `textbook`
**Reacciones totales del modelo:** 95

## 1. Resumen Ejecutivo

- **Rendimiento ajustado (Fitness):** `8.6498 mmol/gDW/h`
- **Genes knock-out:** `b2463, b1854, b0726`

## 2. Validacion Genomica (UniProt & Rhea)

| Gen (ID) | Nombre | Funcion Proteica | Rhea ID |
| :--- | :--- | :--- | :--- |
| `b2463` | **maeB** | NADP-dependent malic enzyme | `N/A` |
| `b1854` | **pykA** | Pyruvate kinase II | `N/A` |
| `b0726` | **sucA** | 2-oxoglutarate dehydrogenase E1 component | `N/A` |

## 3. Cinetica Enzimatica (Base de Datos Local)

Los siguientes parametros fueron extraidos de la base de datos local (`cinetica_local.json`).

| Reaccion | Enzima | kcat (1/s) | Km (mM) |
| :--- | :--- | :--- | :--- |
| `AKGDH` | Alfa-cetoglutarato deshidrogenasa | 120.0 | 0.1 |
> **Nota de diseno:** El software prioriza la universalidad sin requerir licencias de pago (ej. BRENDA). El investigador puede agregar valores manualmente al JSON.

## 4. Analisis Termodinamico (eQuilibrator)

| Reaccion | DG' (kJ/mol) | Viable |
| :--- | :--- | :--- |
| `ME2` | -4.63 | SI |
| `PYK` | -29.04 | SI |
| `AKGDH` | -28.98 | SI |

---
*Generado automaticamente por MetaboEngineer Pro V5.8 - Universal*
