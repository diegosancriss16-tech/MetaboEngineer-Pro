# Reporte de Ingenieria Metabolica: EX_etoh_e

**Fecha:** 2026-11-04 19:28
**Organismo Base:** `textbook`
**Reacciones totales del modelo:** 95

## 1. Resumen Ejecutivo

- **Rendimiento ajustado (Fitness):** `0.0000 mmol/gDW/h`
- **Genes knock-out:** `b3916, b2463, b1779, b0726`

## 2. Validacion Genomica (UniProt & Rhea)

| Gen (ID) | Nombre | Funcion Proteica | Rhea ID |
| :--- | :--- | :--- | :--- |
| `b3916` | **pfkA** | ATP-dependent 6-phosphofructokinase isozyme 1 | `N/A` |
| `b2463` | **maeB** | NADP-dependent malic enzyme | `N/A` |
| `b1779` | **gapA** | Glyceraldehyde-3-phosphate dehydrogenase A | `N/A` |
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
| `PFK` | -15.56 | SI |
| `ME2` | -4.63 | SI |
| `GAPD` | -32.57 | SI |
| `AKGDH` | -28.98 | SI |

---
*Generado automaticamente por MetaboEngineer Pro V5.8 - Universal*
