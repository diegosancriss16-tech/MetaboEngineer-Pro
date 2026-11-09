# Reporte de Ingenieria Metabolica: EX_etoh_e

**Fecha:** 2026-11-04 18:10
**Organismo Base:** `iMM904`
**Reacciones totales del modelo:** 1577

## 1. Resumen Ejecutivo

- **Rendimiento ajustado (Fitness):** `15.0195 mmol/gDW/h`
- **Genes knock-out:** `YKR097W, YEL038W, YJL045W, YOR377W`

## 2. Validacion Genomica (UniProt & Rhea)

| Gen (ID) | Nombre | Funcion Proteica | Rhea ID |
| :--- | :--- | :--- | :--- |
| `YKR097W` | **PCK1** | Phosphoenolpyruvate carboxykinase (ATP) | `N/A` |
| `YEL038W` | **UTR4** | Enolase-phosphatase E1 | `N/A` |
| `YJL045W` | **Desconocido** | Succinate dehydrogenase [ubiquinone] flavoprotein subunit 2, mitochondrial | `N/A` |
| `YOR377W` | **ATF1** | Alcohol O-acetyltransferase 1 | `N/A` |

## 3. Cinetica Enzimatica (Base de Datos Local)

*No se encontraron parametros cineticos para estas reacciones en la base de datos local.*

> **Nota de diseno:** El software prioriza la universalidad sin requerir licencias de pago (ej. BRENDA). El investigador puede agregar valores manualmente al JSON.

## 4. Analisis Termodinamico (eQuilibrator)

| Reaccion | DG' (kJ/mol) | Viable |
| :--- | :--- | :--- |
| `PPCK` | -6.12 | SI |
| `DKMPPD2` | -227.83 | SI |
| `SUCD1m` | 38.60 | NO |
| `OHACT5` | -3.14 | SI |
| `OHACT2` | -3.14 | SI |
| `OHACT1` | -323.19 | SI |

---
*Generado automaticamente por MetaboEngineer Pro V5.8 - Universal*
