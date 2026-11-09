# Reporte de Ingenieria Metabolica: EX_succ_e

**Fecha:** 2026-11-04 19:23
**Organismo Base:** `iML1515`
**Reacciones totales del modelo:** 2712

## 1. Resumen Ejecutivo

- **Rendimiento ajustado (Fitness):** `9.7663 mmol/gDW/h`
- **Genes knock-out:** `b1702, b1676, b0728`

## 2. Validacion Genomica (UniProt & Rhea)

| Gen (ID) | Nombre | Funcion Proteica | Rhea ID |
| :--- | :--- | :--- | :--- |
| `b1702` | **ppsA** | Phosphoenolpyruvate synthase | `N/A` |
| `b1676` | **pykF** | Pyruvate kinase I | `N/A` |
| `b0728` | **sucC** | Succinate--CoA ligase [ADP-forming] subunit beta | `N/A` |

## 3. Cinetica Enzimatica (Base de Datos Local)

*No se encontraron parametros cineticos para estas reacciones en la base de datos local.*

> **Nota de diseno:** El software prioriza la universalidad sin requerir licencias de pago (ej. BRENDA). El investigador puede agregar valores manualmente al JSON.

## 4. Analisis Termodinamico (eQuilibrator)

| Reaccion | DG' (kJ/mol) | Viable |
| :--- | :--- | :--- |
| `PPS` | -14.53 | SI |
| `PYK` | -29.04 | SI |
| `SUCOAS` | -2.77 | SI |

---
*Generado automaticamente por MetaboEngineer Pro V5.8 - Universal*
