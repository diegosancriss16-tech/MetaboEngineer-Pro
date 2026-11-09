# Reporte de Ingenieria Metabolica: EX_succ_e

**Fecha:** 2026-11-04 19:04
**Organismo Base:** `iML1515`
**Reacciones totales del modelo:** 2712

## 1. Resumen Ejecutivo

- **Rendimiento ajustado (Fitness):** `9.7663 mmol/gDW/h`
- **Genes knock-out:** `b4101, b0595, b2423, b1124, b2093`

## 2. Validacion Genomica (UniProt & Rhea)

| Gen (ID) | Nombre | Funcion Proteica | Rhea ID |
| :--- | :--- | :--- | :--- |
| `b4101` | **phnG** | Alpha-D-ribose 1-methylphosphonate 5-triphosphate synthase subunit PhnG | `N/A` |
| `b0595` | **entB** | Enterobactin synthase component B | `N/A` |
| `b2423` | **cysW** | Sulfate transport system permease protein CysW | `N/A` |
| `b1124` | **potC** | Spermidine/putrescine transport system permease protein PotC | `N/A` |
| `b2093` | **gatB** | PTS system galactitol-specific EIIB component | `N/A` |

## 3. Cinetica Enzimatica (Base de Datos Local)

*No se encontraron parametros cineticos para estas reacciones en la base de datos local.*

> **Nota de diseno:** El software prioriza la universalidad sin requerir licencias de pago (ej. BRENDA). El investigador puede agregar valores manualmente al JSON.

## 4. Analisis Termodinamico (eQuilibrator)

| Reaccion | DG' (kJ/mol) | Viable |
| :--- | :--- | :--- |
| `ICHORT` | -58.68 | SI |
| `MOBDabcpp` | -43.67 | SI |
| `SULabcpp` | -43.67 | SI |
| `TSULabcpp` | -43.67 | SI |
| `SLNTabcpp` | -43.67 | SI |
| `SELabcpp` | -43.67 | SI |
| `SPMDabcpp` | -43.67 | SI |
| `PTRCabcpp` | -43.67 | SI |

---
*Generado automaticamente por MetaboEngineer Pro V5.8 - Universal*
