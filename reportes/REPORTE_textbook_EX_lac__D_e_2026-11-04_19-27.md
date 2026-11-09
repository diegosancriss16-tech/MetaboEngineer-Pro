# Reporte de Ingenieria Metabolica: EX_lac__D_e

**Fecha:** 2026-11-04 19:27
**Organismo Base:** `textbook`
**Reacciones totales del modelo:** 95

## 1. Resumen Ejecutivo

- **Rendimiento ajustado (Fitness):** `12.9361 mmol/gDW/h`
- **Genes knock-out:** `b0727, b0724, b0729, b0721`

## 2. Validacion Genomica (UniProt & Rhea)

| Gen (ID) | Nombre | Funcion Proteica | Rhea ID |
| :--- | :--- | :--- | :--- |
| `b0727` | **sucB** | Dihydrolipoyllysine-residue succinyltransferase component of 2-oxoglutarate dehydrogenase complex | `N/A` |
| `b0724` | **sdhB** | Succinate dehydrogenase iron-sulfur subunit | `N/A` |
| `b0729` | **sucD** | Succinate--CoA ligase [ADP-forming] subunit alpha | `N/A` |
| `b0721` | **sdhC** | Succinate dehydrogenase cytochrome b556 subunit | `N/A` |

## 3. Cinetica Enzimatica (Base de Datos Local)

Los siguientes parametros fueron extraidos de la base de datos local (`cinetica_local.json`).

| Reaccion | Enzima | kcat (1/s) | Km (mM) |
| :--- | :--- | :--- | :--- |
| `AKGDH` | Alfa-cetoglutarato deshidrogenasa | 120.0 | 0.1 |
> **Nota de diseno:** El software prioriza la universalidad sin requerir licencias de pago (ej. BRENDA). El investigador puede agregar valores manualmente al JSON.

## 4. Analisis Termodinamico (eQuilibrator)

| Reaccion | DG' (kJ/mol) | Viable |
| :--- | :--- | :--- |
| `AKGDH` | -28.98 | SI |
| `SUCDi` | -2667.48 | SI |
| `SUCOAS` | -2.00 | SI |

---
*Generado automaticamente por MetaboEngineer Pro V5.8 - Universal*
