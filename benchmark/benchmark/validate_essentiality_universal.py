import argparse
import sys
import pandas as pd
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from cobra.io import read_sbml_model
from cobra.flux_analysis import single_gene_deletion
import metabo_engineer as me

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--model", required=True)
    p.add_argument("--data", required=True)
    p.add_argument("--default-medium", action="store_true")
    p.add_argument("--threshold", type=float, default=0.01)
    args = p.parse_args()

    print(f"Cargando modelo: {args.model}...")
    model = read_sbml_model(args.model)
    model.solver = 'glpk'

    if not args.default_medium:
        model = me.configurar_ambiente(model, "aerobio")
        model.solver = 'glpk'
    else:
        print("Usando medio por defecto.")

    wt_growth = model.slim_optimize()
    if wt_growth is None or wt_growth < 0.01:
        print(f"❌ El modelo no crece (WT = {wt_growth}).")
        return
    
    threshold_value = wt_growth * args.threshold
    print(f"  WT max growth: {wt_growth:.4f} 1/h")
    print(f"  Umbral de esencialidad: < {threshold_value:.4f} 1/h")

    print("\n--- Ejecutando single_gene_deletion (esto tarda)... ---")
    deletion_results = single_gene_deletion(model)
    
    predicted_essential = set()
    for idx, row in deletion_results.iterrows():
        gene_id = list(row['ids'])[0] if isinstance(row['ids'], (set, frozenset, list)) else row['ids']
        if row['growth'] < threshold_value:
            predicted_essential.add(gene_id)

    print(f"  Predichos esenciales (IDs del modelo): {len(predicted_essential)}")

    print(f"\n--- Cargando datos experimentales ---")
    df = pd.read_csv(args.data)
    print(f"  Columnas del CSV: {df.columns.tolist()}")

    # Normalizar columnas
    df.columns = [c.lower() for c in df.columns]
    gene_col = 'gene'
    ess_col = 'essential'

    if gene_col not in df.columns or ess_col not in df.columns:
        print("❌ El CSV debe tener columnas 'gene' y 'essential'.")
        return

    df[gene_col] = df[gene_col].astype(str).str.strip()
    df['is_essential_exp'] = df[ess_col].astype(int) == 1
    
    exp_essential_names = set(df[df['is_essential_exp']][gene_col])
    
    # MAPEO CLAVE: Convertir nombres de Keio (thrL) a IDs del modelo (b0001)
    gene_name_to_id = {g.name: g.id for g in model.genes}
    
    # Traducir los esenciales experimentales a IDs del modelo
    exp_essential_ids = set()
    for name in exp_essential_names:
        if name in gene_name_to_id:
            exp_essential_ids.add(gene_name_to_id[name])
        elif name in model.genes: # Si ya es un ID
            exp_essential_ids.add(name)
            
    # Traducir los NO esenciales experimentales a IDs del modelo
    exp_nonessential_names = set(df[~df['is_essential_exp']][gene_col])
    exp_nonessential_ids = set()
    for name in exp_nonessential_names:
        if name in gene_name_to_id:
            exp_nonessential_ids.add(gene_name_to_id[name])
        elif name in model.genes:
            exp_nonessential_ids.add(name)

    model_genes = set(model.genes.list_attr("id"))
    
    # Evaluar solo los genes que están en el modelo Y en el dataset
    all_exp_genes = exp_essential_ids.union(exp_nonessential_ids)
    common_genes = model_genes.intersection(all_exp_genes)
    
    if not common_genes:
        print("❌ No hay genes en común entre el modelo y el CSV (incluso después del mapeo).")
        return

    print(f"  Genes en común evaluados: {len(common_genes)}")

    pred_common = predicted_essential.intersection(common_genes)
    exp_ess_common = exp_essential_ids.intersection(common_genes)
    exp_noness_common = exp_nonessential_ids.intersection(common_genes)

    tp = len(pred_common.intersection(exp_ess_common))
    fp = len(pred_common.intersection(exp_noness_common))
    fn = len(exp_ess_common - pred_common)
    tn = len(exp_noness_common - pred_common)

    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0

    print("\n=== RESULTADOS GEN-A-GEN (REALES) ===")
    print(f"  TP: {tp}, FP: {fp}, FN: {fn}, TN: {tn}")
    print(f"  Precisión: {precision:.2%}")
    print(f"  Recall:    {recall:.2%}")
    print(f"  F1-Score:  {f1:.2%}")

if __name__ == "__main__":
    main()