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
    p.add_argument("--rich-medium", action="store_true", help="Abrir todos los exchanges (LB) para igualar Keio")
    p.add_argument("--threshold", type=float, default=0.01)
    args = p.parse_args()

    print(f"Cargando modelo: {args.model}...")
    model = read_sbml_model(args.model)
    model.solver = 'glpk'

    if args.rich_medium:
        print("Usando medio RICO (todos los exchanges abiertos excepto glucosa) para igualar Keio...")
        for ex in model.exchanges:
            if ex.id != 'EX_glc__D_e':
                ex.bounds = (-1000, 1000)
    elif not args.default_medium:
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

    print("\n--- Ejecutando single_gene_deletion... ---")
    deletion_results = single_gene_deletion(model)
    
    predicted_essential = set()
    for idx, row in deletion_results.iterrows():
        gene_id = list(row['ids'])[0] if isinstance(row['ids'], (set, frozenset, list)) else row['ids']
        if row['growth'] < threshold_value:
            predicted_essential.add(gene_id)

    print(f"  Predichos esenciales (IDs del modelo): {len(predicted_essential)}")

    print(f"\n--- Cargando datos experimentales ---")
    df = pd.read_csv(args.data)
    df.columns = [c.lower() for c in df.columns]
    
    df['gene'] = df['gene'].astype(str).str.strip()
    df['is_essential_exp'] = df['essential'].astype(int) == 1
    
    exp_essential_names = set(df[df['is_essential_exp']]['gene'])
    exp_nonessential_names = set(df[~df['is_essential_exp']]['gene'])
    
    gene_name_to_id = {g.name: g.id for g in model.genes}
    
    exp_essential_ids = set()
    for name in exp_essential_names:
        if name in gene_name_to_id: exp_essential_ids.add(gene_name_to_id[name])
        elif name in model.genes: exp_essential_ids.add(name)
            
    exp_nonessential_ids = set()
    for name in exp_nonessential_names:
        if name in gene_name_to_id: exp_nonessential_ids.add(gene_name_to_id[name])
        elif name in model.genes: exp_nonessential_ids.add(name)

    model_genes = set(model.genes.list_attr("id"))
    
    all_exp_genes = exp_essential_ids.union(exp_nonessential_ids)
    common_genes = model_genes.intersection(all_exp_genes)
    
    if not common_genes:
        print("❌ No hay genes en común.")
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
