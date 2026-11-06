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

    print(f"  Predichos esenciales: {len(predicted_essential)}")

    print(f"\n--- Cargando datos experimentales ---")
    df = pd.read_csv(args.data)
    
    # Imprimir las primeras filas y columnas para diagnosticar
    print(f"  Columnas del CSV: {df.columns.tolist()}")
    print(f"  Primeras 3 filas:\n{df.head(3).to_string()}")

    # Lógica universal para encontrar la columna del gen
    gene_col = next((c for c in df.columns if c.lower() in ['gene', 'gene_id', 'id', 'b_number', 'locus']), df.columns[0])
    
    # Lógica universal para encontrar la columna de esencialidad
    ess_col = next((c for c in df.columns if c.lower() in ['essential', 'essentiality', 'is_essential', 'status']), None)

    df = df.dropna(subset=[gene_col])
    df[gene_col] = df[gene_col].astype(str).str.strip()
    
    # Si hay columna de esencialidad, la usamos. Si no, asumimos que TODOS los genes en el CSV son esenciales.
    if ess_col:
        print(f"  Usando columna de esencialidad: '{ess_col}'")
        df['is_essential_exp'] = df[ess_col].apply(lambda x: str(x).lower() in ['1', 'true', 'yes', 'essential', 'e'])
    else:
        print("  ⚠️ No hay columna de esencialidad explícita. Asumiendo que TODOS los genes en el CSV son esenciales.")
        df['is_essential_exp'] = True
    
    exp_essential = set(df[df['is_essential_exp']][gene_col])
    print(f"  Genes experimentales esenciales (en el CSV): {len(exp_essential)}")
    
    model_genes = set(model.genes.list_attr("id"))
    common_genes = model_genes.intersection(set(df[gene_col]))
    
    if not common_genes:
        print("❌ No hay genes en común entre el modelo y el CSV.")
        print(f"   Ejemplo de gen en modelo: {list(model_genes)[0]}")
        print(f"   Ejemplo de gen en CSV: {list(df[gene_col])[0]}")
        return

    pred_common = predicted_essential.intersection(common_genes)
    exp_common = exp_essential.intersection(common_genes)

    tp = len(pred_common.intersection(exp_common))
    fp = len(pred_common - exp_common)
    fn = len(exp_common - pred_common)
    tn = len(common_genes) - tp - fp - fn

    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0

    print("\n=== RESULTADOS GEN-A-GEN ===")
    print(f"  Genes en común evaluados: {len(common_genes)}")
    print(f"  TP: {tp}, FP: {fp}, FN: {fn}, TN: {tn}")
    print(f"  Precisión: {precision:.2%}")
    print(f"  Recall:    {recall:.2%}")
    print(f"  F1-Score:  {f1:.2%}")

if __name__ == "__main__":
    main()