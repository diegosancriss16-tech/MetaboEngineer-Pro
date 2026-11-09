import pandas as pd
from pathlib import Path

def main():
    data_dir = Path("data")
    
    print("Cargando listas de genes...")
    df_all = pd.read_csv(data_dir / "keio_complete.csv")
    df_essential = pd.read_csv(data_dir / "keio_essential.csv")
    
    df_all.columns = [c.lower() for c in df_all.columns]
    df_essential.columns = [c.lower() for c in df_essential.columns]
    
    all_genes = set(df_all['gene'].astype(str).str.strip())
    essential_genes = set(df_essential['gene'].astype(str).str.strip())
    
    print(f"Genes totales en Keio: {len(all_genes)}")
    print(f"Genes esenciales en Keio: {len(essential_genes)}")
    
    records = []
    for gene in all_genes:
        is_essential = 1 if gene in essential_genes else 0
        records.append({"gene": gene, "essential": is_essential})
    
    df_labeled = pd.DataFrame(records)
    out_file = data_dir / "keio_labeled.csv"
    df_labeled.to_csv(out_file, index=False)
    
    print(f"\n✅ Dataset etiquetado guardado en: {out_file}")
    print(f"   Total de genes: {len(df_labeled)}")
    print(f"   Esenciales (1): {df_labeled['essential'].sum()}")
    print(f"   No esenciales (0): {len(df_labeled) - df_labeled['essential'].sum()}")

if __name__ == "__main__":
    main()
