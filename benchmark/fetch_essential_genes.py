import pandas as pd
import requests
from pathlib import Path

def main():
    urls = {
        "iML1515": "https://raw.githubusercontent.com/MetNetComp/DBgDel/main/initial_remaining_gene/GE_genes_iML1515.csv",
        "iMM904": "https://raw.githubusercontent.com/MetNetComp/DBgDel/main/initial_remaining_gene/GE_genes_iMM904.csv"
    }
    out_dir = Path("data")
    out_dir.mkdir(exist_ok=True)

    for model_name, url in urls.items():
        print(f"Descargando datos de {model_name}...")
        try:
            df = pd.read_csv(url)
            out_file = out_dir / f"essential_{model_name}.csv"
            df.to_csv(out_file, index=False)
            print(f"  ✅ Guardado en: {out_file}")
            print(f"     Columnas: {df.columns.tolist()}")
        except Exception as e:
            print(f"  ❌ Error: {e}")

if __name__ == "__main__":
    main()