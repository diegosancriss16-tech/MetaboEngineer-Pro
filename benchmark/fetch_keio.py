import pandas as pd
import requests
import io
from pathlib import Path

def main():
    print("Descargando dataset de esencialidad (Keio Collection / COBRApy)...")
    # URL estable del repositorio oficial de COBRApy
    url = "https://raw.githubusercontent.com/opencobra/cobrapy/main/tests/data/ecoli_essential_genes.csv"
    
    try:
        response = requests.get(url)
        response.raise_for_status()
        
        df = pd.read_csv(io.StringIO(response.text))
        
        out_path = Path("data")
        out_path.mkdir(exist_ok=True)
        out_file = out_path / "keio_monk_essential.csv"
        df.to_csv(out_file, index=False)
        
        print(f"✅ Dataset guardado en: {out_file}")
        print(f"   Total de genes evaluados: {len(df)}")
        
        # Imprimir las columnas para saber cómo se llaman
        print(f"   Columnas encontradas: {df.columns.tolist()}")
        
    except Exception as e:
        print(f"❌ Error: {e}")

if __name__ == "__main__":
    main()