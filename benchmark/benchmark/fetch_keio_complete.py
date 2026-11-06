import pandas as pd
import requests
from pathlib import Path

def main():
    print("Descargando dataset COMPLETO de esencialidad (Monk 2017)...")
    # Repositorio oficial de Monk et al. 2017 (E. coli)
    url = "https://raw.githubusercontent.com/opencobra/cobrapy/develop/tests/data/ecoli_essential_genes.csv"
    
    try:
        response = requests.get(url)
        response.raise_for_status()
        df = pd.read_csv(url)
        
        out_path = Path("data")
        out_path.mkdir(exist_ok=True)
        out_file = out_path / "keio_monk_complete.csv"
        df.to_csv(out_file, index=False)
        
        print(f"✅ Guardado en: {out_file}")
        print(f"   Columnas: {df.columns.tolist()}")
        print(f"   Total de genes: {len(df)}")
        
        # Verificar si hay negativos
        if 'essential' in df.columns:
            print(f"   Esenciales: {df['essential'].sum()}")
            print(f"   No esenciales: {len(df) - df['essential'].sum()}")
        
    except Exception as e:
        print(f"❌ Error: {e}")
        print("Intentando URL alternativa...")
        # Fallback: intentar desde el repositorio de Monk
        url2 = "https://raw.githubusercontent.com/opencobra/cobrapy/main/tests/data/ecoli_essential_genes.csv"
        try:
            df = pd.read_csv(url2)
            df.to_csv(out_file, index=False)
            print(f"✅ Guardado desde URL alternativa: {out_file}")
        except Exception as e2:
            print(f"❌ Error alternativo: {e2}")

if __name__ == "__main__":
    main()