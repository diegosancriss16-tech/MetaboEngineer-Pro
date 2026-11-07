import requests
import gzip
import shutil
from pathlib import Path

def descargar_modelo(model_id, url):
    print(f"Descargando {model_id}...")
    r = requests.get(url, allow_redirects=True)
    r.raise_for_status()
    
    ruta = Path(f"{model_id}.xml")
    
    # Si el servidor devuelve un archivo comprimido (gz)
    if r.headers.get('Content-Type') == 'application/gzip' or url.endswith('.gz'):
        with open(f"{model_id}.xml.gz", 'wb') as f:
            f.write(r.content)
        with gzip.open(f"{model_id}.xml.gz", 'rb') as f_in:
            with open(ruta, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        Path(f"{model_id}.xml.gz").unlink() # Borrar el .gz
    else:
        with open(ruta, 'wb') as f:
            f.write(r.content)
            
    print(f"✅ Guardado en: {ruta} ({ruta.stat().st_size} bytes)")

# Intentamos con Corynebacterium glutamicum
descargar_modelo("iCW773", "http://bigg.ucsd.edu/static/models/iCW773.xml")
