import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import inspect
import metabo_engineer as me

print("=== FUNCIONES PUBLICAS DE metabo_engineer ===")
funcs = [n for n, obj in inspect.getmembers(me, inspect.isfunction)
         if not n.startswith("_")]
for f in sorted(funcs):
    try:
        sig = str(inspect.signature(getattr(me, f)))
        print(f"{f}{sig}")
    except (ValueError, TypeError):
        print(f"{f}(?)")

print("\n=== CLASES PUBLICAS ===")
clases = [n for n, obj in inspect.getmembers(me, inspect.isclass)
          if not n.startswith("_") and obj.__module__ == me.__name__]
for c in clases:
    print(f"class {c}")

print("\n=== MODULOS INTERNOS (archivos .py del proyecto) ===")
import pathlib
base = pathlib.Path(me.__file__).parent
for p in sorted(base.glob("*.py")):
    print(f"  {p.name}")
