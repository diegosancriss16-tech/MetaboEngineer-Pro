import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))
import metabo_engineer as me
from cobra.io import read_sbml_model

def main():
    model = read_sbml_model("iML1515.xml")
    model = me.configurar_ambiente(model, "anaerobio")
    model.solver = 'glpk'
    biomass = next(r.id for r in model.reactions if "biomass" in r.id.lower())

    PASOS = {
        '1-KO  lactato':           ['b1380'],
        '4-KO  +etanol':           ['b1380','b1241','b1478','b0351'],
        '7-KO  +acetato':          ['b1380','b1241','b1478','b0351','b2296','b1849','b3115'],
        '9-KO  +formato':          ['b1380','b1241','b1478','b0351','b2296','b1849','b3115','b0904','b2492'],
        '13-KO +PFL (final)':      ['b1380','b1241','b1478','b0351','b2296','b1849','b3115','b0904','b2492','b0903','b3114','b3951','b3952'],
    }

    print(f"{'Intermedio del valle':24s}{'crecimiento':>12s}{'¿viable (>=0.05)?':>20s}")
    print("-" * 58)
    for n, genes in PASOS.items():
        m = model.copy()
        for g in genes:
            if g in m.genes: m.genes.get_by_id(g).knock_out()
        mu = m.slim_optimize()
        viable = "SÍ" if mu >= 0.05 else "NO (-1e6 para el GA)"
        print(f"{n:24s}{mu:12.4f}{viable:>20s}")

if __name__ == "__main__":
    main()
