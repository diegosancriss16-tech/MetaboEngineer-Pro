import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import metabo_engineer as me
from cobra.io import read_sbml_model

def main():
    model = read_sbml_model("iML1515.xml")
    model = me.configurar_ambiente(model, "anaerobio")
    model.solver = 'glpk'
    
    KJ060 = ['ldhA', 'adhE', 'ackA', 'focA', 'pflB']
    id_por_nombre = {g.name: g.id for g in model.genes if g.name}
    
    print("="*60)
    print("DIAGNÓSTICO DE MAPEO KJ060")
    print("="*60)
    
    for n in KJ060:
        if n not in id_por_nombre:
            print(f"❌ {n} NO ESTÁ en el modelo.")
            continue
            
        gid = id_por_nombre[n]
        g = model.genes.get_by_id(gid)
        print(f"\n🔹 {n} -> {gid} (afecta a {len(g.reactions)} reacciones):")
        
        for r in g.reactions:
            print(f"   {r.id:15s} bounds={str(r.bounds):15s} GPR: {r.gene_reaction_rule}")

    # Comprobar qué se bloquea realmente
    print("\n" + "="*60)
    print("REACCIONES BLOQUEADAS AL APLICAR LOS KO")
    print("="*60)
    
    ref = {r.id: r.bounds for r in model.reactions}
    m = model.copy()
    
    for n in KJ060:
        if n in id_por_nombre:
            m.genes.get_by_id(id_por_nombre[n]).knock_out()
            
    blocked = 0
    for r in m.reactions:
        if r.bounds != ref[r.id]:
            print(f"   {r.id:15s}: {ref[r.id]} -> {r.bounds}")
            blocked += 1
            
    print(f"\nTotal de reacciones bloqueadas: {blocked}")

if __name__ == "__main__":
    main()
