# -*- coding: utf-8 -*-
r"""
PRUEBAS ADVERSARIALES Fase 1 - v0.2.1 (Actualizado para reconocer S3)
"""
import sys, time, traceback, io
from contextlib import redirect_stdout
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import cobra
from optlang import available_solvers
if "glpk" in available_solvers:
    cobra.Configuration().solver = "glpk"
    import cobra.util.solver as _cus
    _cus.solvers.pop("gurobi", None)
import metabo_engineer as me

BIOMASA = "BIOMASS_Ec_iML1515_core_75p37M"

def preparado():
    m0 = me.cargar_modelo_universal("iML1515")
    m = m0[0] if isinstance(m0, (tuple, list)) else m0
    try:
        m.solver = "glpk"
    except Exception:
        pass
    r = me.configurar_ambiente(m, "anaerobio")
    return r if isinstance(r, cobra.Model) else m

RESULTADOS = []

def clasificar(excepcion, retorno):
    if excepcion is not None:
        # NUEVO: Si es un ValueError, es un rechazo controlado (S3)
        if isinstance(excepcion, ValueError):
            return "S3", f"RECHAZO LIMPIO: {str(excepcion)[:70]}"
        # Si es otra excepción (KeyError, etc.), es un Crash (S0)
        return "S0", f"CRASH: {type(excepcion).__name__}: {str(excepcion)[:70]}"
    
    if retorno is None:
        return "S2", "devolvio None (fallo atrapado)"
        
    return "S1", f"SILENCIO: devolvio {retorno!r} sin advertir"

def adversarial(id_, descripcion, prediccion, func, es_valido=False):
    t0 = time.time()
    excepcion, retorno, salida = None, None, ""
    try:
        buf = io.StringIO()
        with redirect_stdout(buf):
            retorno = func()
        salida = buf.getvalue()
    except KeyboardInterrupt:
        raise
    except Exception as e:
        excepcion = e
    
    # Si el caso era válido (ej. A06) y funcionó, es un S3 (aceptado correctamente)
    if es_valido and excepcion is None:
        sev, detalle = "S3", "ACEPTACIÓN LIMPIA: procesó la entrada correctamente"
    else:
        sev, detalle = clasificar(excepcion, retorno)
        
    if sev == "S1" and any(k in salida.lower() for k in ("no valido", "error", "no existe", "invalid", "warning", "no se encontro")):
        sev = "S2"
        
    RESULTADOS.append((id_, descripcion, prediccion, sev, detalle, time.time()-t0))
    icono = {"S3": "[S3 OK ]", "S2": "[S2 ok ]", "S1": "[S1 !! ]", "S0": "[S0 X  ]"}[sev]
    print(f"{icono} {id_} {descripcion}")
    print(f"       predicho: {prediccion} | obtenido: {sev} ({detalle[:90]})")

if __name__ == "__main__":
    print("=" * 66)
    print(" ADVERSARIALES FASE 1: ENTRADAS INVALIDAS")
    print("=" * 66 + "\n")
    m = preparado()
    print(f"modelo listo: {len(m.reactions)} rxns, anaerobio\n" + "-" * 66)

    adversarial("A01", "quimico inexistente EX_foo_e", "S3",
        lambda: me.evaluar_bacteria(m.copy(), ["b3731"], "EX_foo_e", BIOMASA))
    
    adversarial("A02", "gen inexistente b9999", "S3",
        lambda: me.evaluar_bacteria(m.copy(), ["b9999"], "EX_succ_e", BIOMASA))
    
    adversarial("A03", "lista de genes vacia []", "S3",
        lambda: me.evaluar_bacteria(m.copy(), [], "EX_succ_e", BIOMASA))
    
    adversarial("A04", "reaccion de biomasa falsa", "S3",
        lambda: me.evaluar_bacteria(m.copy(), ["b3731"], "EX_succ_e", "BIOMASS_falsa_XX"))
    
    adversarial("A05", "ambiente 'marciano'", "S3",
        lambda: me.configurar_ambiente(m.copy(), "marciano"))
    
    # A06 es válido (mayúsculas), así que esperamos S3
    adversarial("A06", "ambiente 'ANAEROBIO' (mayusculas)", "S3",
        lambda: me.configurar_ambiente(m.copy(), "ANAEROBIO"), es_valido=True)
    
    adversarial("A07", "gen con espacio 'b3731 '", "S3",
        lambda: me.evaluar_bacteria(m.copy(), ["b3731 "], "EX_succ_e", BIOMASA))
    
    adversarial("A08", "quimico sin prefijo 'succ' en mapa", "S2",
        lambda: me.generar_mapa_universal(modelo_id="iML1515", quimico="succ"))
    
    adversarial("A09", "modelo 'iNoExisto999'", "S3",
        lambda: me.cargar_modelo_universal("iNoExisto999"))
    
    adversarial("A10", "gen E.coli (b3731) en levadura iMM904", "S3",
        lambda: (lambda mm: me.evaluar_bacteria(
            (lambda r: r if isinstance(r, cobra.Model) else mm)(
                me.configurar_ambiente(mm, "anaerobio") if True else mm),
            ["b3731"], "EX_succ_e", "BIOMASS_SC5_notrace"))(
            (lambda crudo: crudo[0] if isinstance(crudo, (tuple, list)) else crudo)(
                me.cargar_modelo_universal("iMM904"))))

    print("\n" + "=" * 66)
    print(" RESUMEN FASE 1")
    print("=" * 66)
    conteo = {}
    for r in RESULTADOS:
        conteo[r[3]] = conteo.get(r[3], 0) + 1
    print(f"{'ID':<5}{'Descripcion':<38}{'Pred.':<8}{'Real':<6}")
    print("-" * 66)
    for id_, desc, pred, sev, _, dur in RESULTADOS:
        print(f"{id_:<5}{desc[:37]:<38}{pred:<8}{sev:<6}")
    print("-" * 66)
    print(f"S3 (rechazo limpio): {conteo.get('S3',0)} | S2 (degrada): {conteo.get('S2',0)} | "
          f"S1 (SILENCIO): {conteo.get('S1',0)} | S0 (crash): {conteo.get('S0',0)}")
    s0s1 = conteo.get("S0", 0) + conteo.get("S1", 0)
    print(f"\nTICKETS A ABRIR (S0+S1): {s0s1}")
    sys.exit(0)
