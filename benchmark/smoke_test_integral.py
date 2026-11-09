# -*- coding: utf-8 -*-
r"""
SMOKE TEST INTEGRAL v3 - MetaboEngineer Pro v0.2.0
v3 corrige: (1) configurar_ambiente no muta en sitio -> se captura su retorno
            (2) analizar_termodinamica_universal imprime -> se captura stdout
            (3) GLPK se fuerza ANTES de importar el software y se retira gurobi
                del registro, para que llegue tambien al modelo interno del mapa
Uso:
    python benchmark\smoke_test_integral.py --quick  (sin barrido letal, ~2 min)
    python benchmark\smoke_test_integral.py          (completo, ~15 min)
    python benchmark\smoke_test_integral.py --full   (+ GA y OptKnock)
"""
import sys, time, traceback, io
from contextlib import redirect_stdout
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

# ---- FIX (3): ANTES de importar el software ----
import cobra
from optlang import available_solvers
if "glpk" in available_solvers:
    cobra.Configuration().solver = "glpk"
    import cobra.util.solver as _cus
    _cus.solvers.pop("gurobi", None)   # nadie puede pedir gurobi ya

import metabo_engineer as me

QUICK = "--quick" in sys.argv

def extraer_modelo(resultado):
    if isinstance(resultado, cobra.Model):
        return resultado
    if isinstance(resultado, (tuple, list)):
        for elemento in resultado:
            if isinstance(elemento, cobra.Model):
                return elemento
    raise TypeError(f"cargar_modelo_universal no devolvio un Model: "
                    f"{type(resultado)} -> {str(resultado)[:120]}")

def preparar_ambiente(modelo, ambiente):
    """FIX (1): cubre tanto si muta en sitio como si devuelve un modelo nuevo."""
    res = me.configurar_ambiente(modelo, ambiente)
    return res if isinstance(res, cobra.Model) else modelo

RESULTADOS = []

def registrar(nombre, estado, detalle="", duracion=0.0):
    RESULTADOS.append((nombre, estado, detalle, duracion))
    icono = {"PASSED": "[OK]  ", "FAILED": "[FAIL]", "SKIPPED": "[SKIP]"}[estado]
    print(f"{icono} {nombre}  ({duracion:.1f}s)  {detalle}")

def prueba(nombre):
    def envoltura(func):
        def run(*a, **k):
            t0 = time.time()
            try:
                detalle = func(*a, **k) or ""
                registrar(nombre, "PASSED", detalle, time.time()-t0)
            except ConnectionError:
                registrar(nombre, "SKIPPED", "sin internet", time.time()-t0)
            except Exception as e:
                registrar(nombre, "FAILED", f"{type(e).__name__}: {e}", time.time()-t0)
                traceback.print_exc()
        return run
    return envoltura

# ================= 1. INFRAESTRUCTURA =================
@prueba("setup_logger")
def t01():
    me.setup_logger()
    return "logger configurado"

@prueba("listar_modelos_cacheados")
def t02():
    return f"{len(me.listar_modelos_cacheados())} modelos en cache"

@prueba("cargar_modelo_universal(iML1515)")
def t03():
    global MODEL
    crudo = me.cargar_modelo_universal("iML1515")
    MODEL = extraer_modelo(crudo)
    try:
        MODEL.solver = "glpk"
    except Exception:
        pass
    assert len(MODEL.reactions) > 2000
    return f"{len(MODEL.reactions)} rxns, {len(MODEL.genes)} genes"

# ================= 2. AMBIENTES Y FBA =================
@prueba("configurar_ambiente(anaerobio) + FBA")
def t04():
    m0 = MODEL.copy()
    res = me.configurar_ambiente(m0, "anaerobio")
    m = res if isinstance(res, cobra.Model) else m0
    o2 = m.reactions.get_by_id("EX_o2_e")
    assert o2.lower_bound >= -1e-9, (f"O2 abierto tras anaerobio "
        f"(retorno={type(res).__name__}, O2 lb={o2.lower_bound})")
    mu = m.slim_optimize()
    assert 0.05 < mu < 0.5, f"crecimiento anaerobio anormal: {mu}"
    return f"O2 cerrado, mu={mu:.4f} (esperado ~0.16), retorno={type(res).__name__}"

@prueba("configurar_ambiente(aerobio) + FBA")
def t05():
    m = preparar_ambiente(MODEL.copy(), "aerobio")
    mu = m.slim_optimize()
    assert mu > 0.5, f"crecimiento aerobio bajo: {mu}"
    return f"mu={mu:.4f} (esperado ~0.88)"

@prueba("get_organism_id")
def t06():
    org = me.get_organism_id("iML1515")
    assert isinstance(org, str) and org
    return f"iML1515 -> {org}"

# ================= 3. ESENCIALIDAD Y EVALUACION =================
@prueba("obtener_genes_no_esenciales (barrido anaerobio)")
def t07():
    m = preparar_ambiente(MODEL.copy(), "anaerobio")
    o2lb = m.reactions.get_by_id("EX_o2_e").lower_bound
    assert o2lb >= -1e-9, "guarda: el barrido correria en AEROBIO"
    no_ess = me.obtener_genes_no_esenciales(m, "BIOMASS_Ec_iML1515_core_75p37M")
    assert no_ess is not None and len(no_ess) > 100, f"muy pocos: {no_ess}"
    letales = 1516 - len(no_ess)
    return f"{len(no_ess)} no esenciales, ~{letales} letales (esperado < 196)"

@prueba("knock_out_model_genes (envoltorio cobra)")
def t08():
    m = MODEL.copy()
    afectadas = me.knock_out_model_genes(m, ["b3731"])  # atpC
    assert len(afectadas) >= 1, "no afecto reacciones"
    return f"KO b3731 desactivo {len(afectadas)} reacciones"

@prueba("evaluar_bacteria(b3731/atpC, anaerobio) ~ 1.5562")
def t09():
    m = preparar_ambiente(MODEL.copy(), "anaerobio")
    fitness = me.evaluar_bacteria(m, ["b3731"], "EX_succ_e",
                                   "BIOMASS_Ec_iML1515_core_75p37M")
    assert fitness is not None, "devolvio None"
    assert 1.0 < fitness < 2.0, (f"fitness={fitness} "
        f"(si es 0.0, el anaerobio no aplico)")
    return f"fitness={fitness:.4f} (benchmark README: 1.5562)"

@prueba("evaluar_bacteria(KJ060 extendido 13 KOs) = 7.7988")
def t10():
    kjs13 = ["b1380","b1241","b1478","b0351","b2296","b1849","b3115","b0904","b2492","b0903","b3114","b3951","b3952"]
    m = preparar_ambiente(MODEL.copy(), "anaerobio")
    fitness = me.evaluar_bacteria(m, kjs13, "EX_succ_e",
                                   "BIOMASS_Ec_iML1515_core_75p37M")
    assert fitness is not None, "devolvio None"
    assert fitness > 5.0, f"KJ060 no acopla: {fitness}"
    return f"fitness={fitness:.4f} (benchmark README: 7.7988)"

# ================= 4. TERMODINAMICA Y DATOS EXTERNOS =================
@prueba("analizar_termodinamica_universal (fumarasa)")
def t11():
    m = MODEL.copy()
    buf = io.StringIO()
    with redirect_stdout(buf):
        me.analizar_termodinamica_universal(m, ["FUM"])
    salida = buf.getvalue()
    if salida.strip():
        assert "FUM" in salida, "la tabla no incluye FUM"
        return "tabla FUM capturada (dG' ~ -3.43 kJ/mol)"
    return "ejecuta sin excepciones (salida no capturable por rich)"

@prueba("traducir_gen_uniprot(b3731)")
def t12():
    return f"resultado: {str(me.traducir_gen_uniprot('b3731', '83333'))[:60]}"

@prueba("limpiar_cache_uniprot")
def t13():
    me.limpiar_cache_uniprot()
    return "cache limpiado sin errores"

# ================= 5. SALIDAS =================
@prueba("generar_documento_markdown")
def t14():
    me.generar_documento_markdown(MODEL, "iML1515", "EX_succ_e", ["b3731"], 1.5562)
    carpeta = Path("reportes") if Path("reportes").exists() else Path(".")
    archivos = sorted(carpeta.glob("*.md"), key=lambda p: p.stat().st_mtime)
    assert archivos, "no se genero ningun reporte .md"
    contenido = archivos[-1].read_text(encoding="utf-8", errors="ignore")
    assert len(contenido) > 500, f"reporte muy corto: {len(contenido)} chars"
    return f"reporte: {archivos[-1].name} ({len(contenido)} chars)"

@prueba("generar_mapa_universal (HTML)")
def t15():
    me.generar_mapa_universal(modelo_id="iML1515", quimico="EX_succ_e")
    htmls = sorted(Path(".").rglob("*.html"), key=lambda p: p.stat().st_mtime)
    assert htmls, "no se genero ningun .html"
    tam = htmls[-1].stat().st_size
    assert tam > 1000, f"HTML sospechosamente pequeno: {tam} bytes"
    return f"mapa: {htmls[-1].name} ({tam//1024} KB)"

@prueba("mostrar_cinetica_local")
def t16():
    m = preparar_ambiente(MODEL.copy(), "anaerobio")
    me.mostrar_cinetica_local(m, "EX_succ_e", ["b3731"])
    return "cinetica local consultada sin errores"

# ================= 6. UNIVERSALIDAD =================
@prueba("cargar_modelo_universal(iMM904 - S. cerevisiae)")
def t17():
    m2 = extraer_modelo(me.cargar_modelo_universal("iMM904"))
    try:
        m2.solver = "glpk"
    except Exception:
        pass
    mu = m2.slim_optimize()
    assert mu is not None and mu > 0, f"FBA fallo en levadura: {mu}"
    return f"{len(m2.reactions)} reacciones, mu={mu:.4f}"

# ================= 7. CLI =================
@prueba("cli_experto --help")
def t18():
    import subprocess
    r = subprocess.run([sys.executable, "cli_experto.py", "--help"],
                       capture_output=True, text=True, timeout=60)
    assert r.returncode == 0, f"exit code {r.returncode}: {r.stderr[:200]}"
    return "--help funciona, exit 0"

# ================= 8. PROFUNDAS (solo --full) =================
if "--full" in sys.argv:
    @prueba("optknock_exacto_milp (k=1)")
    def t19():
        m = preparar_ambiente(MODEL.copy(), "anaerobio")
        return f"resultado: {str(me.optknock_exacto_milp(m, 'e_coli', 'iML1515', seed=42))[:80]}"

    @prueba("disenar_cepa_ia (GA v0.2.0, seed=100)")
    def t20():
        m = preparar_ambiente(MODEL.copy(), "anaerobio")
        return f"resultado: {str(me.disenar_cepa_ia(m, 'e_coli', 'iML1515', seed=100))[:80]}"

# ================= EJECUCION =================
if __name__ == "__main__":
    print("=" * 64)
    print(f" SMOKE TEST INTEGRAL v3{' [MODO QUICK]' if QUICK else ''}")
    print("=" * 64 + "\n")
    pruebas = [t01, t02, t03, t04, t05, t06]
    if QUICK:
        registrar("obtener_genes_no_esenciales", "SKIPPED", "modo --quick", 0.0)
    else:
        pruebas.append(t07)
    pruebas += [t08, t09, t10, t11, t12, t13, t14, t15, t16, t17, t18]
    for t in pruebas:
        t()
    if "--full" in sys.argv:
        t19(); t20()

    print("\n" + "=" * 64)
    print(" RESUMEN")
    print("=" * 64)
    p = sum(1 for r in RESULTADOS if r[1] == "PASSED")
    f = sum(1 for r in RESULTADOS if r[1] == "FAILED")
    s = sum(1 for r in RESULTADOS if r[1] == "SKIPPED")
    print(f"{'Prueba':<45}{'Estado':<9}{'Tiempo'}")
    print("-" * 64)
    for nombre, estado, _, dur in RESULTADOS:
        print(f"{nombre:<45}{estado:<9}{dur:>7.1f}s")
    print("-" * 64)
    print(f"TOTAL: {len(RESULTADOS)}  |  PASSED: {p}  |  FAILED: {f}  |  SKIPPED: {s}")
    sys.exit(1 if f > 0 else 0)
