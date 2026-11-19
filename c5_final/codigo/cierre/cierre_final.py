"""(17/11) the final closure of n = 9 after the last results of the cloud are imported (importar_nube.bat): the
laboratory writes cierre_9.json and resultados_cierre.zip when ITS part ends, perhaps before all the cloud sessions are
in.  This recomputes, from all the records (the lab's parts and the imported c5_9_parte_nubeXX.jsonl), the closure of
n = 9 (the same function closure() of cierre.py, untouched), writes work/cierre_9.json (the earlier one is kept as
cierre_9.json.antes_<date>) and packs resultados_cierre.zip again.  It only reads records; it computes nothing.
Run with the closing stopped or finished.  usage: python cierre_final.py   (cierre_final.bat)"""
import os, sys, json, time, shutil
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cierre as C

if C.running():
    print("el cierre esta en marcha: ejecutar cierre_final cuando haya terminado (cierre_fin.txt) o este parado")
    sys.exit(1)
out = 0
for cell in C.CFG["cells"]:
    if cell["n"] != 9:
        continue
    cl = C.closure(cell)
    fn = os.path.join(C.WORK, "cierre_9.json")
    if os.path.exists(fn):
        shutil.copyfile(fn, fn + ".antes_" + time.strftime("%Y%m%d_%H%M%S"))
    with open(fn + ".tmp", "w", encoding="utf-8") as f:
        json.dump(cl, f, indent=1)
    os.replace(fn + ".tmp", fn)
    txt = C.closure_text(cl)
    C.log("cierre_final (tras importar la nube): n = 9\n" + txt)          # (log prints it too)
    allc = all(d["closed"] and d.get("closed_c11", True) for d in cl["consts"].values())
    print(f"C5 n = 9 en cualquier orden: {'CERRADA (C7 y C11, las tres constantes)' if allc else 'AUN NO CERRADA'}")
    out = 0 if allc else 2
C.pack()
print("resultados_cierre.zip rehecho")
sys.exit(out)
