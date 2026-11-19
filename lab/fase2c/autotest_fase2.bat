@echo off
rem Autotest de FASE 2c (5-15 min, prioridad baja). Debe decir AUTOTEST OK antes de lanzar nada.
rem Solo escribe dentro de esta carpeta (work\autotest) y copia ..\c7\res_9_0.txt a py\cert tras comprobar su hash.
cd /d %~dp0
echo Autotest de FASE 2c (5-15 min)...
start "" /BELOWNORMAL /WAIT /B cmd /c "py -3.13 -B "%~dp0fase2.py" autotest > "%~dp0autotest_fase2.log" 2>&1"
type autotest_fase2.log
