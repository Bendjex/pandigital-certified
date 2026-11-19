@echo off
rem FASE 2c: lanza el orquestador en segundo plano (prioridad baja). Si ya habia empezado, reanuda.
rem No busca ventanas por titulo: el orquestador se identifica por su PID y su hora de creacion (fase2_pid.txt).
cd /d %~dp0
if exist fase2_fin.txt (echo FASE 2c ya terminada: ver resultados_fase2c.zip & exit /b 1)
py -3.13 -B "%~dp0fase2.py" vivo && (echo FASE 2c ya esta en marcha. Estado: estado_fase2.bat & exit /b 1)
if not exist work\autotest\AUTOTEST_OK.txt (echo Antes hay que pasar autotest_fase2.bat con AUTOTEST OK. No se lanza nada. & exit /b 1)
echo Lanzado %date% %time% >> fase2_inicio.txt
start "orquestador fase2c" /BELOWNORMAL /MIN cmd /c "py -3.13 -B "%~dp0fase2.py" run >> "%~dp0fase2_consola.log" 2>&1"
echo FASE 2c lanzada en segundo plano. Estado: estado_fase2.bat   Log: fase2.log
