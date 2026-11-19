@echo off
rem C7: motor A, 9 digits, any order, windows for pi, gamma, delta (w_c7_9_0.txt). Low priority.
rem Writes c7_cand.txt (candidate paths) and c7.log (one line per window). Only inside this folder.
cd /d %~dp0
if "%1"=="" (set OMP_NUM_THREADS=%NUMBER_OF_PROCESSORS%) else (set OMP_NUM_THREADS=%1)
echo Lanzado %date% %time% con %OMP_NUM_THREADS% hilos > c7_inicio.txt
start "C7" /BELOWNORMAL /MIN cmd /c "%~dp0motora.exe multi 9 0 res_9_0.txt w_c7_9_0.txt c7_cand.txt 2> c7.log & echo FIN %%date%% %%time%% > c7_fin.txt"
echo C7 lanzado en segundo plano con prioridad baja. Estado: type c7.log
