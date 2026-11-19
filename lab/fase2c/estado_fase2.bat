@echo off
rem Estado de FASE 2c. SOLO LECTURA: porcentaje de cada paso, horas restantes estimadas, memoria libre, avisos y procesos.
cd /d %~dp0
py -3.13 -B "%~dp0fase2.py" estado
