@echo off
rem Estado de FASE 2. SOLO LECTURA: ultimo bloque de progreso (porcentajes por tarea), avisos y procesos.
cd /d %~dp0
py -3.13 -B "%~dp0fase2.py" estado
