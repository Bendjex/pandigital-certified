@echo off
rem Para FASE 2c (SOLO si Claude o Daniel lo piden): para por PID (y hora de creacion) el orquestador y los procesos que lanzo.
rem No toca ninguna otra ventana ni programa. Lo hecho se conserva: lanzar_fase2.bat reanuda.
cd /d %~dp0
py -3.13 -B "%~dp0fase2.py" parar
