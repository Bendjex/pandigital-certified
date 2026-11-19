@echo off
rem Estado de C7. SOLO LECTURA.
cd /d %~dp0
echo === C7  %date% %time% ===
if exist c7_inicio.txt type c7_inicio.txt
if exist c7.log (powershell -NoProfile -Command "Get-Content c7.log | Select-String -Pattern 'progress','candidates','tables size 5' | Select-Object -Last 4 | ForEach-Object { $_.Line }")
powershell -NoProfile -Command "Get-Process motora -ErrorAction SilentlyContinue | Format-Table Id,PriorityClass,@{n='CPU_h';e={[math]::Round($_.CPU/3600,1)}},@{n='RAM_GB';e={[math]::Round($_.WorkingSet64/1GB,2)}},@{n='Hilos';e={$_.Threads.Count}} -AutoSize | Out-String"
if exist c7_fin.txt (echo TERMINADO: & type c7_fin.txt)
