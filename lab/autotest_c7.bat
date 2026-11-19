@echo off
rem C7 autotest: motor A on 6 digits, any order, 5 windows; compares the output hash with Linux.
rem Writes only inside this folder, low priority. About 1-3 minutes and 2 GB of RAM.
cd /d %~dp0
set OMP_NUM_THREADS=2
echo Autotest C7...
start "" /BELOWNORMAL /WAIT /B cmd /c "%~dp0motora.exe multi 6 0 res_6_0.txt w_ident6.txt ident6_out.txt 2> ident6.log"
certutil -hashfile ident6_out.txt SHA256 | findstr /v ":" > ident6_hash.txt
set /p H=<ident6_hash.txt
set /p E=<esperado_ident6.txt
set H=%H: =%
if /i "%H%"=="%E%" (echo AUTOTEST OK: salida identica a Linux) else (echo AUTOTEST FALLA: %H% & echo esperado %E%)
