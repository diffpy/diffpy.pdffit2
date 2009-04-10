@echo off
setlocal
set PYTHON="%~dp0..\python.exe"
set PDFFIT2_SCRIPT="%~dp0pdffit2"
%PYTHON% -i %PDFFIT2_SCRIPT% %*
