@echo off
setlocal
set PYTHONINSPECT=1
set PYTHON="%~dp0..\python.exe"
set PDFFIT2_SCRIPT="%~dp0pdffit2"
%PYTHON% %PDFFIT2_SCRIPT% %*
