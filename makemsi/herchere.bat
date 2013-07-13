@ECHO OFF
set _P=%~1%
set PATH=%_P%;%PATH%
%comspec% /K setherc.bat