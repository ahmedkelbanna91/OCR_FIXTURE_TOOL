@Rem Edited by Banna

@ECHO OFF
echo:
echo [93m================================"CLEAN"=====================================[0m
echo:

pause

    ECHO [96m"Cleaning input"[0m
    DEL /S /Q %~dp0input\*.*
    ECHO [96m"Cleaning output"[0m
    DEL /S /Q %~dp0output\*.*
    ECHO [96m"Cleaning Temp"[0m
    DEL /Q %~dp0\DO_NOT_OPEN_OR_DELETE\model\input\*.*
    DEL /Q %~dp0\DO_NOT_OPEN_OR_DELETE\model\output\*.*

timeout 2