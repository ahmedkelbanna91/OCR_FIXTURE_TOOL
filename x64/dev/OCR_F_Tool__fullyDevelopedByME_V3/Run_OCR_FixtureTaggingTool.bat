@Rem Edited by Banna
@echo off
setlocal EnableDelayedExpansion
echo [93m============================="Created by Banna"===============================[0m
echo [93m========================="OCR_FixtureTaggingTool V3"==========================[0m

cd /D %~dp0
SET output_path=%~dp0output
IF not exist %output_path% MD %output_path%


echo:
cd /D %~dp0

SET inCNT=0
SET outCNT=0
SET doneCNT=0

ECHO       [96m"Cleaning output"[0m
IF EXIST "output\*.stl" DEL /S /Q %~dp0output\*.* >nul
echo:

cd /D %~dp0

set /p CaseID="-     [93mWhat is the Case ID?[0m "
set /p upper="-     [93mHow many UPPER?[0m "
set /p Upassive="-     [93mHow many UPPER PASSIVE?[0m "
set /p Uretainer="-     [93mHow many UPPEER RETAINER?[0m "
set /p Utemplate="-     [93mHow many UPPER TEMPLATE?[0m "

set /p lower="-     [93mHow many LOWER?[0m "
set /p Lpassive="-     [93mHow many LOWER PASSIVE?[0m "
set /p Lretainer="-     [93mHow many LOWER RETAINER?[0m "
set /p Ltemplate="-     [93mHow many LOWER TEMPLATE?[0m "

if not defined upper set "upper=0"
if not defined Upassive set "Upassive=0"
if not defined Uretainer set "Uretainer=0"
if not defined Utemplate set "Utemplate=0"

if not defined lower set "lower=0"
if not defined Lpassive set "Lpassive=0"
if not defined Lretainer set "Lretainer=0"
if not defined Ltemplate set "Ltemplate=0"


set /a inCNT=%upper% + %Uretainer% + %Utemplate% + %Upassive% + %lower% + %Lretainer% + %Ltemplate% + %Lpassive%
echo:
echo       [92m%inCNT%[0m  [93mFixtures will be created with OCR Tag[0m
echo:
echo       USERNAME:  %USERDOMAIN%\%USERNAME%
set start_Time=%TIME%
echo       [93mStart Time: %start_Time%[0m

echo:
echo [93m================================"Adding Tag"==================================[0m
echo:

cd /D %~dp0DO_NOT_OPEN_OR_DELETE

for /L %%i in (1,1,%upper%) do call :ProcessModel %%i UN "UPPER"
for /L %%i in (1,1,%Upassive%) do call :ProcessModel %%i UP "UPPER PASSIVE"
for /L %%i in (1,1,%Uretainer%) do call :ProcessModel %%i UR "UPPER RETAINER"
for /L %%i in (1,1,%Utemplate%) do call :ProcessModel %%i UT "UPPER TEMPLATE"

for /L %%i in (1,1,%lower%) do call :ProcessModel %%i LN "LOWER"
for /L %%i in (1,1,%Lpassive%) do call :ProcessModel %%i LP "LOWER PASSIVE"
for /L %%i in (1,1,%Lretainer%) do call :ProcessModel %%i LR "LOWER RETAINER"
for /L %%i in (1,1,%Ltemplate%) do call :ProcessModel %%i LT "LOWER TEMPLATE"

cd /D %~dp0

echo:
echo [93m================================="Finished"===================================[0m
echo [93m=================================="REPORT"====================================[0m

FOR /F "delims=" %%i IN ('dir /b /s %output_path%\*.stl') DO SET /a outCNT+=1
IF %outCNT% GTR 0 (
    echo:
    echo       [92m%outCNT%/%inCNT%[0m  OCR Tag combined with fixture [92msuccessfully[0m.
    echo:
)

echo       [92m%outCNT%[0m  STL in "output" - Fixture with OCR Tag


echo:
set end_time=%time%
set /a hours=%end_time:~0,2%-%start_time:~0,2%
set /a minutes=%end_time:~3,2%-%start_time:~3,2%
set /a seconds=%end_time:~6,2%-%start_time:~6,2%
set /a total_seconds=%hours%*3600 + %minutes%*60 + %seconds%
set /a minutes=%total_seconds% / 60
set /a seconds=%total_seconds% %% 60
echo:
echo       USERNAME:  %USERDOMAIN%\%USERNAME%
echo:
echo       [93mElapsed:  %minutes% minutes %seconds% seconds[0m
echo:
echo       Start  :  %start_Time%
echo       Finish :  %time%
echo:
pause


:ProcessModel
set "num=%1"
set "LABEL=%2"
set "type=%3"
if %1 lss 10 (
    set "num=0%1"
)
set "ID=!CaseID!!LABEL!!num!"
CALL OCR_FIXTURE_TOOL_V3.exe -O "%output_path%\!ID!_F.stl" -N "!ID!" -D "-1.0"
IF exist "%output_path%\!ID!_F.stl" (
    SET /a doneCNT+=1
    echo       ++  !doneCNT!/%inCNT% - !ID!_F.stl exported %type%    OCR Tag added successfully.
) ELSE (
    echo       --  !doneCNT!/%inCNT% - !ID!_F.stl export %type%    OCR Tag adding failed.
)
echo:
goto :eof