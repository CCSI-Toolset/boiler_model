@echo off
REM Assemble the zip file

set VERSION=2016.02.00
set WORK_DIR=.jenkins_work
set PROJECT_FULL_NAME="Boiler Model for Windows"
set PROJECT_NAME=BoilerModel_Windows_%VERSION%
set CCSI_COMMON=.ccsi_common

mkdir %WORK_DIR%\%PROJECT_NAME% || goto :error
cd %WORK_DIR%\%PROJECT_NAME% || goto :error
copy ..\..\Setup_Boiler_Model\Setup_Boiler_Model\Express\SingleImage\DiskImages\DISK1\setup.exe || goto :error
copy ..\..\%CCSI_COMMON%\LEGAL || goto :error

powershell -Command "(Get-Content ..\..\%CCSI_COMMON%\CCSI_TE_LICENSE.txt) | ForEach-Object { $_ -replace  '\[SOFTWARE NAME & VERSION\]', '%PROJECT_FULL_NAME% v.%VERSION%' } | Set-Content CCSI_TE_LICENSE_%PROJECT_NAME%.txt"

cd .. || goto :error
7za a -bd CCSI_%PROJECT_NAME%.zip %PROJECT_NAME% || goto :error
cd .. || goto :error
move %WORK_DIR%\CCSI_%PROJECT_NAME%.zip . || goto :error
rmdir /Q /S %WORK_DIR% || goto :error
goto :EOF

:error
echo Failed with error #%errorlevel%.
exit /b %errorlevel%
