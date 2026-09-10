@echo off
setlocal
call "C:\Zsq\Programs\Visual Studio\Visual Studio2022\Community\VC\Auxiliary\Build\vcvarsall.bat" x64 >nul 2>&1
if errorlevel 1 (
    echo [ERROR] vcvarsall failed
    exit /b 1
)
set "PROJECT_ROOT=%~dp0"
cmake -S "%PROJECT_ROOT%." -B "%PROJECT_ROOT%build-tools" -G "NMake Makefiles" ^
    -DOCC_ROOT="%OCC_ROOT%"
if errorlevel 1 exit /b 1
cmake --build "%PROJECT_ROOT%build-tools"
exit /b %errorlevel%
