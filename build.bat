@echo off
setlocal
call "C:\Zsq\Programs\Visual Studio\Visual Studio2022\Community\VC\Auxiliary\Build\vcvarsall.bat" x64 >nul 2>&1
if errorlevel 1 (
    echo [ERROR] vcvarsall failed
    exit /b 1
)
cmake -S "C:\Zsq\Projects\conns-surf\occ_SGK\guided-coons-surf" -B "C:\Zsq\Projects\conns-surf\occ_SGK\guided-coons-surf\build-sgk" -G "NMake Makefiles"
if errorlevel 1 exit /b 1
cmake --build "C:\Zsq\Projects\conns-surf\occ_SGK\guided-coons-surf\build-sgk"
exit /b %errorlevel%
