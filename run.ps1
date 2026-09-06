# run.ps1 —— 一键编译 + 运行 guided_coons_surf
# 用法：在 PowerShell 里执行  .\run.ps1
Set-Location $PSScriptRoot

# 1. OpenCASCADE DLL 路径（不设置会报找不到 TKernel.dll）
$env:PATH = "C:\Zsq\Develop\OpenCASCADE-7.7.0-vc14-64\opencascade-7.7.0\win64\vc14\bin;$env:PATH"

# 2. 首次运行 / 改过 CMakeLists.txt 时需要配置，可手动执行：
#    cmake -S . -B build -G "Visual Studio 17 2022" -A x64

# 3. 编译
Write-Host "===== 编译 =====" -ForegroundColor Cyan
cmake --build build --config Release
if ($LASTEXITCODE -ne 0) {
    Write-Host "编译失败，请检查上方报错" -ForegroundColor Red
    exit 1
}

# 4. 运行
Write-Host "===== 运行 =====" -ForegroundColor Cyan
.\build\Release\guided_coons_surf.exe

Write-Host "===== 完成 =====" -ForegroundColor Green
