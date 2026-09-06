# GuidedCoonsSurf

带引导线的扩展 Coons 曲面生成算法。

## 概述

在经典 Coons 曲面（四条边界线插值）的基础上，引入内部**引导线（guide curves）**作为约束条件，通过迭代优化使曲面精确包络引导线，同时保持曲面的光顺性。

### 核心算法

1. **Coons 曲面构造** — 由四条 B 样条边界线通过双线性混合生成初始张量积曲面
2. **引导线投影** — 将引导线采样点投影到当前曲面，计算偏移向量
3. **偏移曲面拟合** — 以最小二乘 + 双向等参线光顺能量拟合偏移曲面
4. **曲面叠加** — 在控制点层面叠加偏移量，修正曲面形状
5. **自适应节点细化** — 在误差超差区域插入节点，增加局部自由度
6. **迭代收敛** — 重复步骤 2-5，直至所有采样点偏差满足容差要求

### 曲线光顺（CurveFair）

独立的 B 样条曲线光顺模块，通过**弧长参数化三阶导数能量极小化**消除曲线上的不必要波动：

$$E = \int \left\| \frac{d^3 C(s)}{ds^3} \right\|^2 ds$$

- 构造弧长参数映射函数 $f(s) = t$，使得复合函数 $C(f(s))$ 满足弧长参数化
- 通过链式法则将对弧长的三阶导数转换为对曲线参数 $t$ 的导数
- 30 点 Gauss-Legendre 数值积分构建稀疏能量矩阵
- 带端点切向约束的 KKT 系统求解，在光顺度和形状保真度之间平衡

### 曲面光顺

采用**双向等参线光顺**策略：对 U/V 两族等参线分别计算曲线光顺能量，聚合为曲面光顺矩阵，避免了薄板能量过度扁平化的问题。

## 目录结构

```
.
├── include/                    # 头文件
│   ├── CurveFair.h
│   ├── GuidedCoonsSurfGenerator.h
│   └── KnotUpdate.h
├── src/                        # 源文件
│   ├── main.cpp
│   ├── CurveFair.cpp
│   ├── GuidedCoonsSurfGenerator.cpp
│   └── KnotUpdate.cpp
├── data/
│   ├── input/                  # 测试输入数据 (.brep)
│   ├── output/                 # 主输出 (STEP)
│   └── coons/                  # 迭代调试导出 (STEP)
├── CMakeLists.txt
└── README.md
```

## 依赖

| 依赖 | 版本 | 用途 |
|------|------|------|
| [OpenCASCADE](https://dev.opencascade.org/) | 7.7（Windows）/ brew 最新（macOS） | 几何内核（B 样条曲线/曲面、STEP 读写） |
| [Eigen](https://eigen.tuxfamily.org/) | 3.4 | 线性代数（稀疏矩阵、KKT 系统求解），header-only |
| CMake | ≥ 3.10 | 构建系统 |

> **路径说明**：`CMakeLists.txt` 中 OpenCASCADE 和 Eigen 路径为 Windows 硬编码路径，请按本机实际安装位置修改。数据目录通过编译宏 `GUIDED_COONS_DATA_DIR` 注入，与运行时工作目录无关。

## 构建（Windows）

### 环境要求

- Visual Studio 2019 / 2022（含 C++ 桌面开发工作负载）
- [CMake](https://cmake.org/download/) ≥ 3.10
- OpenCASCADE 7.7.0 Windows 版（已解压，如 `C:/Zsq/Develop/OpenCASCADE-7.7.0-vc14-64`）
- Eigen 3.4.0（header-only，解压即可，如 `C:/Zsq/Develop/eigen-3.4.0`）

### 配置与编译

```powershell
# 在项目根目录执行（需将 cmake 加入 PATH，或用 Visual Studio 开发者终端）
cmake -S . -B build -G "Visual Studio 17 2022" -A x64
cmake --build build --config Release
```

也可直接用 Visual Studio 打开 `build/guided_coons_surf.sln` 编译。

### 运行

```powershell
.\build\Release\guided_coons_surf.exe
```

程序读取 `data/input/` 下的 `.brep` 边界线和引导线文件，生成曲面并导出为 STEP 格式（`data/output/`）。

## 构建（macOS）

### 安装依赖

```bash
brew install opencascade eigen cmake
```

### 编译与运行

```bash
mkdir build && cd build
cmake ..
make -j$(sysctl -n hw.logicalcpu)
./guided_coons_surf
```

> **注意**：macOS 下需将 `CMakeLists.txt` 中 OpenCASCADE / Eigen 路径改为 Homebrew 实际安装路径。

## 输入/输出

- **输入**：B 样条曲线（`.brep` 或 `.step` 格式）
  - `1_boundary.brep` — 四条边界曲线
  - `1_internal.brep` — 内部引导线
- **输出**：B 样条曲面（`.step` 格式），导出至 `data/output/`；迭代中间结果导出至 `data/coons/`

## 模块说明

### GuidedCoonsSurfGenerator

带引导线的 Coons 曲面生成器，核心入口：
- `Perform()` — 执行算法主流程
- `ConstructCoonsSurf()` — 构造初始 Coons 曲面
- `ConstructSurfWithGuideCrvs()` — 引导线迭代修正

### CurveFair

B 样条曲线光顺工具，可独立使用：
- `Perform()` — 执行光顺
- `ComputeEnergyMatrix()` / `ComputeSparseEnergyMatrix()` — 计算三阶导数能量矩阵
- `GetTempFairCurveWithTangentConstraint()` — 带端点切向约束的优化求解

### KnotUpdate

自适应节点更新，支持多种策略：
- `PARAM_BASED_BY_INTERVAL_ERROR` — 按区间误差分布插入节点
- `MID_KNOT_BY_SINGLE_ERROR` — 在最大误差处插入节点
- `UNIFORM_UP_DATE` — 均匀加细

## License

MIT
