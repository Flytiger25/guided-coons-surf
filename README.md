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
│   └── input/                  # 测试输入数据 (.brep)
├── CMakeLists.txt
└── README.md
```

## 依赖

| 依赖 | 版本 | 用途 |
|------|------|------|
| [OpenCASCADE](https://dev.opencascade.org/) | 7.9 | 几何内核（B 样条曲线/曲面、STEP 读写） |
| [Eigen](https://eigen.tuxfamily.org/) | 3.4 | 线性代数（稀疏矩阵、KKT 系统求解） |
| CMake | ≥ 3.10 | 构建系统 |

## 构建（macOS）

### 安装依赖

```bash
brew install opencascade eigen cmake
```

### 编译

```bash
mkdir build && cd build
cmake ..
make -j$(sysctl -n hw.logicalcpu)
```

> **注意**：`CMakeLists.txt` 中 OpenCASCADE 和 Eigen 路径硬编码为 Homebrew 默认路径。如果版本不同，请根据实际安装路径修改。

### 运行

```bash
cd build
./guided_coons_surf
```

程序读取 `data/input/` 下的 `.brep` 边界线和引导线文件，生成曲面并导出为 STEP 格式。

## 输入/输出

- **输入**：B 样条曲线（`.brep` 或 `.step` 格式）
  - `1_boundary.brep` — 四条边界曲线
  - `1_internal.brep` — 内部引导线
- **输出**：B 样条曲面（`.step` 格式），导出至 `data/` 目录

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
