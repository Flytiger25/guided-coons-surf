# SGK 替换 OCC 迁移日志

> 目标：把 guided-coons-surf 项目对 OpenCASCADE (OCC) 的全部依赖替换为国产几何内核 SGK。
> 工作方式：逐模块推进、每步验证、每步 git 提交可回溯。
> 分支：`sgk-migration`（基线：main @ `5a90f4f`）

---

### 第 0 步：建立 git 基线并初始化日志

- 时间：2026-09-06 18:14
- 改动文件：`.gitignore`（新增忽略 `run.log`）；提交当前 OCC 版本改动（`.vscode/*`、`CMakeLists.txt`、`README.md`、`src/*.cpp`、`src/main.cpp`、`run.ps1`）
- 验证方式：`git status` 干净、`git log` 有基线提交、新分支 `sgk-migration` 已从 main 切出
- 结果：✅ 成功
- commit：`5a90f4f`（main 基线）

---

### 第 1 步：SGK 冒烟测试

- 时间：2026-09-07 11:30
- 改动文件：新增 `test/sgk_smoke.cpp`、`test/CMakeLists.txt`、`test/build_smoke.bat`；`.gitignore` 加 `test/build-smoke/`
- 验证方式：编译 + 运行 `sgk_smoke.exe`，逐一验证 init/许可证、B 样条曲线/曲面构造求值、点投影、STEP 写读
- 结果：✅ 成功（全部 API 通过）

**关键结论（直接影响后续迁移）：**
1. **编译必须用 `/MD` + Release**（NMake 生成器 + `vcvarsall` 环境，绕开本机崩溃的 MSBuild）。若用 `/MDd`（Debug），`_ITERATOR_DEBUG_LEVEL` 与 SGK 的 Release DLL 不匹配，跨 DLL 传 `std::vector` 会 `vector too long` 异常或段错误。
2. **点投影用 `BSplineSurface::CalcNearestPnt(pnt, uvParam)`** 返回最近点 + (u,v)，等价 OCC `GeomAPI_ProjectPointOnSurf`。⚠️ `GeomProject::PntSrfProject` **声明了但未实现**（抛"待完善"异常），不可用。
3. STEP 写需要一个**带边界 Loop 的完整 Face**（`TopoBuilder::MakeEdge`×4 → `MakeCoedge` → `MakeLoop` → `FaceAddLoop`）；直接 `MakeFace(msrf)` 无 Loop 会报 "No loop in face"。
4. 许可证有效（`sggk::init()` 成功，STEP 读写正常，说明 `[DATAEXCHANGE]` 授权可用）。

- commit：`d79dee5`

---

### 第 2 步：数据转换 .brep → .step

- 时间：2026-09-07 11:35
- 改动文件：新增 `data/input/1_boundary.step`、`1_internal.step`；新增 `tools/brep2step.cpp` + `tools/CMakeLists.txt`（OCC 一次性转换工具）；新增 `test/read_step.cpp`（SGK 读取验证）
- 验证方式：`brep2step` 把两个 `.brep` 转 `.step`；`read_step` 用 SGK 读回
- 结果：✅ 成功

**验证数据（SGK 读 .step 结果）：**
- `1_boundary.step`：4 条边（全部 BSplineCurve3D，degree=3，控制点 14/13/21/17）
- `1_internal.step`：26 条边（全部 BSplineCurve3D，degree=3，控制点 52~174）

**结论：** SGK 能正确读取 OCC 导出的 `.step`，边数/曲线类型与 `.brep` 一致，`Edge::GeomCurve()->ToBSpline()` 可提取 B 样条曲线（对应 OCC 的 `GeomConvert::CurveToBSplineCurve`）。

- commit：`4e75c71`

---

### 第 3 步：迁移 KnotUpdate 模块

- 时间：2026-09-07 11:45
- 改动文件：`include/KnotUpdate.h`、`src/KnotUpdate.cpp`（OCC → SGK 类型替换）
- 改动要点：`Handle(Geom_BSplineCurve)&` → `const sggk::BSplineCurve3DPtr&`；`bspline->Value(u)` → `bspline->CalcPoint(u)`；`gp_Pnt::Distance` → `Point3D::DistanceTo`；`Standard_Real/Integer/Boolean` → `double/int/bool`；`Precision::Angular()` → `1e-12`
- 验证方式：新增 `test/test_knotupdate.cpp`，单独编译 `src/KnotUpdate.cpp` + 运行；grep 确认无 OCC include
- 结果：✅ 成功（运行输出 `newKnot=0.5 maxError=0.00375 sequences=9`，节点插入逻辑正确；无 OCC 残留）

- commit：待提交

---

