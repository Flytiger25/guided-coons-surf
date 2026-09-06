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
