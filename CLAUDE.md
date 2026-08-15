# Hadwiger–Nelson 下界 6 — 会话必读

找有限点集 V ⊂ R²，其单位距离图不可 5 染色 ⇒ χ(R²) ≥ 6。截至目前已知界仍是 5 ≤ χ ≤ 7，**本仓库尚未声称获得下界 6**。Python + SymPy 精确代数坐标 + 多 SAT 求解器（CaDiCaL/Kissat 等外装）+ NetworkX。起点是 Heule 529 点图（`data/heule529/`）。

## 红线

1. **数值只能用于预筛**；边的判定必须走 SymPy 精确代数（`exact_geometry.py`）。曾发生 `np.int64` 溢出 + `evalf` 返回伪造坐标的严重 bug（复盘在 research_log），这是本项目最大的教训。
2. `Point.approximate` 必须用 `maxn=100000` 保护精度。
3. solver 的 UNSAT 必须配 DRAT/LRAT 证书才算数。
4. 独立验证器（`*_verifier.py`）刻意不 import 搜索器，保持隔离。

## 现状与未决（2026-08-13 盘点，接手先搞清这三件）

1. research_log 结尾记录一批增强列边子图 CNF（14.5 万变量）在多个求解器上 `UNKNOWN/RUNNING`——**这些进程的最终归宿没有记录**【待确认：跑完了没有、结果在哪】。
2. 2026-08-10 的 `artifacts/cegis_38029_projective_pairs_color3/` 工件在 README 和 research_log 里**都没有对应记述**——最新一轮实验未入日志。
3. `artifacts/`（22 个 cegis_* 目录，共 7.5GB）**没有索引**：哪个是当前 checkpoint、哪个是已废弃审计件（README 提到 1428 候选/2489 点工件已废弃），需要一份清单。

## 怎么跑

README 无「复现」节、无 Makefile——流水线入口按阶段：`cegis.py` → `closure_cegis.py` → `d4_*_cegis.py` 系列；测试 `python -m pytest tests/`（61 个）。外部求解器安装方式未记录【待补】。**本项目无 git**（third_party 下的 .git 是第三方自带）。
