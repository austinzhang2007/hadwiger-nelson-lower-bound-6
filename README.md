# Hadwiger–Nelson 下界 6：精确几何与 SAT 搜索

本项目寻找有限点集 \(V\subset\mathbb R^2\)，其**完整单位距离图**

\[
E(V)=\{\{u,v\}\subset V:\|u-v\|^2=1\}
\]

不可 5 染色。只有同时通过精确几何验证和可独立检查的 SAT
不可满足证明，才构成“平面色数至少为 6”的有限证书。数值近似、
求解器单独报告 `UNSAT`、或只验证给定边表都不是数学证明。

截至 2026-08-25，已知界仍是
\(5\leq\chi(\mathbb R^2)\leq 7\)。本仓库尚未声称获得下界 6。

2026-08-25 更新：49,820 点、462,566 条已认证列边图已由相位归一化
Kissat walk 找到并独立验证合法 5-染色。其后缺失色 4 全投影扫描得到
3,561 个精确强制点和 617 条精确点对；施加涉及的 977 点后得到 50,797
点、467,091 边图，仍有独立验证的 5-染色。当前正在补扫新增交互边。
详见 `research_log_20260825.md`。

## 当前实验结论

固定种子 `20260725` 已完成两阶段 CEGIS。首轮在第 67 个新反例处
停滞；第二阶段允许任意大小的候选点联合构形：

- 从 Heule 529 点的 94148 对可交单位圆中心生成 187846 个交点记录；
- 数值去重得 151971 个外部交点，1314 个通过“数值度至少 3”预筛；
- 532 个度至少 3 的候选恢复为精确代数坐标并逐邻居符号认证；
- 候选间另有 707 条严格单位距离边；
- 从 128 个规范 5 染色开始，连续得到 67 个新的 5 染色反例；
- 度 3 联合构形从模型 194 继续，再得到 67 个新反例；首个新核是一个
  三点 V 型列表阻塞器；
- 最终全部 532 个候选可同时加入：完整图有 1061 点、5397 边和显式
  5 染色，色类大小 `[253,252,218,186,152]`。

独立几何验证器已确认 1061 点图无重复点、恰有 5397 条单位边且边表
无缺失/伪边；显式 5 染色在全部边上冲突数为 0。因此当前列举的
度至少 3 单位圆交点池已被严格排除。数值预筛并不是对所有理论圆交点
或新代数域的完备枚举，这仍不是下界 6 证明。

第三阶段以这 1061 点作圆心再闭包一层。一次独立全图审计发现旧的
10 位数值去重键会合并坐标极接近但精确不同的交点，并使保留候选漏掉
154 条基图—候选单位边。修复为“14 位预筛＋精确点去重＋完整二分边
重认证”后：

- 320457 对可交圆心，639720 个交点记录；
- 548898 个数值外部交点，3714 个数值度至少 3 种子；
- 1449 个精确候选，候选间 2452 条精确单位边；
- 完整闭包图为 2510 点、13578 边；
- 独立重建缺边 0、伪边 0；
- 显式 5 染色色类 `[578,604,503,450,375]`，全部边冲突 0。

所以第二层闭包也已得到精确失败证书。旧的 1428 候选/2489 点工件仅保留
为数值去重缺陷的审计记录，不作为完整图证书。

第四阶段离开 Moser ring，使用四次非交换扩张中的精确平移
Heule-529 副本。主候选池有 24 个平移向量；CEGIS 连续选入 11 个
列表着色阻塞副本后，剩余候选全部可延拓。随后已对 55 对副本完成
代数基分离：三对有跨副本单位边，数目为 648、655、651；其余 52 对
严格无边。独立 `copy_interaction_verifier.py` 从作者模板重建全部坐标和
边，得到**完整单位距离图** 8329 点、53332 边，缺边 0、伪边 0。
CaDiCaL 找到并独立验证了该完整图的 5-染色，色类大小
`[1534,1625,1684,1806,1680]`。因此此前的 `UNKNOWN/TIMEOUT` 已由
精确 SAT 反例解决，不能作为下界 6 证据。

截至 2026-08-04，在同一非交换 32 维 D4 代数中继续运行阻塞 CEGIS：

- 保留全部候选—基图及候选间单位边的完整版本达到 10229 点、71130
  边，仍有严格验证的 5-染色；1900 个新点中最大基图单位度为 24，
  最大已加入候选交互度为 11；
- 为加速组合搜索，随后只列每个单点阻塞器所需的五条已认证单位边。
  单点池在 11316 点、76565 边时对当前染色重复，图仍可 5 染；
- 升级为“缺一种颜色的四邻居强制点＋同强制色单位边”点对。第一轮
  当前染色给出 1057 个精确强制点对；合并加入 799 点和 4511 边后，
  12115 点、81076 边图仍可 5 染（92.34 秒）。再运行三轮后达到
  14421 点、94796 边、6092 个候选，仍有验证通过的 5-染色；
- 允许所有候选也作圆心的首轮投影闭包扫描了 49174936 对近圆心，
  得到 1366 个第五色数值命中；首批 25 个候选全部精确恢复并加入后，
  14446 点、94921 边图仍可 5 染。继续批量闭包曾达到23051点；独立
  重放确认84614条新增边全部精确，但发现低精度根式消去导致5对重复
  顶点。提高求值保护精度并合并全部邻接约束后，当前干净检查点为
  23046点、137927边，重新求解仍可5染。

后两项是已列明边子图搜索工件；它们的 SAT 模型是严格失败证书，若将来
已列明子图变成 UNSAT，仍须补齐所有单位边、生成并检查 DRAT/LRAT，
才能宣布下界 6。`d4_projective_closure_cegis.py` 进一步允许这些新代数
点本身充当圆心，而不是把候选生成永久限制在原始 8329 点。

2026-08-10 后续实验已经解决上述 29046 点交互图的求解停滞：

- 随机增量边阶梯找到 29046 点、247392 条已认证边的完整
  5-染色，所有列边独立验证通过；
- 经“投影单点闭包→候选间边增广”三轮，依次得到 29546、
  30046、30983 点检查点，每个均有逐边验证的 5-染色；
- 一轮基图强制点对和后续闭包/交互增广达到 31966 点、276361
  条已认证列边，仍 SAT；
- 改为在**全部已认证投影点**上按缺失颜色扫描四色强制点对。
  依次施加缺失色 0、1、2、3 的 602、1478、2758、4531 条精确单位点对，
  得到 42943 点、329638 边检查点；CaDiCaL 用 7.44 秒找到新的合法
  5-染色；
- 缺失色 4 的全投影精确扫描得到 13765 个强制点、6758 条精确
  单位点对，涉及 6877 个候选。全部施加后为 49820 点、363904 条
  已认证列边，仍有合法 5-染色，色类为
  `[7137,9643,10399,11160,11481]`；
- 对 49820 点的点对交互做内存有界的 KD-tree 分块预筛，找到并
  严格认证 98662 条新单位边，目标子图共 462566 条边。正在用
  更细的增量边阶梯求它的 5-染色；所有被中断的难前缀均只记为
  `UNKNOWN_INTERRUPTED`，不当作 UNSAT。

因此最新严格结论仍是：当前列边图 **SAT**，下界 6 **未证明**，
没有可验证的 5-UNSAT DRAT/LRAT 证书。

最终 195 个反例模型上的最佳单点覆盖 133 个。最佳相邻候选对覆盖
138 个，其中 130 个依赖候选间单位边才联合阻塞。全部数据见
`artifacts/cegis/heule529-unit-circles.json`。

## 坐标与图格式

项目格式 `udg-exact-v1` 是 UTF-8 JSON：

```json
{
  "schema": "udg-exact-v1",
  "name": "example",
  "coordinate_format": "sympy-radical-v1",
  "vertices": [
    {"id": 1, "x": "0", "y": "0"},
    {"id": 2, "x": "1", "y": "0"}
  ],
  "edges": [[1, 2]]
}
```

`sympy-radical-v1` 只接受整数、四则运算、括号和 `sqrt(...)`；
解析器不调用 Python `eval`。表达式表示实代数数并由 SymPy
精确化简。若搜索产生不能用根式简洁表示的坐标，后续格式必须记录：

- 原始整系数最小多项式；
- 唯一指定实根的有理隔离区间；
- 各坐标在该数域基上的有理系数。

小数只能作为搜索缓存，不能进入证书。

作者的 Mathematica `.vtx` 文件由受限转换器读取：
`Sqrt[...]` 转换为 `sqrt(...)`，花括号只用作坐标分隔。

## 有限证书规范

成功证书目录必须包含：

1. `graph.json`：所有顶点的无歧义精确代数坐标及完整边表；
2. `graph.edge`：同一图的 DIMACS edge 文件；
3. `coloring-5.cnf`：每顶点恰取一色、每边端点异色、已记录的颜色
   对称破除子句；
4. `coloring-5.drat` 或 `coloring-5.lrat`：不可 5 染色证明；
5. 求解器、证明生成器和验证器的版本、命令行与 SHA-256；
6. `graph_verifier.py`：不读取搜索缓存，从坐标重新生成所有单位距离边，
   检查与两份图文件完全一致，并调用独立 DRAT/LRAT 检查器。

几何验收同时检查：

- 顶点两两不同；
- 每条声明边的平方距离严格等于 1；
- 每个平方距离严格等于 1 的点对都出现在边表；
- 不以浮点容差决定任何最终边。

## 代码边界

- `exact_geometry.py`：受限代数表达式、精确点、平方距离、圆交点及作者
  坐标导入；
- `coloring_sat.py`：图染色 CNF、颜色先出现顺序对称破除、增量模型枚举、
  DIMACS/证明输出；
- `graph_verifier.py`：独立读取证书并从坐标重建完整图；
- `cegis.py`：数值预筛、精确恢复、着色阻塞覆盖与增广反例循环；
- `configuration_cegis.py`：候选间精确单位边、任意大小列表染色阻塞核
  和第二阶段续跑；
- `closure_cegis.py`：检查点式多层圆交点闭包、全候选池直解和完整图
  精确审计；
- `translated_copy_cegis.py`：非 Moser-ring 二次扩张中的 Heule-529
  平移副本、两锚点跨边定理和多副本 CEGIS；
- `translated_copy_verifier.py`：不导入搜索模块，独立重建平移副本
  已列明边子图、验证代数扩张及显式染色；
- `d4_exact.py`：固定 32 维代数基的有理系数运算、齐次坐标相等和单位
  距离判定，并支持用已有分式代数点继续恢复新圆交点；
- `d4_blocker_cegis.py`、`d4_pair_cegis.py`：精确五色单点阻塞与同强制色
  相邻点对 CEGIS；
- `d4_projective_closure_cegis.py`：把所有已加入候选也作为后续圆心的
  投影坐标闭包 CEGIS；
- `d4_projective_pair_scan.py`、`d4_projective_pair_apply.py`：在整个投影
  D4 检查点上恢复四色强制点、筛选并严格认证同强制色单位
  点对，再将稀疏点对证书施加到图中；
- `incremental_edge_ladder.py`：先加入当前染色已满足的所有目标边，
  再随机加入一条当前冲突边并增量重求，用来获得直接冷启动
  容易停滞的大图 SAT 反例；
- `tests/`：Moser spindle 和作者公开 Heule 529 顶点 5 色单位距离图。

## 可复现数据

公开基图来自 Marijn Heule 的作者仓库
[`marijnheule/CNP-SAT`](https://github.com/marijnheule/CNP-SAT)，本地固定提交
`bb414955a6ef5f49f7df2b245b1e778aa67c068a`。使用：

- `third_party/CNP-SAT/vtx/529.vtx`；
- `third_party/CNP-SAT/edge/529.edge`；
- `third_party/CNP-SAT/cnf/529-4.cnf`；
- `third_party/CNP-SAT/proof/529-4-sbp.drat`。

DRAT 检查器来自作者维护的 `marijnheule/drat-trim`，固定提交
`2e3b2dc0ecf938addbd779d42877b6ed69d9a985`。

运行：

```bash
.venv/bin/python -m pytest -q
.venv/bin/python graph_verifier.py \
  --coordinates third_party/CNP-SAT/vtx/529.vtx \
  --edge-file third_party/CNP-SAT/edge/529.edge
third_party/drat-trim/drat-trim \
  data/heule529/529-4-sbp.cnf \
  third_party/CNP-SAT/proof/529-4-sbp.drat
.venv/bin/python cegis.py --colorings 128 --iterations 100 \
  --seed 20260725 \
  --output artifacts/cegis/heule529-unit-circles.json
.venv/bin/python configuration_cegis.py --iterations 100 --min-degree 3 \
  --output artifacts/cegis/heule529-degree3-configurations.json
.venv/bin/python closure_cegis.py --direct-audit --min-degree 3 \
  --output artifacts/cegis/heule1061-second-closure-recertified.json
.venv/bin/python graph_verifier.py \
  --graph-json \
  artifacts/cegis/heule1061-second-closure-recertified-complete.graph.json
.venv/bin/python translated_copy_cegis.py --primary-d4-cegis \
  --output-graph artifacts/cegis/heule2510-primary-d4-cegis.graph.json \
  --output-report artifacts/cegis/heule2510-primary-d4-cegis.json
.venv/bin/python translated_copy_verifier.py \
  --graph artifacts/cegis/heule2510-primary-d4-cegis.graph.json \
  --report artifacts/cegis/heule2510-primary-d4-cegis.json \
  --base-graph \
    artifacts/cegis/heule1061-second-closure-recertified-complete.graph.json \
  --source-vtx third_party/CNP-SAT/vtx/529.vtx \
  --source-edge third_party/CNP-SAT/edge/529.edge
```

注意：作者仓库漏收了与 `529-4-sbp.drat` 配套的 529 SBP CNF。直接拿
普通 CNF 检查会失败。`data/heule529/README.md` 给出跨同库实例对照、
只补三条单位子句的机械重建步骤、哈希和回归测试。

项目内虚拟环境可用：

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
```
