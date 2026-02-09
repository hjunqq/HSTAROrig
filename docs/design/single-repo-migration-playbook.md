# 单仓双轨迁移执行手册

## 1. 结论

当前阶段采用单仓迁移，不拆新仓库。

原因：
- 需要长期并行维护 legacy 与 next。
- 新旧结果对比、回归守门、回退开关在同仓库最容易执行。
- 迁移中期之前拆仓会增加同步与发布复杂度。

## 2. 仓库结构建议

在现有仓库内增加以下逻辑分区（不是强制一次到位）：

- `HSTAR/`：legacy 主代码（持续可发布）
- `next/`：新架构代码（按切片迁移）
- `migration/baseline/`：算例清单、容差、对比结果模板
- `docs/design/`：架构与迁移文档

## 3. 分支策略（单仓）

- `main`：稳定分支（生产可发布）
- `migration/<milestone>-<slice>`：迁移切片分支
- `hotfix/legacy-*`：legacy 紧急修复

规则：
1. 迁移改动只允许单切片提交，禁止打包大改。
2. `hotfix` 先修 `legacy`，再补到 `next` 适配层（同日完成）。
3. PR 必带“老跑 + 影子跑”差异报告。

## 4. 开关与回退策略

- 保留 legacy 为默认执行路径。
- 对每个迁移切片提供一个开关（编译时或运行时）。
- 一旦回归失败，立即切回 legacy 路径，不做线上试错。

## 5. M1 执行步骤（本周）

1. 建立黄金算例清单
- 编辑 `migration/baseline/cases.csv`。
- 至少覆盖：线弹性、弹塑性、接触、温度相关场景。

2. 固化容差
- 编辑 `migration/baseline/tolerances.yaml`。
- 分别设置位移、应力、残差、迭代步、耗时阈值。

3. 建立对比输出格式
- 每次对比都输出同一格式：`case_id, metric, legacy, next, delta, threshold, result`。

4. 设置放行规则
- 硬失败：位移/应力/残差超阈值。
- 软失败：耗时超阈值但数值通过（可临时放行并跟踪优化）。

## 6. 里程碑推进节奏

- M1：基线和门禁
- M2：adapter/facade/shadow-run
- M3：输入/错误/输出横切迁移
- M4：统一 newton driver（先静力线弹性）
- M5：首批元素与材料迁移
- M6：灰度切换与治理固化

## 7. 何时考虑拆仓

仅在满足以下条件后评估：
- `next` 已可独立构建与发布
- 70% 以上核心路径已迁移并稳定
- 团队需要独立发布节奏与权限边界

在此之前，坚持单仓双轨。

## 8. M1 工具命令（结果对比）

1. 复制模板：
- `Copy-Item migration/baseline/results-template.csv migration/baseline/results.csv`

2. 填入每个算例的新旧指标后运行：
- `python migration/baseline/compare_results.py --results migration/baseline/results.csv`

3. 查看输出：
- 明细：`migration/baseline/report.csv`
- 结论：`migration/baseline/summary.md`
