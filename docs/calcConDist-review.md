# calcConDist() 审查与修复说明

审查日期：2026-11-02。上游 main 基线：
`c0b5eecfa0e6bdec00795af92eb88f6d5e821e45`。
源码：https://github.com/Muchmorepig/pigRutils/blob/c0b5eecfa0e6bdec00795af92eb88f6d5e821e45/R/mf_calcConDist.R

## 已确认问题

1. chisq/fisher 提前 return，使 bySample 永不可达。改为先分派 bySample，
   并使用 lapply 保留具名列表，避免 apply 简化或丢失矩阵维度和标签。
2. freq 写死输入列名，长宽表混用，N 在转换后语义不再可靠；重建样本×类型
   计数表，再按每个观察到的样本总细胞数归一化。
3. freq 只检查一侧样本数、没有对缺失类型补零，并且未明确独立重复单位。
   现在两侧都检查、样本等权，默认关闭推断，显式 opt-in 只表示调用者声明
   独立性；函数无法从 sample ID 自动验证独立性。
4. 缺少方法、标志位、列名、缺失元数据和阈值验证；现在提前报错，避免静默
   丢弃 NA 而改变分析分母。不会修改调用者 data.table。
5. 零边际、单条件、常数样本会引发无意义检验或数值错误；现在返回 NA 和
   可检查的状态。保留卡方近似警告，不把错误自动替换成其他检验。
6. Fisher 的过滤原来会改变“其他类型”背景；现在 min.rowSum 仅筛报告目标。
   非数字类型名通过显式标签匹配，不靠隐式数字索引。
7. 移除脚本顶层 library(data.table)，使用限定命名空间，避免 source 时改变
   搜索路径。DESCRIPTION 补充本函数的 Imports 与测试 Suggests。

## 接口与行为约定

前七个参数顺序和默认值不变，末尾新增 `sample.inference = FALSE`。
不改变 pooled chisq 的矩阵形态以及 Fisher 的四个主要字段。
chisq 新增 attributes：counts、expected、p.value、status、test.warning；
Fisher 长表保留 V1/V2、p.value、OR、p.adj，增加 count/status。
原来无法工作的 bySample 现在返回“样本名 -> 与 pooled 同结构结果”的列表。
旧的 bySample 代码原本不可达，不能视为稳定返回契约。

freq 保留 res、dist.logFC.tb、dist.FDR.tb、dist.charSig.tb，并增加 sample.freq.tb。
输出历史列名 majorCluster/loc/sampleID 固定，输入列名完全由参数控制。
res 增加 mean.x/mean.y、n.x/n.y、status.t/status.w。
即使不可检验，也保留描述性结果，不再丢弃小样本条件。

**有意的行为变化**：freq 不默认输出推断 p 值；需要显式声明独立重复。
Fisher 过滤后背景不同于旧实现；BH 只覆盖报告的目标×条件，不含被过滤类型。
有零边际时保留维度；chisq 在正边际子矩阵检验并回填，其余位置 NA。
chisq 至少需要两个正行边际及两个正列边际，不退化为拟合优度检验。
单个 Fisher 2×2 有任一零边际就返回 NA。真零格不加伪计数，OR 可以为 0/Inf。

factor 未使用的类型水平保留：freq 比例为零，0/0 的 logFC 为 NA；
未使用的条件水平只在细胞计数结果中保留为零边际，不作为 freq 的真实条件。
未使用的样本水平不生成虚构样本。
缺失/空白元数据被拒绝，使用者应先决定如何处理，而不是让 table 隐式丢行。
宽表保留原始条件名，包括与标识列相撞的名称；推荐长表规避重复列名。

## 三种统计目标

| 目标 | 数据单位 | 当前函数能力 | 结论范围 |
|---|---|---|---|
| 描述捕获细胞组成 | 细胞计数与每个库比例 | freq 默认模式 | 本次样本、当前 QC 和注释下的观察结果 |
| 检验捕获细胞计数的关联 | 细胞 | chisq / Fisher | 独立细胞抽样假设下的关联，不能估计生物学重复间变异 |
| 比较条件间组成 | 独立生物学样本 | freq 显式 opt-in 的简单非配对检验 | 依赖独立性、设计、分布假设和足够重复；不能自动去除混杂 |

freq 的 Welch t 检验针对均值，Wilcoxon 检验针对分布/秩位置，二者不检验
完全相同的零假设。比例有边界、稀疏与组成依赖，检验的可计算性不等于模型适用。
Wilcoxon 使用 `exact=FALSE` 的含 ties 校正的近似；两个样本一组时功效尤其有限。
两条件的互为反向结果仍保留并纳入 BH，以保持输出形式；不是两次独立验证。
多条件是 condition-vs-rest，参考样本来自不同条件时，异质性可能使简单检验不适用。

## 地钱 Tak-1 / Tak-2，11 天，每组一个样本

建议调用（替换为实际元数据字段）：

```r
r <- calcConDist(seu@meta.data,
  colname.cluster = "celltype", colname.sample = "sample_id",
  colname.condition = "sex", method = "freq")
r$sample.freq.tb
r$res[, c("majorCluster", "loc", "mean.x", "mean.y", "logFC", "n.x", "n.y"),
      with = FALSE]
```

可以报告：两个已测文库中某类型的细胞数、占比、占比差（mean.x - mean.y）及
log2 比例比；描述为“在本次 Tak-1/Tak-2 样本中观察到……”，作为后续验证候选。
不能报告：由细胞数得到的 p/FDR 证明可推广的雄雌群体差异、估计群体间变异、
把未检出解释为某性别绝对缺失，或把比例变化解释为绝对细胞数改变。
chisq 的 observed/expected、Fisher OR 与频率 logFC 是不同效应量，不能混称。

Tak-1 和 Tak-2 的比较还包含遗传背景差异；单一雄株系与单一雌株系不能
单独识别性别与株系效应。即使同一株系增加独立培养/取材重复，也应将结论
限定为这两个株系；推广为性别效应还需更广泛株系/遗传设计。
文库批次若与性别完全混杂，统计模型也无法分离它们。
随机把细胞拆成几个“样本”是伪重复。下采样、细胞 bootstrap 可以评估捕获层面
稳定性，但不能提供生物学重复间方差，批量 RNA 重复也不能代替组成分析的单细胞重复。

补充独立培养和取材的重复，并尽可能平衡解离、建库、测序批次后，才适合进行
样本层面的条件推断。实际重复数应根据变异、效应量与功效确定；代码的 >=2
门槛只是阻止显然无自由度的检验，不是建议实验只做两个重复。
正式分析可按问题采用样本级组成模型、beta-binomial/多项式模型，或 Milo 等
差异丰度框架；不能期望换模型自动弥补当前 n=1 的设计局限。

## 检验与依赖范围

`Rscript tools/test-calcConDist.R` source 单一函数并运行 testthat，避免加载包内
与本修复无关的 Seurat、clusterProfiler、SummarizedExperiment 等依赖。
同时提供标准 tests/testthat.R，可在完整依赖环境运行包测试。
仓库原 DESCRIPTION 使用非标准 `import:`，License 也是占位文本，NAMESPACE
存在多个未在规范 Imports 声明的全包导入；这次只声明本函数需要的 data.table。
这些全包问题不能通过单函数测试证明已修复，完整 R CMD check 不应宣称通过。

## 参考

- R chisq.test：https://stat.ethz.ch/R-manual/R-patched/library/stats/html/chisq.test.html
- R fisher.test：https://stat.ethz.ch/R-manual/R-patched/library/stats/html/fisher.test.html
- R wilcox.test：https://stat.ethz.ch/R-manual/R-patched/library/stats/html/wilcox.test.html
- Dann et al., Milo 原始论文：https://doi.org/10.1038/s41587-021-01033-z
- Büttner et al., scCODA 原始论文：https://doi.org/10.1038/s41467-021-27150-6
- Zimmerman et al., 伪重复研究：https://doi.org/10.1038/s41467-021-21038-1

## 本次验证结果

- R 4.3.3，data.table 1.14.10，testthat 3.2.1。
- 13 组 testthat 测试全部通过；无失败、无未预期警告、无跳过。
- chisq/Fisher/Welch/Wilcoxon 结果与 R stats 数值核对通过。
- 自定义列名、具名逐样本结果、2 vs 3 样本、1 vs 3 样本、全零类别、
  未使用 factor 水平、单条件、单类型、Inf/-Inf、常数数据、缺失元数据、
  多条件 condition-vs-rest、data.table 无原地修改均有回归覆盖。
- Rd 与 DESCRIPTION 解析通过，git diff --check 通过。
- 补丁在精确基线的干净工作树上 git apply --check 通过，应用后逐文件一致。
- 本次没有用户实际单细胞元数据，Tak-1/Tak-2 场景用相同设计结构的模拟计数验证；
  没有重新计算真实地钱细胞比例，也没有运行完整包 R CMD check。
