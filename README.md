# pigRutils

`pigRutils` 是我个人日常分析里积累的一组 R 小工具，主要服务于生物信息学数据处理与可视化，尤其偏向单细胞、GO 富集分析和一些高频辅助操作。

这个仓库首先是“自己会反复用到的脚本工具箱”，然后才是一个轻量 R 包。因此它的特点很直接：

- 函数覆盖面杂，但都来自真实分析场景
- 重点在提高日常处理效率，而不是做成大而全框架
- 部分函数依赖特定生信对象或分析习惯，使用前建议先看函数参数

## 安装

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github("Muchmorepig/pigRutils")
```

如果你只想拿其中某几个函数，也可以直接克隆仓库后按脚本方式 source 使用。

## 适用内容

目前仓库里的函数大致分成几类：

- GO 富集分析与结果导出
- 单细胞降维图和质控图绘制
- 矩阵标准化、平滑和按组汇总
- 一些分析流程里的零散辅助函数

## 常用函数

### GO 富集分析

- `pGO()`: 对一组基因做 GO enrichment
- `batch_GO()`: 对多个基因集批量做 GO 分析
- `plotGO()`: 绘制 GO 结果
- `gobar()`: 生成 GO bar plot
- `export_GO()`: 导出 GO 结果到文件

### 单细胞可视化

- `embSCdim()`: 绘制 Seurat 对象的降维图
- `qc_plot()`: 绘制常见单细胞 QC 图
- `emb_top_middle_low_cells()`: 按指定基因集合信号把细胞分层并展示
- `umapForArchR()`: 为 ArchR 对象绘制 UMAP

### 数值与矩阵处理

- `mf_zscore()`: 按行做 Z-score
- `mf_minmax()`: 做 min-max 归一化
- `mf_mean()`: 按分组计算矩阵均值
- `smooth_byRow()`: 按行平滑数据
- `calc_RowMeansByColumnBlocks()`: 按列分块后计算行均值

### 其他辅助函数

- `select_colors()`: 常用配色方案
- `calculate_pdf_dimensions()`: 根据图数量估算 PDF 尺寸
- `doPCA()`: 快速执行 PCA 并整理结果
- `calcConDist()`: 比较不同条件下细胞类型分布

## 示例

下面是一个批量 GO 分析并导出结果的简单例子。

```r
library(clusterProfiler)
library(org.At.tair.db)
library(readxl)
library(dplyr)
library(purrr)
library(pigRutils)

orgdb <- org.At.tair.db

ff <- "./dataLib/marker/marker_wilcoxon_pts.xlsx"
sheets <- excel_sheets(ff)

go_data <- map(sheets, ~ read_excel(ff, sheet = .x) %>%
  filter(logfoldchanges >= 1) %>%
  pull(names))

names(go_data) <- sheets

res <- batch_GO(go_data, orgdb = orgdb, keyType = "TAIR")
export_GO(res, file = "./go_allcluster.xlsx")
```

## 使用说明

- 这是个人工具仓库，函数风格不一定完全统一
- 某些函数依赖 `Seurat`、`ArchR`、`clusterProfiler` 等包
- 某些函数面向我自己的分析习惯封装，未必适合所有项目

如果你准备长期复用，建议先阅读对应函数源码和参数说明，再纳入自己的流程。

## 仓库定位

这个项目不是通用 R 包模板，更接近一个持续增长的个人分析工具箱。保留它的目的，是把日常分析里重复出现的代码沉淀下来，减少复制粘贴，方便后续直接调用。

## `calcConDist()`：细胞组成与统计单位

输入一行一个细胞，可为 `data.frame` 或 `data.table`，列名由参数指定。
推荐先输出每个样本的细胞数、比例与描述性效应量：

```r
composition <- calcConDist(
  seu@meta.data,
  colname.cluster = "celltype",
  colname.sample = "sample_id",
  colname.condition = "sex",
  method = "freq"
)
composition$sample.freq.tb  # 每样本×细胞类型；未检出类型补 0
composition$res            # 条件均值、logFC、样本数和检验状态
```

- `freq` 始终先按样本计算比例，再对样本等权求均值；不按细胞数加权。
  它要求每个样本只对应一个条件。不会给不存在的样本×条件组合补零。
- 默认 `sample.inference = FALSE`：p 值和 FDR 为 `NA`。仅当样本 ID 确实代表
  **独立生物学重复**时，显式设置 `sample.inference = TRUE`，启用非配对
  Welch t 检验和 Wilcoxon 秩和检验。每个对比的两侧都至少需要 2 个样本；
  这是计算门槛，不是充分的实验设计或功效保证。
- 多于两个条件时，每个条件对其余所有样本比较，不是所有条件两两比较。
  多个参考条件按样本数加权。两条件的正反向结果都保留，BH 按全部结果行校正；
  t 与 Wilcoxon 各自校正，`char.sig` 沿用 t 检验 FDR。
- `chisq` 返回 observed/expected 比值矩阵，不是 fold change。
  `fisher` 对每个细胞类型×条件做 2×2“该类型/其他类型 × 该条件/其他条件”检验。
  两者都基于**细胞计数**，不能把细胞数量当作生物学重复数。
- `bySample = TRUE` 对 `chisq/fisher` 返回以样本名命名的列表，每个元素保持
  相应方法的返回结构。样本内只有一个条件时无法检验条件差异，结果为 `NA`。
  这与跨生物学重复的组间比较不同。对 `freq`，此参数不改变始终按样本计算的规则。
- `min.rowSum` 仅筛选 Fisher 报告目标（严格 `>`），检验背景保留所有细胞。
  这修正了旧实现中先删细胞类型、再改变背景的问题。
- 不添加伪计数。`logFC = log2(mean.x / mean.y)`，零分母可以产生 `Inf/-Inf`，
  `0/0` 为 `NA`。退化检验、常数数据和数值异常通过 `status` 或
  `status.t/status.w` 说明，不伪造 p=1。
- 为兼容下游代码，freq 输出仍使用 `sampleID`、`majorCluster`、`loc` 等固定字段名；
  它们不要求输入有这些列。宽表首列是类型名，其后是条件名；条件名与首列同名时
  会保留重复列名，应按位置访问或使用长表 `res/count.dist`。

**Tak-1 与 Tak-2 的 11 天数据（各一个样本）**：使用默认 `freq` 做描述比较，
报告“本次捕获并通过 QC 的细胞中，该类型比例为……”及比例差或 logFC。
不能据此声称群体层面的显著性别差异。Tak-1/Tak-2 遗传背景、建库批次和性别
也可能混杂；随机拆细胞、下采样或技术重复不能解决独立重复缺失。
细胞回收比例还受解离、捕获、QC 和注释影响，不等同于组织内真实细胞丰度。

有独立重复时，`freq` 的两个单变量检验仅是简单探索；配对设计、批次协变量、
过度离散及组成数据约束需要适当模型。正式差异丰度分析可进一步考虑样本级
计数模型或组成模型，不能靠本函数的 p 值自动处理这些问题。

独立运行该函数的回归测试（无需安装其他生信依赖）：

```sh
Rscript -e 'install.packages(c("data.table", "testthat"))'
Rscript tools/test-calcConDist.R
```

详细审查、兼容性变化与统计解释见 [审查报告](docs/calcConDist-review.md)。
