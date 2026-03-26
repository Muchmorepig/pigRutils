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
