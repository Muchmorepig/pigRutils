# pigRutils

## install
```R
devtools::install_github("Muchmorepig/pigRutils")
```
## 安装
```R
devtools::install_github("Muchmorepig/pigRutils")
```

## Example
之前保存好的单细胞FindMarker结果xlsx，每个Sheet对应一个Cluster，以下操作包括读取后GO分析，然后保存
```R
# 加载工具
library(clusterProfiler)
library(org.At.tair.db)
library(readxl)
library(dplyr)
library(purrr)
library(pigRutils)

orgdb <- org.At.tair.db

# 读取文件
ff <- "./dataLib/marker/marker_wilcoxon_pts.xlsx"
sheets <- excel_sheets(ff)

go_data <- map(sheets, ~ read_excel(ff, sheet = .) %>%
  filter(logfoldchanges >= 1) %>%
  pull(names))
names(go_data) <- sheets

# GO
res <- batch_GO(go_data, orgdb = orgdb, keyType = "TAIR")
# 保存GO结果
export_GO(res, file = "./go_allcluister.xlsx")

```
