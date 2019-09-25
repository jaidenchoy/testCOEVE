## ----setup, include=FALSE, echo=FALSE------------------------------------
# knitr::knit_hooks$set(optipng = knitr::hook_optipng)
# knitr::opts_chunk$set(optipng = '-o7')

knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(fig.align = "center")
knitr::opts_chunk$set(fig.width = 12)
knitr::opts_chunk$set(fig.height = 6)

library(immunarch)
source("../R/testing.R")
# data(test)
# immdata = test_make_db(twb, data.frame(Sample = c("A", "B", "C", "D"), Twin = c("TwA", "TwA", "TwB", "TwB"), stringsAsFactors = F))
immdata = load_test_data()
# immdata$data = immdata$data[c(1:4, 12:16)]
# immdata$meta = immdata$meta[c(1:4, 12:16), ]

## ----gene-usage----------------------------------------------------------
gene_stats()

## ------------------------------------------------------------------------
# Next four function calls are equal. "hs" is from the "alias" column.
imm_gu = geneUsage(immdata$data, "hs.trbv")
# imm_gu = geneUsage(immdata$data, "HomoSapiens.trbv")
# imm_gu = geneUsage(immdata$data, "hs.TRBV")
# imm_gu = geneUsage(immdata$data, "HomoSapiens.TRBV")

imm_gu

## ---- message=F, fig.width=15, fig.height=8------------------------------
imm_gu = geneUsage(immdata$data, "hs.trbv", .norm = T, .ambig = "exc")

vis(imm_gu, .plot = "hist", .grid = T)

## ---- message=F, warning=FALSE, fig.width=12, fig.height=5---------------
vis(imm_gu, .plot = "hist", .grid = F, .by = "Status", .meta = immdata$meta)

## ---- message=F, warning=FALSE, fig.width=12, fig.height=5---------------
vis(imm_gu, .by = "Status", .meta = immdata$meta, .plot = "box")

## ---- message=F, fig.width=15, fig.height=8, warning=FALSE---------------
imm_gu = geneUsage(immdata$data, "hs.trbv", .norm = T, .ambig = "exc")

vis(imm_gu, .plot = "tree")

## ---- warning=FALSE, fig.width=12, fig.height=5--------------------------
imm_gu = geneUsage(immdata$data, "hs.trbv", .norm = T, .ambig = "exc")

imm_gu_js = geneUsageAnalysis(imm_gu, .method = "js", .verbose = F)
imm_gu_cor = geneUsageAnalysis(imm_gu, .method = "cor", .verbose = F)

gridExtra::grid.arrange(vis(imm_gu_js, .title = "Gene usage JS-divergence", .leg.title = "JS", .text.size=1.5), vis(imm_gu_cor, .title = "Gene usage correlation", .leg.title = "Cor", .text.size=1.5), ncol = 2)

## ---- warning=FALSE, fig.width=10, fig.height=4--------------------------
imm_gu_js[is.na(imm_gu_js)] = 0

vis(geneUsageAnalysis(imm_gu, "cosine+hclust", .verbose = F))

#vis(geneUsageAnalysis(imm_gu, "js+dbscan", .verbose = F))

## ---- message=F, warning=F, fig.width=12, fig.height=4-------------------
imm_cl_pca = geneUsageAnalysis(imm_gu, "js+pca+kmeans", .verbose = F)
imm_cl_mds = geneUsageAnalysis(imm_gu, "js+mds+kmeans", .verbose = F)
imm_cl_tsne = geneUsageAnalysis(imm_gu, "js+tsne+kmeans", .perp = .01, .verbose = F)

grid.arrange(vis(imm_cl_pca, .plot = "clust"), vis(imm_cl_mds, .plot = "clust"), vis(imm_cl_tsne, .plot = "clust"), ncol = 3)

## ---- message=F, warning=F, fig.width=8, fig.height=4--------------------
imm_cl_pca2 = geneUsageAnalysis(imm_gu, "js+pca+kmeans", .k = 3, .verbose = F)
vis(imm_cl_pca2)

## ----spectr, fig.width=12, fig.height=4----------------------------------
p1 = vis(spectratype(immdata$data[[1]], .quant = "id", .col = "aa", .gene = "v"))
p2 = vis(spectratype(immdata$data[[1]], .quant = "count", .col = "aa", .gene = "v"))

grid.arrange(p1, p2, ncol = 2)

