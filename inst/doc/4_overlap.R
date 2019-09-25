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

## ----overlap, message=F, warning=FALSE, fig.width=12, fig.height=7-------
imm_ov1 = repOverlap(immdata$data, .method = "public", .verbose = F)
imm_ov2 = repOverlap(immdata$data, .method = "morisita", .verbose = F)

grid.arrange(vis(imm_ov1), vis(imm_ov2, .text.size=1.5), ncol = 2)

vis(imm_ov1, "heatmap2")

## ---- warning=TRUE, fig.width=14, fig.height=8---------------------------
warning("TODO")

## ----overlap-1, warning=F, fig.width=8, fig.height=5---------------------
# Apply different analysis algorithms to the matrix of public clonotypes:
# "mds" - Multi-dimensional Scaling
repOverlapAnalysis(imm_ov1, "mds")
# "tsne" - t-Stochastic Neighbor Embedding
repOverlapAnalysis(imm_ov1, "tsne")

# Visualise the results
vis(repOverlapAnalysis(imm_ov1, "mds"))

## ----overlap-2, warning=F, fig.width=10, fig.height=5--------------------
# Apply different analysis algorithms to the matrix of public clonotypes:
# "mds" - Multi-dimensional Scaling
repOverlapAnalysis(imm_ov1, "mds")
# "tsne" - t-Stochastic Neighbor Embedding
repOverlapAnalysis(imm_ov1, "tsne")

# Visualise the results
vis(repOverlapAnalysis(imm_ov1, "mds"))

# Clusterise the MDS resulting components using K-means
vis(repOverlapAnalysis(imm_ov1, "mds+kmeans"))

## ---- warning=TRUE, fig.width=14, fig.height=8---------------------------
# Pass "nt" as the second parameter to build the public repertoire table using CDR3 nucleotide sequences
pr.nt = pubRep(immdata$data, "nt", .verbose = F)
pr.nt

## ---- warning=TRUE, fig.width=14, fig.height=8---------------------------
# Pass "aa+v" as the second parameter to build the public repertoire table using CDR3 aminoacid sequences and V alleles
# In order to use only CDR3 aminoacid sequences, just pass "aa"
pr.aav = pubRep(immdata$data, "aa+v", .verbose = F)
pr.aav

## ---- eval=FALSE, warning=TRUE, fig.width=14, fig.height=8---------------
#  # You can also pass the ".coding" parameter to filter out all noncoding sequences first:
#  pr.aav.cod = pubRep(immdata$data, "aa+v", .coding=T)

## ---- eval=FALSE, warning=TRUE, fig.width=14, fig.height=8---------------
#  # Create a public repertoire with coding-only sequences using both CDR3 amino acid sequences and V genes
#  pr = pubRep(immdata$data, "aa+v", .coding = T, .verbose = F)
#  
#  # Apply the filter subroutine to leave clonotypes presented only in healthy individuals
#  pr1 = pubRepFilter(pr, immdata$meta, c(Status = "C"))
#  
#  # Apply the filter subroutine to leave clonotypes presented only in diseased individuals
#  pr2 = pubRepFilter(pr, immdata$meta, c(Status = "MS"))
#  
#  # Divide one by another
#  pr3 = pubRepApply(pr1, pr2)
#  
#  # Plot it
#  p = ggplot() + geom_jitter(aes(x = "Treatment", y = Result), data=pr3)
#  p

