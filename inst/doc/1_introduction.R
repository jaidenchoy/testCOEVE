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

## ----eval=FALSE----------------------------------------------------------
#  install.packages("devtools", dependencies = T)
#  devtools::install_local("path/to/your/folder/with/immunarch.tar.gz", dependencies=T)

## ----eval=FALSE----------------------------------------------------------
#  # Load the data to the package
#  immdata = repLoad("path/to/your/folder/with/repertoires")
#  # If you folder contains metadata.txt file, immdata will have two elements:
#  # - immdata$data with a list of parsed repertoires
#  # - immdata$meta with the metadata file
#  
#  # Compute and visualise overlap statistics
#  ov = repOverlap(immdata$data)
#  vis(ov)
#  
#  # Cluster samples using K-means algorithm applied to the number of overlapped clonotypes
#  # and visualise the results
#  ov.kmeans = repOverlapAnalysis(ov, .method = "kmeans")
#  vis(ov.kmeans)
#  
#  # Compute and visualise gene usage with samples, grouped by their disease status
#  gu = geneUsage(immdata$data)
#  vis(gu, .by="Status", .meta=immdata$meta)
#  
#  # Compute Jensen-Shannon divergence among gene distributions of samples,
#  # cluster samples using the hierarchical clustering and visualise the results
#  gu.clust = geneUsageAnalysis(gu, .method = "js+hclust")
#  vis(gu.clust)
#  
#  # Compare diversity of repertoires and visualise samples, grouped by two parameters
#  div = repDiversity(immdata$data, .method = "chao1")
#  vis(div, .by=c("Status", "Treatment"), .meta=immdata$meta)
#  
#  # Manipulate the visualisation of diversity estimates to make the plot publication-ready
#  div.plot = vis(div, .by=c("Status", "Treatment"), .meta=immdata$meta)
#  fixVis(div.plot)

## ----eval=FALSE----------------------------------------------------------
#  data(immdata)

