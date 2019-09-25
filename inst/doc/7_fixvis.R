## ----setup, include=FALSE, echo=FALSE------------------------------------
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(fig.align = "center")
knitr::opts_chunk$set(fig.width = 12)
knitr::opts_chunk$set(fig.height = 6)

library(immunarch)

## ----fixvis1, eval=F-----------------------------------------------------
#  data(immdata)
#  gu = geneUsage(immdata$data)
#  p = vis(gu)
#  fixVis(p)

## ----fixvis2, eval=F-----------------------------------------------------
#  fixVis()

## ----fixvis3, eval=T-----------------------------------------------------
knitr::include_app("https://immunomind.shinyapps.io/fixvis/", height = "800px")

