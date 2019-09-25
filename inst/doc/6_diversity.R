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

## ----diversity, fig.width=10, fig.height=4, warning=FALSE, message=FALSE----
# Compute statistics and visualise them
# Chao1 diversity measure
div_chao = repDiversity(immdata$data, "chao1")

# Hill numbers
div_hill = repDiversity(immdata$data, "hill")

# D50
div_d50 = repDiversity(immdata$data, "d50")

# Ecological diversity measure
div_div = repDiversity(immdata$data, "div")

p1 = vis(div_chao)
p2 = vis(div_chao, .by=c("Status", "Sex"), .meta=immdata$meta)
p3 = vis(div_hill, .by=c("Status", "Sex"), .meta=immdata$meta)

p4 = vis(div_d50)
p5 = vis(div_d50, .by="Status", .meta=immdata$meta)
p6 = vis(div_div)

gridExtra::grid.arrange(p1, p2, ncol = 2)
gridExtra::grid.arrange(p3, p6, ncol = 2)
gridExtra::grid.arrange(p4, p4, ncol = 2)

## ---- warning=F, fig.width=12, fig.height=4.5----------------------------
imm_raref = repDiversity(immdata$data, "raref", .verbose = F)

grid.arrange(vis(imm_raref), vis(imm_raref, .by="Status", .meta=immdata$meta), ncol=2)

