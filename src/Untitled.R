library(flowCore)
library(magrittr)

try(setwd("~/Documents/FACS/"))

files <- dir("FACS OV-90/",full.names = T)

fcs1 <- read.FCS(files[1])


