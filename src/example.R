# Example of using analyzeFCS functions

setwd("~/Documents/FACS/")
source("src/analyzeFCS.R")



fcsfiles <- dir(path = "FACS OV-90/",pattern = "*.fcs")
fcsfiles # view the file names


fcsSet <- read.flowSet(fcsfiles,path = "FACS OV-90") # reads a set of FCS files
fcsSet <- fcsSet@frames %>% as.list # transforms the enviroment in a list of fcs objects
typeof(fcsSet)

# define gates on a particular experiment
defineGates(fcsSet$`OV-90_22072021_OV-90 Cell Line.fcs`) -> gates 

# analyze a single experiment
# generates 4 individual plots
analyzeFCS(fcsSet$`OV-90_22072021_OV-90 Cell Line.fcs`,gates) -> results 

# same analysis, but combining the four plots
par(mfrow=c(2,2))
analyzeFCS(fcsSet$`OV-90_22072021_OV-90 Cell Line.fcs`,gates) -> results 

# analyze all FCS using the gates defined in line 15
par(mfrow=c(2,2))
analyzeFCS_set(fcsSet,gates) -> results_all

# plot speed vs peak at G0/1

library(ggplot2)
install_if_not_installed("ggrepel")
library(ggrepel)

ggplot(results_all) +
  aes(G01_at,speed) +
  geom_point(fill="lightblue",size=3,pch=21) +
  geom_smooth(method="lm",formula=y ~ I(1/x)) +
  geom_label_repel(label=results_all$sampleName,alpha=.8) +
  xlab("Position of G0/1 Peak") +
  ylab("Flow speed (events/s)")

summary(lm(speed ~ I(1/G01_at),results_all))

pdf("pdf_test.pdf",8,8)
par(mfrow=c(2,2))
analyzeFCS_set(fcsSet,gates) -> results_all
dev.off()

defineGates(fcs1)

