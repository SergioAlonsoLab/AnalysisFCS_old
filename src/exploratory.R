# exploratory

setwd("~/Documents/FACS/")
source("src/analyzeFCS.R")

library(ggrepel)
library(umap)


cutDensity <- function(x,cut=0.05,...) {
  x <- na.omit(x)
  density(x,from=max(x)*(cut),to=max(x)*(1-cut),...)
}

fcsfiles <- dir(path = "Ciclo Celular DAPI/",pattern = "*.fcs")
fcsfiles # view the file names

fcsSet <- read.flowSet(fcsfiles,path = "Ciclo Celular DAPI/") # reads a set of FCS files
fcsSet <- fcsSet@frames %>% as.list # transforms the enviroment in a list of fcs objects
typeof(fcsSet)

colpalette <- c("darkblue","lightblue","yellow","orange","red")

peaks <- function(x) {
  which(diff(sign(diff(x)))<0) + 1
}

valleys <- function(x) {
  which(diff(sign(diff(x)))>0) +1
}


sizeAndDAPI <- function(fcs,plot.results=T) {
  x <- exprs(fcs) %>% as.data.frame()
  sampleName <- fcs@description$`TUBE NAME`
  time <- x$Time * as.numeric(fcs@description$`$TIMESTEP`)
  events <- nrow(x)
  densitySize <- cutDensity(x$`FSC-A`,0.1,adjust=1/2)
  densityPBA <- cutDensity(x$`Pacific Blue-A`,0.1,adjust=1/2)
  peaksSize <- peaks(densitySize$y)
  peaksPBA <- peaks(densityPBA$y)
  
  highestPeakSize <- which.max(densitySize$y[peaksSize])
  highestPeakPBA <- which.max(densityPBA$y[peaksPBA])
  
  modalSize <- densitySize$x[peaksSize][highestPeakSize]
  modalPBA <- densityPBA$x[peaksPBA][highestPeakPBA]
  
  if(plot.results) {
    par(mfrow=c(1,2))
    plot(densitySize,main=sampleName,xlab="FSA-A")
    points(densitySize$x[peaksSize],densitySize$y[peaksSize],pch=19)
    abline(v=modalSize,lty=2)
    
    plot(densityPBA,main=sampleName,xlab="Pacific Blue-A")
    points(densityPBA$x[peaksPBA],densityPBA$y[peaksPBA],pch=19)
    abline(v=modalPBA,lty=2)
  }
  
  invisible(list(sampleName=sampleName,
                 events=events,
                 time=max(time),
                 speed=events/max(time),
                 modalSize=modalSize,
                 modalPBA=modalPBA))
  
}


results <- sapply(fcsSet,sizeAndDAPI,plot.results=F) %>% t %>% as.data.frame
for(i in 1:ncol(results)) results[,i] <- unlist(results[,i])
results$date <- sapply(fcsSet,keyword,"$DATE") %>% unlist 
results$P11V <- sapply(fcsSet,keyword,"$P11V") %>% unlist %>% as.numeric()



# add OV90

fcsSet2 <- read.flowSet(path = "FACS OV-90/",pattern = "*.fcs")
fcsSet2 <- fcsSet2@frames %>% as.list()

results2 <- sapply(fcsSet2,sizeAndDAPI,plot.results=F) %>% t %>% as.data.frame
for(i in 1:ncol(results2)) results2[,i] <- unlist(results2[,i])

results2$date <- sapply(fcsSet2,keyword,"$DATE") %>% unlist 
results2$P11V <- sapply(fcsSet2,keyword,"$P11V") %>% unlist %>% as.numeric()

results$cells <- "LS174T"
results2$cells <- "OV90"


results3 <- rbind(results,results2)

g0 <- ggplot(results) + 
  geom_point(aes(fill=unlist(foo2$minRange)),pch=21,size=3) 

g0 + aes(events,time)
g0 + aes(events,speed)
g0 + aes(events,modalSize)
g0 + aes(speed,modalSize)
g0 + aes(events,modalPBA)
g0 + aes(speed,modalPBA)
g0 + aes(modalSize,modalPBA) 


g1 <- ggplot(as.data.frame(exprs(fcsSet[[128]])))
dotColor <- densCols(g1$data$`FSC-A`,g1$data$`SSC-A`,colramp = myPalette)
g1 + aes(x=`FSC-A`,y=`SSC-A`) +  geom_point(pch=19,col=dotColor) + geom_density2d(col="black") 
g1 + aes(x=`FSC-A`,y=`SSC-A`) + geom_hex(bins=150)

#g1 + aes(`Pacific Blue-A`,`Pacific Blue-W`) + geom_density2d() + geom_point(pch=".")
dotColor <- densCols(g1$data$`Pacific Blue-A`,g1$data$`Pacific Blue-W`,colramp = myPalette)
g1 + aes(`Pacific Blue-A`,`Pacific Blue-W`) + geom_point(col=dotColor) + geom_density2d()




# PCA & umap test

foo <- as.data.frame(fcsSet$`120219_C++ mismo vPMT 4++.fcs`@exprs)


mySmoothScatter(foo$`Pacific Blue-A`,foo$`Pacific Blue-W`)

umap(foo[,1:12]) -> umap0
mySmoothScatter(umap0$layout[,1],umap0$layout[,2],mask=!selected)
gate3 <- locator(); polygon(gate3)
selected <- pointsInGate(umap0$layout[,1],umap0$layout[,2],gate3)


mySmoothScatter(foo$`Pacific Blue-A`,foo$`Pacific Blue-W`,mask = !selected)

pca0 <- prcomp(foo[,1:12])
mySmoothScatter(pca0$x[,1],pca0$x[,2])
mySmoothScatter(pca0$x[,1],pca0$x[,3])

mySmoothScatter(foo$`FSC-A`,foo$`Pacific Blue-W`)

gate3 <- locator(); polygon(gate3,col=alpha("yellow",.2))
selected <- pointsInGate(pca0$x[,1],pca0$x[,3],gate3)
mySmoothScatter(foo$`Pacific Blue-A`,foo$`Pacific Blue-W`,mask = !selected)

drawGate <- function(x,y,newplot=T) {
  if(newplot) mySmoothScatter(x,y)
  g0 <- locator()
  polygon(g0,col=alpha("yellow",.2))
  inGate <- pointsInGate(x,y,g0)
  return(list(gate=g0,inGate=inGate))
} 

