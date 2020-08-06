# Intall required libraries, if nescessary

install_if_not_installed <- function(pkgname) {
  if(!requireNamespace(pkgname)) install.packages(pkgname) 
  else cat(pkgname,"already installed\n")
}

install_if_not_installed("magrittr")
install_if_not_installed("ggplot2")
install_if_not_installed("mixtools")
install_if_not_installed("sp")
install_if_not_installed("BiocManager")

if(!requireNamespace("flowCore")) BiocManager::install("flowCore")

library(flowCore)
library(magrittr)
library(mixtools)

alpha <- ggplot2::alpha

myPalette <- colorRampPalette(c("darkblue","lightblue","yellow","orange","red"))

mySmoothScatter <- function(x,y,cex=.3,mask=NULL,...) {
  require(ggplot2)
  colors <- densCols(x,y,colramp = myPalette)
  if(!is.null(mask)) colors[mask] <- "grey"
        
  plot(x,y,pch=19,col=colors,cex=cex,xlab=NA,ylab=NA,...)
}

pointsInGate <- function(x,y,gate) {
  
  sp::point.in.polygon(x,y,gate$x,gate$y)
  
}


defineGates <- function(fcs1) {
  par(mfrow=c(1,1))
  x <- exprs(fcs1) %>% as.data.frame()
  
  # Define gate 1
  g1_defined <- g2_defined <- F
  
  while(!g1_defined) {
    mySmoothScatter(x$`FSC-A`,x$`SSC-A`)
    mtext("Define Gate 1",3,-1.5)
    title(description(fcs1)[["TUBE NAME"]],xlab="FSC-A",ylab="SSC-A")
    gate1 <- locator()
    polygon(gate1)
    
    n=sum(pointsInGate(x$`FSC-A`,x$`SSC-A`,gate1))
    print(n)
    mtext(sprintf("Total events: %i\nIn Gate_1: %i",nrow(x),n),3,-3.5)
    
    if(menu(choices = c("Gate_1 OK","Redraw Gate_1"),graphics = F)==1) g1_defined <- T
  }
  
  # Define gate 2
  
  while(!g2_defined) {
    mySmoothScatter(x$`Pacific Blue-A`,x$`Pacific Blue-W`)
    mtext("Define Gate 2",3,-1.5)
    title(description(fcs1)[["TUBE NAME"]],xlab="Pacific Blue-A",ylab="Pacific Blue-W")
    gate2 <- locator()
    polygon(gate2)
    
    n=sum(pointsInGate(x$`Pacific Blue-A`,x$`Pacific Blue-W`,gate2))

    mtext(sprintf("Total events: %i\nIn Gate_2: %i",nrow(x),n),3,-3.5)
    
    if(menu(choices = c("Gate_2 OK","Redraw Gate_2"),graphics = F)==1) g2_defined <- T
  }
  
 list(gate1=gate1,gate2=gate2)
  
}


analyzeFCS <- function(fcs1,gates,sampleName=NULL) {
  
  if(is.null(sampleName)) sampleName <- fcs1@description$`TUBE NAME`
  
  x <- exprs(fcs1) %>% as.data.frame()
  
  par(mar=c(5,5,2,5))
  
  x$Time <- x$Time * as.numeric(fcs1@description$`$TIMESTEP`) # adjust time to seconds
  
  timebins <- round(x$Time,0)+1
  speed <- data.frame(table(timebins))
  speed$timebins <- as.numeric(speed$timebins)

  mySmoothScatter(x$Time,x$`FSC-A`)
  title(sampleName,xlab="Time",ylab="FSC-A")
  
  par(new=T)
  plot(as.numeric(speed$timebins),
       speed$Freq,
       type="S",
       axes=F,
       bty="none",
       xlab=NA,
       ylab=NA,
       col="orange",
       lwd=1)
  
  abline(h=mean(speed$Freq),lty=2)
         
  axis(4)
  mtext("events/s",4,3)
  
  par(mar=c(5,5,2,1))
  
  mySmoothScatter(x$`FSC-A`,x$`SSC-A`)
  title(sampleName,xlab="FSC-A",ylab="SSC-A")
  polygon(gates$gate1)
  inGate1 <- pointsInGate(x$`FSC-A`,x$`SSC-A`,gates$gate1)
  mtext(sprintf("in gate: %i  ",sum(inGate1)),1,line = -1.5,adj = 1,cex=par()$cex*.8)
  mtext(sprintf("  events: %i",length(inGate1)),3,line = -1.5, adj = 0,cex=par()$cex*.8)
  
  mySmoothScatter(x$`Pacific Blue-A`,x$`Pacific Blue-W`)
  title(sampleName,xlab="Pacific Blue-A",ylab="Pacific Blue-W")
  polygon(gates$gate2)
  inGate2 <- pointsInGate(x$`Pacific Blue-A`,x$`Pacific Blue-W`,gates$gate2)
  mtext(sprintf("in gate: %i  ",sum(inGate2)),1,line = -1.5,adj = 1,cex=par()$cex*.8)
  
  fal <- function(x) c(x[1],x,x[length(x)]) # firs-all-last (fal)
  
  cellcycle <- density(x$`Pacific Blue-A`[inGate1 & inGate2])
  
  plot(cellcycle,main=sampleName,xlab="Pacific Blue-A")
  polygon(fal(cellcycle$x),c(0,cellcycle$y,0),col=alpha("red",.2))
  
  priorMu <- cellcycle$x[which.max(cellcycle$y)] * c(1,2)
  abline(v=priorMu,col="lightblue")
  
  
  good_fit <- F
  tries <- 0
  
  while(!good_fit & tries < 10) {
    fitted <- normalmixEM(x$`Pacific Blue-A`[inGate1 & inGate2],k=2,mu = priorMu)
    if(abs(fitted$mu[1] - priorMu[1])/priorMu[1] < .05) good_fit <- T
    tries <- tries + 1
  }
  
  fit1 <- dnorm(cellcycle$x,fitted$mu[1],fitted$sigma[1]) * fitted$lambda[1]
  lines(cellcycle$x,fit1)
  polygon(fal(cellcycle$x),c(0,fit1,0),col = alpha("green",.2))
  abline(v=fitted$mu[1],lty=2)

  results <- list(sampleName=sampleName,
                  events=nrow(x),
                  time=max(x$Time),
                  speed=mean(speed$Freq),
                  selected=cellcycle$n,
                  G01=fitted$lambda[1],
                  G01_at=fitted$mu[1])
  
  mtext(sprintf("Events: %i\nTime: %1.1f s\nSpeed: %1.1f events/s\nIn gates: %i\nG0/1: %1.1f%%\nat: %1.1f",
                results$events,results$time,results$speed,results$selected,results$G01*100,results$G01_at),
        side=3,
        line=-6,
        adj=.9,
        cex=par()$cex*.9)
  
  invisible(results)
}


analyzeFCS_set <- function(fcsSet,gates) {
  results <- sapply(fcsSet,analyzeFCS,gates) %>% t %>% as.data.frame
  for(i in 1:ncol(results)) results[,i] <- unlist(results[,i])
  invisible(results)
}


# ggplot version 


