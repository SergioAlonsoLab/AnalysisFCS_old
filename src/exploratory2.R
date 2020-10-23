setwd("/imppc/labs/mplab/salonso/AnalysisFCS/")

source("src/analyzeFCS.R")

fcsSet <- read.flowSet(path="data")
fcsSet <- fcsSet@frames %>% as.list
fcsSet[[4]] <- NULL
names(fcsSet) <- names(fcsSet) %>% gsub("23102020_LS174T_SORTER_","",.) %>% gsub(".fcs","",.)


png("output/exploratory1.png",600,1600)

par(mfrow=c(5,2))

for(i in 1:4) {
  
  x <- exprs(fcsSet[[i]]) %>% as.data.frame
  
  mySmoothScatter(x$`Pacific Blue-A` / 1000,x$`Pacific Blue-W` / 1000)
  grid()
  title(names(fcsSet)[i],xlab="Pacific Blue A (x1000)",ylab="Pacific Blue W (x1000)")
  
  
}

k <- 25
m <- list()
for(i in c(1:4)) {
  
  x <- exprs(fcsSet[[i]]) %>% as.data.frame

  table(cut(x$`Pacific Blue-A`,k),
        cut(x$`Pacific Blue-W`,k)) -> foo
  foo <- foo[2:(k-1),2:(k-1)]
  image(log(foo+1))
  grid()
  title(names(fcsSet)[i])
  m[[i]] <- foo / sum(foo)
}

m <- sapply(m,c) %>% as.data.frame
names(m) <- names(fcsSet)

dim(m)

m <- m[apply(m,1,function(x) any(x!=0)),]
dim(m)
plot(hclust(dist(t(log(m+1e-10))),method = "average"),xlab=NA)
plot(hclust(dist(t(log(m+1e-10))),method = "complete"),xlab=NA)


dev.off()

sum(m$Tube_001 - m$Tube_004) 

dist(t(log(m+1e-10)))

foo <- exprs(fcsSet$Tube_004) %>% as.data.frame

ggplot(foo) + aes(`Pacific Blue-A`,`Pacific Blue-W`) + geom_point(aes(col=`FSC-A`),size=.5)
ggplot(foo) + aes(`Pacific Blue-A`,`FSC-A`) + geom_point(aes(col=`Pacific Blue-W`),size=.5)


