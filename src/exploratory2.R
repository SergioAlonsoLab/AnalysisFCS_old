setwd("/imppc/labs/mplab/salonso/AnalysisFCS/")

source("src/analyzeFCS.R")

fcsSet <- read.flowSet(path="data")
fcsSet <- fcsSet@frames %>% as.list
fcsSet[[4]] <- NULL
names(fcsSet) <- names(fcsSet) %>% gsub("23102020_LS174T_SORTER_","",.) %>% gsub(".fcs","",.)


filtered <- list()

for(i in 1:4) {
  
  x <- as.data.frame(exprs(fcsSet[[i]]))
  g1 <- gating(x$`FSC-A`,x$`SSC-A`)
  
  filtered[[i]] <- x[g1,]
  
}

names(filtered) <- names(fcsSet)

par(mfrow=c(2,2))
for(i in 1:4) {
  x <- filtered[[i]]
  mySmoothScatter(x$`Pacific Blue-A`,x$`Pacific Blue-W`)
  title(names(filtered)[[i]])
}


gates <- list()
par(mfrow=c(1,1))
for(i in 1:6) {
  plotGraph <- ifelse(i==1,T,F)
  gates[[i]] <- gating(filtered$Tube_001$`Pacific Blue-A`,filtered$Tube_001$`Pacific Blue-W`,plotGraph = plotGraph)
}

y <- filtered$Tube_001[NULL,]

for(i in 1:6) {
  
  x <- filtered$Tube_001[gates[[i]],]
  x$Gate <- paste("Gate",i)
  
  y <- rbind(y,x)
  
  
}

ggplot(y) + geom_point(aes(`Pacific Blue-A`,`Pacific Blue-W`,col=Gate),alpha=.5,size=2)
ggplot(y) + geom_density(aes(`FSC-A`,fill=Gate),alpha=.5) + ggtitle("Particle size")
ggplot(y) + geom_density(aes(`SSC-A`,fill=Gate),alpha=.5) + ggtitle("Particle granularity")

ggplot(y) + geom_point(aes(`FSC-A`,`SSC-A`,col=Gate),alpha=.4,size=2) + geom_density2d(aes(`FSC-A`,`SSC-A`))

ggplot(y) + geom_violin(aes(Gate,`SSC-A`,fill=Gate)) 
ggplot(y) + geom_violin(aes(Gate,`FSC-A`,fill=Gate))

ggplot(y) + geom_violin(aes(Gate,(`Pacific Blue-A` / `Pacific Blue-W`) ,fill=Gate))

# centroids ------

foo <- hclust(dist(filtered$Tube_001[,c("Pacific Blue-A","Pacific Blue-W")]))
cutree(foo,k=6)

# -----

png("output/exploratory1.png",600,1600)

par(mfrow=c(5,2))

for(i in 1:4) {
  
  x <- exprs(fcsSet[[i]]) %>% as.data.frame
  
  mySmoothScatter(x$`Pacific Blue-A` / 1000,x$`Pacific Blue-W` / 1000)
  grid()
  title(names(fcsSet)[i],xlab="Pacific Blue A (x1000)",ylab="Pacific Blue W (x1000)")
  
  
}

k <- 100
m <- list()
for(i in c(1:4)) {
  
  x <- filtered[[i]] %>% as.data.frame

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

foo <- exprs(fcsSet$Tube_001) %>% as.data.frame
gate1 <- gating(foo$`FSC-A`,foo$`SSC-A`)
foo <- foo[gate1,]

gate2 <- gating(foo$`Pacific Blue-A`,foo$`Pacific Blue-W`)
foo <- foo[gate2,]

ggplot(foo) + geom_density(aes(`Pacific Blue-A`),fill="lightblue")

ggplot(foo) + aes(`Pacific Blue-A`,`Pacific Blue-W`) + geom_point(aes(col=`FSC-A`),size=.5)






# some tests


# UMAP

library(umap)

u1 <- umap(y[,c("FSC-A","SSC-A")])

y$UMAP1 <- u1$layout[,1]
y$UMAP2 <- u1$layout[,2]

ggplot(y) + geom_point(aes(UMAP1,UMAP2,color=Gate))

pca1 <- prcomp(y[,1:10])
y$PCA1 <- pca1$x[,1]
y$PCA2 <- pca1$x[,2]

ggplot(y) + geom_point(aes(PCA1,PCA2,color=Gate))




# 3D 

library(rgl)

plot3d(x=foo$`Pacific Blue-A`,y=foo$`Pacific Blue-W`,z=foo$`FSC-A`,type = "s",col="lightblue",radius=5e3
       )

mySmoothScatter(foo$`Pacific Blue-A`,foo$`Pacific Blue-W`)


download.file("https://teams.microsoft.com/l/file/C2E7DC57-0FAA-4600-ABA0-0738DB6C24BC?tenantId=966ff61f-fe29-4b94-9886-29307ab0aba7&fileType=xlsx&objectUrl=https%3A%2F%2Figtp.sharepoint.com%2Fsites%2FCancerGeneticsEpigenetics%2FDocumentos%20compartidos%2FMaster%20Maria%20Jose%20Alonso%2FFinal_Documents%2Fareasdate_update26102020.xlsx&baseUrl=https%3A%2F%2Figtp.sharepoint.com%2Fsites%2FCancerGeneticsEpigenetics&serviceName=teams&threadId=19:22d6aed8d78847abb5ab71b027c0553b@thread.skype&groupId=644fe3ed-f37d-4111-a6f9-7fa4d5a75299","output/test.xlsx")
gdata::read.xls("output/test.xlsx")
