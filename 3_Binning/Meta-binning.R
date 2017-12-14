###############################################################################
# Meta-binning #
# A R-script to binning metagenomic sequence based KDE and normal mixture model
###############################################################################

# This is a professional version with some modifications by liqi (liqi at ihb.ac.cn).
# The original author of the script is minleiR. An operable and detailed step-by-step guide is provided, please see http://mingleir.github.io/meta-binning/?from=singlemessage&isappinstalled=0


setwd ("/home/your_workspace")

# A inputfile 'contigs.info' is needed when using this program. It including 4 column, that is contig name, contig length, GC content and coverage

ctgs <- read.table("contig.info", header = T, as.is= T)
ctgs <- subset(ctgs, (length >= 1000) & (coverage >= 3))
ctgs$coverage <- sqrt(ctgs$coverage)
logit <- function(x) 0.5 * log( x / (1-x) )
ctgs$gc <- logit(ctgs$gc/100)
temp.ctgs <- ctgs

my_weight_fun <- function(df){
  temp <- as.numeric(as.matrix(df[, c("coverage", "gc")]))
  temp2 <- rep(temp, times= rep(df[,"length"], 2L))
  dim(temp2) <- c(sum(df[, "length"]), 2L)
  return(temp2)
}

points <- my_weight_fun(ctgs)
u.cov <- quantile(points[,1], 0.999)
points <- points[points[,1] <= u.cov, ]
ctgs <- ctgs[ctgs$coverage <= u.cov,]

library(RColorBrewer)
pal.col <-c("lightgrey", brewer.pal(4,"Set1"), brewer.pal(8,"Dark2"), brewer.pal(12,"Paired"))


# Compute the KDE using the weighted data and create contour graph
library(KernSmooth)
bw <- c(dpik(ctgs$coverage), dpik(ctgs$gc))
dens <- bkde2D(points[, 1:2], bandwidth = bw, gridsize = c(401L, 401L), truncate = TRUE)
contour(dens$x1, dens$x2, dens$fhat, nlevels= 20, xlim = c(0,30), main="name",family="mono")

# Determine the initial group of data based on different level values of contour
# Try some level values manually, and determine the optimal one [group]
library(sp)
ctgs$groupnew <- 0

levelValue <- 0.50
cl <- contourLines(dens$x1, dens$x2, dens$fhat, levels = levelValue)
ctgs$group <- 0
for (i in 1:length(cl)) {
  gi <- as.logical(point.in.polygon(ctgs$coverage, ctgs$gc, cl[[i]]$x, cl[[i]]$y))
  ctgs[gi, "group"] <- i
}

table(ctgs$group)
tapply(ctgs$length, ctgs$group, sum)/10^6
contour(dens$x1, dens$x2, dens$fhat, levels = levelValue, xlim = c(0,30), main="name",family="mono")
palette(pal.col)
text(ctgs$coverage, ctgs$gc, labels = ctgs$group, col= ctgs$group+1, cex = 0.5)

ctgs[ctgs$group==2, ]$groupnew <- 1


levelValue <- 0.80
cl <- contourLines(dens$x1, dens$x2, dens$fhat, levels = levelValue)
ctgs$group <- 0
for (i in 1:length(cl)) {
  gi <- as.logical(point.in.polygon(ctgs$coverage, ctgs$gc, cl[[i]]$x, cl[[i]]$y))
  ctgs[gi, "group"] <- i
}

table(ctgs$group)
tapply(ctgs$length, ctgs$group, sum)/10^6
contour(dens$x1, dens$x2, dens$fhat, levels = levelValue, xlim = c(0,30), main="name",family="mono")
palette(pal.col)
text(ctgs$coverage, ctgs$gc, labels = ctgs$group, col= ctgs$group+1, cex = 0.5)

ctgs[ctgs$group==3, ]$groupnew <- 2
ctgs[ctgs$group==4, ]$groupnew <- 5


levelValue <- 1.00
cl <- contourLines(dens$x1, dens$x2, dens$fhat, levels = levelValue)
ctgs$group <- 0
for (i in 1:length(cl)) {
  gi <- as.logical(point.in.polygon(ctgs$coverage, ctgs$gc, cl[[i]]$x, cl[[i]]$y))
  ctgs[gi, "group"] <- i
}

table(ctgs$group)
tapply(ctgs$length, ctgs$group, sum)/10^6
contour(dens$x1, dens$x2, dens$fhat, levels = levelValue, xlim = c(0,30), main="name",family="mono")
palette(pal.col)
text(ctgs$coverage, ctgs$gc, labels = ctgs$group, col= ctgs$group+1, cex = 0.5)

ctgs[ctgs$group==3, ]$groupnew <- 3


levelValue <- 1.50 
cl <- contourLines(dens$x1, dens$x2, dens$fhat, levels = levelValue)
ctgs$group <- 0
for (i in 1:length(cl)) {
  gi <- as.logical(point.in.polygon(ctgs$coverage, ctgs$gc, cl[[i]]$x, cl[[i]]$y))
  ctgs[gi, "group"] <- i
}

table(ctgs$group)
tapply(ctgs$length, ctgs$group, sum)/10^6
contour(dens$x1, dens$x2, dens$fhat, levels = levelValue, xlim = c(0,26), main="name",family="mono")
palette(pal.col)
text(ctgs$coverage, ctgs$gc, labels = ctgs$group, col= ctgs$group+1, cex = 0.5)

ctgs[ctgs$group==2, ]$groupnew <- 3
ctgs[ctgs$group==1, ]$groupnew <- 4
ctgs[ctgs$group==6, ]$groupnew <- 6
ctgs[ctgs$group==7, ]$groupnew <- 7


#  Combine more than one level value to get the optimal inital group [groupnew]
table(ctgs$groupnew)
tapply(ctgs$length, ctgs$groupnew, sum)/10^6
contour(dens$x1, dens$x2, dens$fhat, levels = 0.5, xlim = c(0,30), main="name",family="mono")
palette(pal.col)
text(ctgs$coverage, ctgs$gc, labels = ctgs$groupnew, col= ctgs$groupnew+1, cex = 0.5)


# Perform the second binning
library(mclust)
z <- unmap(ctgs$groupnew)
Vinv <- hypvol(ctgs[, c("coverage", "gc")], reciprocal = TRUE)
wt <- ctgs$length / max(ctgs$length)
fit.weighted <- me.weighted("VVV", ctgs[, c("coverage", "gc")], z = z, weights = wt, Vinv = Vinv)
ctgs$class <- map(fit.weighted$z)-1
table(ctgs$class)
tapply(ctgs$length, ctgs$class, sum)/10^6


# Read the alignment info against 107 essential genes and get some fragments in group 0 back
ctgs.ess <- read.table("essential_gene.info", header=FALSE, as.is=TRUE)
colnames(ctgs.ess) <- c('count', 'name')
idx <- intersect(ctgs[ctgs$class==0,'name'], ctgs.ess$name)
ctgs0.ess <- ctgs[match(idx, ctgs$name), ]
ctgs0.ess$count <- ctgs.ess[match(idx, ctgs.ess$name), 'count']
ctgs0.dist <- matrix(0, nrow=nrow(ctgs0.ess),ncol= length(unique(ctgs$class)))

# Set a threshold to pick up the non-clustered fragments containing essential single-copy genes
coor.scale <- 3*(IQR(ctgs$coverage)/IQR(ctgs$gc))
ctgs0.ess <- transform(ctgs0.ess, gc.s= gc*coor.scale)
ctgs <- transform(ctgs, gc.s= gc*coor.scale)
each.group <- split(ctgs, f= ctgs$class)

for(i in seq.int(ctgs0.ess$name)){
  tmp.num <-lapply(each.group, function(idx){
    tmp.cen <- kmeans(idx[, c("coverage", "gc.s")], centers=1, iter.max=100)$centers
    tmp.dist <- dist(matrix(c(tmp.cen, ctgs0.ess$coverage[i], ctgs0.ess$gc.s[i]), nrow=2, byrow=T))
    as.matrix(tmp.dist)[2]
  })
  ctgs0.dist[i, ] <- unlist(tmp.num)
}

ctgs0.dist <- ctgs0.dist[, -1]
distlevel <- 1.8

ctgs0.ess$class.dis <- 0
for(i in 1:nrow(ctgs0.dist)){
  if(min(ctgs0.dist[i,]) < distlevel){
    ctgs0.ess$class.dis[i] <- which.min(ctgs0.dist[i,])
  }
}


palette(adjustcolor(col= pal.col, alpha.f = 1))
contour(dens$x1, dens$x2, dens$fhat, cex.lab= 1.5, levels = levelValue, xlim = c(0,30), main="name",family="mono")
points(ctgs$coverage, ctgs$gc, cex= 0.3, pch= 19, col= ctgs$class+1)

palette(adjustcolor(col= pal.col, alpha.f = 0.3))
points(ctgs0.ess$coverage, ctgs0.ess$gc, cex= sqrt(ctgs0.ess$count), pch=19,col= ctgs0.ess$class.dis+1)

x <- ctgs[, c("gc", "coverage", "groupnew")]
x.cen <- matrix(0, nrow= length(unique(x$groupnew)), 3)
m <- 1

for( i in unique(x$groupnew)){
  if (i == 0) next
  x.s <- x[x$groupnew== i, c("gc", "coverage")]
  x.cen[m, ] <- c(kmeans(as.matrix(x.s), 1)$centers, i)
  text(x.cen[m,2], x.cen[m,1], labels = x.cen[m,3], col= "black", cex = 1.0)
  m <- m + 1;
}

ctgs[match(ctgs0.ess$name, ctgs$name), "class"] <- ctgs0.ess$class.dis
table(ctgs$class)
tapply(ctgs$length, ctgs$class, sum)/10^6


# Display the final binning
library(ggplot2)
ctgssort <- ctgs[order(ctgs$class),]
pal.col <-c("lightgrey", "red", "deepskyblue", "seagreen" ,"purple","coral","yellow","brown","magenta","lightseagreen","blue")
ggplot(ctgssort, aes(coverage, gc, col=as.factor(class), size = length)) + geom_point(alpha=0.8) + 
  labs(x= "Read Coverage (Transformed)",y= "GC Content (Transformed)") + xlim(0,30) +
  scale_color_manual(values=pal.col)


# Split the assembly file in each group and export them out in Fasta format
suppressMessages(library(Biostrings))
ctgsSeq <- readDNAStringSet(file.choose())

eachGroup<- split(ctgs$name, ctgs$class)
sapply(names(eachGroup),function(x){seq<- ctgsSeq[eachGroup[[x]]];
  writeXStringSet(seq, filepath= paste("group_", x, ".fa", sep=''))
})




