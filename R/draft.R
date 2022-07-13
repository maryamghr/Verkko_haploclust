setwd("/Users/gharegha/Documents/haploclust_mount/haploclust/pipeline/")
#setwd('../pipeline/')
.libPaths(c(.libPaths(), 'utils/R-packages/'))

.libPaths(c(.libPaths(), 'utils/R-packages/'))
library(dplyr)
#library(SaaRclust)
library(data.table)
library(Rsamtools)
library(doParallel)
library(matrixStats)
library(assertthat)
library(biovizBase)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(dbscan)
library(igraph)

all.sources <- c("../R/calcProbs.R", "../R/countDirectionalReads.R",
                 "../R/dumped_functions.R", "../R/EMclust.R",
                 "../R/export.R", "../R/findClusterPartners.R",
                 "../R/hardClust.R", "../R/helperFuctions.R",
                 "../R/import.R", "../R/importReads.R",
                 "../R/SaaRclust_evaluation_plots.R",
                 "../R/SaaRclust.R", "../R/timedMessage.R",
                 "../R/utils.R", "../R/wrapper_parallel.R",
                 "../R/wrapper.R", "../R/plotting.R")

sapply(all.sources, function(f) source(file.path('..', 'R', f))) %>% invisible()


inputfolder <- '../../HG002/verkko_test_graph/bwa_ss_unitigs'
outputfolder <- '../../HG002/SaaRclust'
input_type <- 'bam'
sex = 'male'
ref.aln.bam <- '../../HG002/verkko_test_graph/ref_aln/asm.r_utg.bam'
input.alignment.files <- list.files(path=inputfolder, pattern='.mdup.bam$', full.names=TRUE)
num.clusters <- 80
EM.iter <- 2
minLib <- 35
minSScov <- 200
HC.only = FALSE
hardclustMethod='hclust'
hard.theta.estim.method='median'
min.cluster.frequency=0.002
numAlignments = 3000
cellNum <- NULL #cells for testing. It should be NULL
log.scale=TRUE
hardclust.file="../../HG002/verkko_test_graph/SaaRclust/Clusters/hard_clusters.RData"
softclust.file="../../HG002/verkko_test_graph/SaaRclust/Clusters/soft_clusters.RData"
MLclust.file="../../HG002/verkko_test_graph/SaaRclust/Clusters/MLclusters.RData"
ss.clust.file="../../HG002/verkko_test_graph/SaaRclust/Clusters/ss_clusters.data"
clust.pairs.file='../../HG002/verkko_test_graph/SaaRclust/Clusters/clust_partners.txt'
wc.cells.file="../../HG002/verkko_test_graph/SaaRclust/Clusters/wc_cells_clusters.data"
chrom.flag.file="../../HG002/verkko_test_graph/ref_aln/chrom_flag.data"
edge.file <- '../../HG002/verkko_test_graph/asm.bp.r_utg.edges.gfa'
store.chrom.flag=TRUE

ss.bam.dir <- '../../ss_bams/'
ss.bam.suffix <- '_haplotagged.bam'
numCPU=2

clust <- runSaaRclust(inputfolder=inputfolder, outputfolder=outputfolder,
                      input_type=input_type, input.alignment.files=input.alignment.files,
                      sex=sex, num.clusters=num.clusters, EM.iter=EM.iter,
                      numAlignments=numAlignments, hardclust.file=hardclust.file,
                      softclust.file=softclust.file, MLclust.file=MLclust.file,
                      ss.clust.file=ss.clust.file, clust.pairs.file=clust.pairs.file,
                      wc.cells.file=wc.cells.file, ref.aln.bam=ref.aln.bam, numCPU=numCPU)



# args <- c(inputfolder, outputfolder, num.clusters, EM.iter, minLib, minSScov,
#           upperQ, hardclustMethod, hard.theta.estim.method, numAlignments, log.scale,
#           hardclust.file, softclust.file, ref.aln.bam, numCPU)
# 
# 
# clust <- runSaaRclust(inputfolder=args[1], outputfolder=args[2], num.clusters=as.numeric(args[3]),
#                       EM.iter=args[4], minLib=as.numeric(args[5]), minSScov=as.numeric(args[6]),
#                       upperQ=as.numeric(args[7]), hardclustMethod=args[8],
#                       hard.theta.estim.method=args[9], numAlignments=as.numeric(args[10]),
#                       HC.only=FALSE, log.scale=args[11], hardclust.file=args[12],
#                       softclust.file=args[13], ref.aln.bam=args[14], numCPU=as.numeric(args[15]))


hard.clust <- clust[[1]]
soft.clust <- clust[[2]]

utg.bam.file <- '../../HG002/verkko_test_graph/ref_aln/asm.r_utg.bam'

hard.clust <- addChromFlag(hard.clust, utg.bam.file, 'hard')
soft.clust <- addChromFlag(soft.clust, utg.bam.file, 'soft')

################################################### hard clust
# Note: The wrapper function misses only chromosome X! Check out the parameters
counts.l <- get(load('../../HG002/verkko_test_graph/SaaRclust/RawData/read_selected_counts.RData'))
chrom.flag <- fread(chrom.flag.file)
# counts.selected <- get(load('../../HG002/verkko_test_graph/SaaRclust/RawData/read_selected_counts.RData'))

num.clusters <- 70 # 300 ----- 80: missing clusters chr22_0, chrX_16, chrY_0, chrY_16
hardclustMethod <- 'hclust'
by_chrom <- FALSE # by_chrom=T does not work well
hardClust.ord <- hardClust(counts.l=counts.l, method=hardclustMethod, 
                           num.clusters=num.clusters, 
                           min.cluster.frequency=min.cluster.frequency, 
                           by_chrom=by_chrom, chrom.flag=chrom.flag)
n.clust <- numFoundClusters(hardClust.ord, chrom.flag, sex)
n.chroms <- numFoundChromosomes(hardClust.ord, chrom.flag, sex)

theta.estim <- estimateTheta(counts.l, hard.clust=hardClust.ord, alpha=alpha, 
                             method=hard.theta.estim.method)

#Merge splitted clusters after hard clustering
k <- 46
if (sex == 'male') {k <- 48}
if (add.garbage.cluster) {k <- k+1}
hardClust.ord.merged <- mergeClusters(hard.clust=hardClust.ord, theta.l=theta.estim, k=k)
n <- numFoundClusters(hardClust.ord.merged, chrom.flag, sex)
n2 <- numFoundChromosomes(hardClust.ord.merged, chrom.flag, sex)

# note: chroms usually get missed only in one direction
# clustering first by chrom and then by direction might work better


################################################## input graph
valid.chroms <- paste0('chr',c(1:22, 'X', 'Y'))
gfa.comp <- get_gfa_components(edge.file)
# gfa.comp <- merge(gfa.comp, chrom.flag, by='rname', all.x=T)
clust.dt <- data.table(rname=names(hardClust.ord), clust=hardClust.ord)
gfa.comp <- merge(gfa.comp, clust.dt, by='rname', all.x=T)
gfa.comp <- merge(gfa.comp, chrom.flag, by='rname', all.x=T)

gfa.comp.remove.na.chrom <- gfa.comp[!is.na(chrom)]

# component compositions in each chromosome
ggplot(gfa.comp.remove.na.chrom, 
       aes(x=factor(chrom, levels = valid.chroms)))+
  geom_bar(aes(fill=factor(component)))+
  xlab('chrom')+
  ylab('component')

# chromosome compositions in each component
ggplot(gfa.comp.remove.na.chrom, 
       aes(x=component))+
  geom_bar(aes(fill=factor(chrom, levels = valid.chroms)), colour='black')+
  ylab('chrom')


# clust compositions in each component
ggplot(gfa.comp.remove.na.chrom, 
       aes(x=component))+
  geom_bar(aes(fill=factor(clust)), colour='black')+
  ylab('clust')


mixed.comp <- gfa.comp[!is.na(chrom) & component==8]
mixed.comp[, table(chrom)]


#################################################### creating graph between clust and components
min.clust.frac <- 0.3
gfa.comp[, clust.count:=.N, by=.(clust, component)]
gfa.comp.clust <- gfa.comp[!is.na(clust)]
gfa.comp.clust[, clust.frac:=clust.count/sum(clust.count), by=component]
gfa.comp.clust <- gfa.comp.clust[, head(.SD,1), by=.(clust, component)]
gfa.comp.clust[clust.frac > min.clust.frac]
gfa.comp.clust[, clust:=paste0('c', clust)]

clust.match.graph <- as.matrix(gfa.comp.clust[, .(clust, component)])
c.match.g <- graph_from_edgelist(clust.match.graph, directed = F)
plot(c.match.g, vertex.size=3)
count_components(c.match.g)
plot_graph_components(in.graph = c.match.g)
#####plot(g, layout = layout_in_circle(g), vertex.size=1, vertex.label=NA)

plot_graph_components(g)


########### PCA analysis of counts data
counts.mat <- Reduce(cbind, counts.l)
rownames(counts.mat) <- rownames(counts.l[[1]])
pca <- prcomp(counts.mat)
pc.dim <- 10
counts.pc <- as.data.table(pca$x[, 1:pc.dim])
counts.pc[, rname:=rownames(counts.mat)]
counts.pc <- merge(counts.pc, chrom.flag, by='rname')

ggplot(counts.pc, aes(x=PC1, y=PC2, colour=chrom=='chr11', shape=factor(flag)))+geom_point()
ggplot(counts.pc[chrom=='chr1'], aes(x=PC1, y=PC2, colour=factor(flag)))+geom_point()

plots <- list()
for (chr in valid.chroms){
  subset.rows <- which(rownames(counts.l[[1]]) %in% chrom.flag[chrom==chr, rname])
  counts.mat.chrom <- counts.mat[subset.rows,]
  rownames(counts.mat.chrom) <- rownames(counts.mat)[subset.rows]
  pca <- prcomp(counts.mat.chrom)
  counts.pc <- as.data.table(pca$x[, 1:2])
  counts.pc[, rname:=rownames(counts.mat.chrom)]
  counts.pc <- merge(counts.pc, chrom.flag, by='rname')
  plt <- ggplot(counts.pc, aes(x=PC1, y=PC2, colour=factor(flag)))+
    geom_point()
  plots[[chr]] <- plt
}
plot_grid(plotlist = plots, ncol = 5, nrow = 5)
####################################################################

# numfoundclusters for hard.clust...
# check random shuffling
# check hard clustering for all reads

 
# counting w/c in bam files
# inputfolder <- '../../HG00733/bwa_ss_unitigs/'
bam.files <- list.files(path=inputfolder, pattern='.mdup.bam$', full.names=TRUE)
counts.l <- count.wc.bam(bam.files, numCPU = 42)

### Perform k-means hard clustering method ###
set.seed(1000) #in order to reproduce hard clustering results
hardClust.ord <- hardClust(counts.l, num.clusters=num.clusters, num.iter = 10)
names(hardClust.ord) <- rownames(counts.l[[1]])

# looking into hard clustering features
num.cells <- ncol(ratios.m)
colnames(ratios.m) <- paste0('cell', 1:num.cells)
ratios.m[, rname:=rownames(counts.l[[1]])]
ratios.m <- merge(ratios.m, chrom.flag, by='rname')
ratios.m[, chrom_flag:=paste0(chrom, '_',  flag), by=.(chrom, flag)]

theta.estim <- estimateTheta(counts.l, hard.clust=hardClust.ord, alpha=alpha)
hardClust.ord.merged <- mergeClusters(hard.clust=hardClust.ord, theta.l=theta.estim, k=47)
#Re-estimate theta parameter after cluster merging
theta.estim <- estimateTheta(counts.l, hard.clust=hardClust.ord.merged, alpha=alpha)

#Initialize theta parameter
theta.param <- theta.estim
#Estimate pi parameter based on # of PB reads in each cluster
readsPerCluts <- table(hardClust.ord.merged)
pi.param <- readsPerCluts/sum(readsPerCluts)

#save hard clustering results into a file
hard.clust <- list(ord=hardClust.ord.merged, unmerged.ord=hardClust.ord, theta.param=theta.param)#, pi.param=pi.param, pb.chr = chr.rows, pb.flag = pb.flag)
destination <- file.path(Clusters.store, outputfilename)

################## eval:
hard.clust <- get(load('../../HG002/verkko_test_graph/SaaRclust/Clusters/hard_clusters.RData'))
soft.clust <- get(load('../../HG002/verkko_test_graph/SaaRclust/Clusters/soft_clusters.RData'))
bam.file <- '../../HG002/verkko_test_graph/ref_aln/asm.r_utg.haplotagged.bam'

#bam.files <- list.files(path=inputfolder, pattern='.mdup.bam$', full.names=TRUE)
#counts.l <- count.wc.bam(bam.files, numCPU = 42)
#save(file='../../HG00733/SaaRclust/read_counts.RData', counts.l)

counts.l <- get(load('../../HG002/verkko_test_graph/SaaRclust/RawData/read_counts.RData'))
#chrom.flag <- getChromFlag(bam.file)
chrom.flag <- fread('../../HG002/verkko_test_graph/SaaRclust/chrom_flags.data')
#hard.clust <- addChromFlag(clust.obj=hard.clust, bam.file=bam.file, clust.type = 'hard')
#soft.clust <- addChromFlag(soft.clust, bam.file=bam.file, clust.type = 'soft')
##valid.chroms <- paste0('chr', c(1:22, "X"))
##expected.clusters <- paste0(rep(valid.chroms, each=2), rep(c('_0','_16'),23))


#numFoundClusters(hard.clust$ord, chrom.flag) 
# hard clust

counts.selected <- get_representative_counts(counts.l=counts.l, num.alignments = 10000, min.ss.cov=200)[[2]]

### hclust
wfrac.matrix <- get.wfrac.matrix(counts.dt=counts.selected)
# hc <- hard.hclust(wfrac.matrix=wfrac.matrix, num.clusters=60, chrom.flag=chrom.flag)
hc <- hard.hclust(counts.l = counts.selected, num.clusters=80, chrom.flag=chrom.flag)

##### assign theta based on median wfrac values
theta.estim <- estimateTheta(counts.selected, hc)

hardClust.ord.merged <- mergeClusters(hard.clust=hc, theta.l=theta.estim, k=47)
clust.to.chrom <- numFoundClusters(ord=hardClust.ord.merged, chrom.flag)
##
## plot heatmap
counts.selected <- merge(counts.selected, chrom.flag, by='rname', all.x=T)
plt.l <- plot_heatmap_wfrac(counts.l=counts.selected, compute_Wfrac=F)
plt.chrom.l <- plot_heatmap_wfrac(counts.l=counts.selected, by_chrom =F)

# num.clusters=47: missing clusters: chr14_0, chr22_0, clust accuracy: 0.9745495
# num.clusters=50: missing clusters: chr22_0, clust accuracy: 0.9877252
# num.clusters=60,80 on 10000 reads yields very good results: 99.8% acc with finding all chrom_dir clsuters, time=18s

clust <- data.table(rname=names(hc), clust=hc)

ord <- clust[, clust]
#theta.estim <- estimateTheta(counts.selected.l, hard.clust=ord, alpha=alpha)



####### theta ...
#theta[, state:=0]
#theta[median.w.frac <= -0.5, state:=-1, by=.(clust,lib)]
#theta[median.w.frac > -0.5 & median.w.frac < 0.5, state:=0, by=.(clust,lib)]
#theta[median.w.frac >= 0.5, state:=1, by=.(clust,lib)]

#theta[, wc.frac:=length(which(state==0))/.N, by=clust]

#wc.frac <- theta[, head(.SD, 1), by=clust]
#setkey(wc.frac, chrom, flag)
#wc.frac

#theta.matrix <- dcast(theta, clust~lib, value.var = "state")
#rownames(theta.matrix) <- theta.matrix[, clust]
#theta.matrix[, clust:=NULL]
#Heatmap(theta.matrix, show_row_names = T, show_column_names = F)

#theta.l <- split(theta[, .(lib, state)], by='lib')
#lapply(theta.l, function(dt) dt[, lib:=NULL]) %>% invisible()

# better than the probabilistic approach--- test it
# TODO: try removing the clusters with high variance
#################

set.seed(1000)
num.clusters <- 60
km = hardClust(counts.selected.l, num.clusters= num.clusters) #, iter.max = 100)
km.chrom = hardClust(counts.selected.l, num.clusters= 40, by_chrom=T)
clust <- data.table(rname=names(km), clust=km)
clust <- merge(clust, chrom.flag, by='rname', all.x=T)
clust.by.chrom <- data.table(rname=names(km.chrom), clust=km.chrom)
clust.by.chrom <- merge(clust.by.chrom, chrom.flag, by='rname', all.x=T)
# num.clusters=50: 7 missing clusters: chr10_0, chr13_0, chr15_16, chr16_0, chr18_0, chr21_0, chr21_16, clust accuracy: 0.9027575
# num.clusters=60: missing clusters: chr21_16, chr22_16, clust accuracy: 0.9834553
# num.clusters=80 on 10000 reads yields very good results: 99.5% acc with finding all chrom_dir clsuters, time=18s

##clust.to.chrom <- numFoundClusters(ord=clust$clust, chrom.flag=chrom.flag) # ord should be named
##numFoundChromosomes(ord=clust.by.chrom$clust, chrom.flag=chrom.flag)

theta.estim <- estimateTheta(counts.selected.l, hard.clust=ord, alpha=alpha) # it has warnings

# test clust=1 (chr4_0)
test <- counts.selected[rname %in% clust[ord==1, rname]]


# theta.wc <- Reduce(c, lapply(theta.estim, function(x) x[,3]))
# table(theta.wc)
# test estimateTheta...
# problem: theta_wc is often the most likely state
# Looking at single-cell counts for cluster 1 (chr17_16) to see strand states of single-cells...
# split(counts.selected.l[[j]], hard.clust)[[1]] for j=1:115 ....
# The data (strand states) does not look clean although the cluster is rather clean containing mainly the true reads!
# double check the unitigs and ss reads mapping orientations (in IGV)
lib1.clust3 <- counts.selected.l[["HG00733.P00IL002"]][which(ord==3)] # only chr6_0
lib1.clust3[, rname:=rownames(counts.selected.l[["HG00733.P00IL002"]])[which(ord==3)]]
chrom.flag[rname %in% lib1.clust3[, rname]] # shows a mixture of chromosomes!!!
#lib1.clust3[, rname] != names(ord)[which(ord==3)]!!!

hardClust.ord.merged <- mergeClusters(hard.clust=ord, theta.l=theta.estim, k=47)
clust.to.chrom <- numFoundClusters(ord=hardClust.ord.merged, chrom.flag=chrom.flag)

theta.estim <- estimateTheta(counts.selected.l, clust.to.chrom, hard.clust=hardClust.ord.merged, alpha=alpha)
compare_theta_with_wfrac(theta=theta.estim, clust.to.chrom, countfiles.dir)
### Note: almost all thetas are wc... I guess kmeans does not work well with sparse data
# cluster only by chrom
abs.ratios <- abs(ratios.m)
abs.ratios[abs.ratios==0]=0.5
chrom.km = hardClust(counts.selected.l, num.clusters= 30)
clust <- data.table(rname=names(chrom.km), ord=chrom.km)
clust <- merge(clust, chrom.flag, by='rname', all.x=T)
#clust <- clust[chrom %in% valid.chroms]
ord = chrom.km

# ord should be named
##numFoundChromosomes(ord, chr=clust$chrom, expected.clusters = valid.chroms)
# cluster by chrom


chr.hists <- list()

for (data.file in list(hard.clust, soft.clust)){
  if ("ord" %in% names(data.file)) {utg.names <- names(data.file$ord)} else{utg.names <- rownames(data.file$soft.pVal)}
  gr.truth <- data.table(name=utg.names, chrom=data.file$chrom, flag=data.file$flag)
  gr.truth <- gr.truth[chrom %in% valid.chroms]
  gr.truth[, chrom_flag:=paste0(chrom, '_', flag), by=.(chrom, flag)]
  # order
  gr.truth[, chr.ord:=gsub('chr','',chrom), by=chrom]
  gr.truth[chrom=='chrX', chr.ord:="23"]
  gr.truth[, chr.ord:=as.numeric(chr.ord)]
  setkey(gr.truth, chr.ord)
  #######
  chr.hists[[1+length(chr.hists)]] <- ggplot(gr.truth)+geom_bar(aes(x=chrom_flag))
}

chr.flag.ord <- paste0(rep(valid.chroms, each=2), rep(c('_0','_16'),23))

## contibait
# contibait: http://bioconductor.org/packages/release/bioc/manuals/contiBAIT/man/contiBAIT.pdf
library(contiBAIT)
strand.freq <- StrandFreqMatrix(counts=as.matrix(ratios.m))
strand.states <- preprocessStrandTable(exampleStrandFreq, lowQualThreshold=0.8)
cl <- clusterContigs(strand.states$strandMatrix)
data("exampleDividedChr")
barplotLinkageGroupCalls(cl, exampleDividedChr)

# todo: PCA dim reduction:
#https://www.youtube.com/watch?v=K5HcVuLra38
pca.ratios <- prcomp(ratios.m)
summary(pca.ratios)
pca <- as.data.table(pca.ratios$x)
pca[,rname:=rownames(counts.l[[1]])]
pca = merge(pca, chrom.flag, by='rname')
pca[,chrom_flag:=paste0(chrom, '_', flag)]
ggplot(pca, aes(x=PC1, y=PC2, color=paste(chrom, flag)))+geom_point()


# plotting binned w/c read counts for large unitigs
unitigs.len <- ...
# importing bam as Granges
counts.utg000445l <- importBams(bamfolder=inputfolder, chromosomes='utg000445l', bin.length=10000)

################################################ plotting
# hard clust heatmap

clust.chrom.count <- data.table(chrom=hard.clust$chrom, clust=hard.clust$ord)[!is.na(chrom)]
clust.chrom.count[, count:=.N, by=.(chrom, clust)]
clust.chrom.count <- clust.chrom.count[, head(.SD, 1), by=.(chrom, clust)]
mat <- dcast(clust.chrom.count, chrom~clust, value.var = 'count', fill=0)
mat <- mat[order(factor(chrom, levels=paste0('chr', c(1:22, 'X'))))]
row.names <- mat[, chrom]
mat[, chrom:=NULL]
mat <- as.matrix(mat)
rownames(mat) <- row.names
heatmap(mat, scale = 'column', Rowv = NA, Colv = NA, col= colorRampPalette(brewer.pal(8, "Blues"))(25))

# input data insights- bar plot of true clusters in all and selected unitigs
ggplot(chrom.flag)+ geom_bar(aes())

# complex Heatmap:
#Heatmap(mat)
################### SS Wfrac in long reads/unitigs
counts <- data.table()
for (i in 1:length(counts.l)){
  counts.dt <- ratios.l[[i]] # counts.l[[i]]
  lib.name=names(counts.l)[i]
  rnames=rownames(counts.l[[i]]) #(counts.dt)
  
  counts.dt[, `:=`(rname=rnames, lib.name=lib.name)]
  
  counts <- rbind(counts, counts.dt)
}

colnames(counts)[1] <- 'W.frac' #counts[, W.frac:=(w-c)/(w+c)]
counts.selected <- counts

rnames <- counts[, unique(rname)]
lib.names <- counts[, unique(lib.name)]
expand.counts <- data.table(expand.grid(rnames, lib.names))
colnames(expand.counts) <- c('rname', 'lib.name')
counts <- merge(counts, expand.counts, by=c('rname', 'lib.name'), all=T)
counts[is.na(counts)] <- 0

counts <- merge(counts, chrom.flag, by='rname')

# subsetting the highest coverage unitigs
counts[, ss.cov:=sum(w+c), by=rname]
setkey(counts, ss.cov)

counts.selected <- counts[, head(.SD, 1), by=rname]
counts.selected <- counts.selected[ss.cov >= min.ss.cov]

selected.rname <- counts.selected[, rname]
if (nrow(counts.selected) > num.alignments){
  selected.rname <- tail(selected.rname, num.alignments)
}
counts.selected <- counts[rname %in% selected.rname]
counts.selected[, avg.w.frac:=mean(W.frac[W.frac!=0]), by=.(lib.name,chrom,flag)]

avg.unitigs.w.frac <- counts.selected[, head(.SD,1), by=.(lib.name,chrom,flag)][, .(lib=lib.name, chrom, flag, W.frac=avg.w.frac)]
avg.unitigs.w.frac[, chrom:=paste0(chrom,'_', flag, '_unitig'), by=.(chrom, flag)]
avg.unitigs.w.frac[, flag:=NULL]

w.frac = rbind()