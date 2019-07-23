library("Seurat")
library("SC3")
library("SingleCellExperiment")
library("mclust")
library("Rtsne")
library("kernlab")
library("ConsensusClusterPlus")
library("class")
library("dbscan")
library("pals")
library("data.table")


Path<-"/home/nikolay/WABI/K_Pietras/Easy_scRNAseq_tSNE_Cluster/Large_DataSets/"
setwd(Path)
#files<-list.files(pattern="*.txt")
files<-c("montoro.txt","muraro.txt","nestorowa.txt")

for(i in 1:length(files))
{
df <- files[i]
file_name <- tail(unlist(strsplit(df,"/")),1)
pdf(paste0("Easy_Output/COMPARE_CLUSTERS_PLOTS_",file_name,".pdf"), paper="a4r", width = 210, height = 297)
sink(paste0("Easy_Output/COMPARE_CLUSTERS_LOG_",file_name),split=TRUE)

# READ DATA AND MAKE TSNE PLOT
print(paste0("READING DATA AND SELECTING OPTIMAL TSNE DIMENSIONALITY REDUCTION: ",file_name))
system(paste0("grep 'OPTIMAL NUMBER OF PRINCIPAL COMPONENTS' Easy_Output/LOG_",file_name," > optPC.txt"))
optPC<-as.numeric(unlist(strsplit(unlist(strsplit(scan("optPC.txt",what="character")[2]," = "))[2],","))[1])
system(paste0("grep 'OPTIMAL MIN SIZE OF CLUSTERS' Easy_Output/LOG_",file_name," > minPts.txt"))
minPts<-as.numeric(unlist(strsplit(unlist(strsplit(scan("minPts.txt",what="character")[2]," = "))[2],","))[1])
system("rm optPC.txt")
system("rm minPts.txt")

expr <- suppressWarnings(as.data.frame(fread(df,sep="\t")))
rownames(expr)<-expr$V1; expr$V1<-NULL;
expr <- expr[grepl("ERCC_",rownames(expr))==FALSE,]
expr <- expr[rowMeans(as.matrix(expr)) >= 1,]
expr<-na.omit(expr)
print(expr[1:5,1:5])
N_cells <- dim(expr)[2]
optPerp <- round(sqrt(N_cells),0)

N_tsne <- 100
tsne_out <- list(length = N_tsne)
KL <- vector(length = N_tsne)
for(k in 1:N_tsne)
{
  tsne_out[[k]]<-Rtsne(t(log10(expr+1)),initial_dims=optPC,verbose=FALSE,check_duplicates=FALSE,
                       perplexity=optPerp,dims=2,max_iter=10000)
  KL[k]<-tail(tsne_out[[k]]$itercosts,1)
  print(paste0("FINISHED ",k," TSNE ITERATION"))
}
names(KL) <- c(1:N_tsne)
opt_tsne <- tsne_out[[as.numeric(names(KL)[KL==min(KL)])]]$Y


par(mfrow=c(3,3))

# 1) HDBSCAN
print("PERFORMING HDBSCAN CLUSTERING")
print("**********************************************************************************************")
res_opt <- hdbscan(opt_tsne, minPts = minPts)
print(res_opt)
print("USING KNN FOR CLASSIFYING HDBSCAN OUTLIERS")
if(length(res_opt$cluster[res_opt$cluster==0]) > 0 & length(res_opt$cluster[res_opt$cluster==1]) > 0)
{
  res_opt$cluster[res_opt$cluster==0]<-class::knn(opt_tsne[which(res_opt$cluster!=0),],
                                                  opt_tsne[which(res_opt$cluster==0),],
                                                  res_opt$cluster[res_opt$cluster!=0], k=3)
}
print(res_opt)
N_clust<-length(sort(unique(res_opt$cluster)))
if(N_clust <= 25){colors <- cols25(N_clust)}else{colors <- rainbow(N_clust)}
names(colors) <- sort(unique(res_opt$cluster))
plot(opt_tsne, main="", col=colors[as.character(res_opt$cluster)], xlab="tSNE1", ylab="tSNE2",cex=0.5,pch=19)
mtext("HDBSCAN")

# 2) KMEANS ON TSNE
print("PERFORMING KMEANS CLUSTERING")
print("**********************************************************************************************")
kmeans_cluster <- kmeans(opt_tsne, N_clust)
print(table(kmeans_cluster$cluster))
plot(opt_tsne,main="", col=colors[as.character(kmeans_cluster$cluster)], xlab="tSNE1", ylab="tSNE2")
mtext("KMEANS")

# 3) GAUSSIAN MIXTURE MODEL ON TSNE
print("PERFORMING GAUSSIAN MIXTURE MODEL CLUSTERING")
print("**********************************************************************************************")
GMM_cluster <- Mclust(opt_tsne, N_clust)
summary(GMM_cluster)
plot(opt_tsne,main="", col=colors[as.character(GMM_cluster$classification)], xlab="tSNE1", ylab="tSNE2")
mtext("GAUSSIAN MIXTURE MODEL")

# 4) HIERARCHICAL CLUSTERING
print("PERFORMING HIERARCHICAL CLUSTERING")
print("**********************************************************************************************")
dist.euc<-dist(opt_tsne)
hcl.euc<-hclust(dist.euc, method="ward.D2")
clusters.euc<-cutree(hcl.euc, N_clust)
print(table(clusters.euc))
plot(opt_tsne,main="", col=colors[as.character(clusters.euc)], xlab="tSNE1", ylab="tSNE2")
mtext("HIERARCHICAL CLUSTERING")

# 5) SPECTRAL CLUSTERING
print("PERFORMING SPECTRAL CLUSTERING")
print("**********************************************************************************************")
sc <- specc(opt_tsne, centers=N_clust)
print(table(as.numeric(sc)))
plot(opt_tsne,main="", col=colors[as.character(as.numeric(sc))], xlab="tSNE1", ylab="tSNE2")
mtext("SPECTRAL CLUSTERING")

# 6) SC3
print("PERFORMING SC3 CLUSTERING")
print("**********************************************************************************************")
sce <- SingleCellExperiment(assays = list(counts = t(opt_tsne + 50)))
rowData(sce)$feature_symbol <- c("tSNE1","tSNE2")
logcounts(sce) <- t(log(opt_tsne + 50))
sce <- sc3(sce, ks = N_clust, gene_filter=FALSE)
print(table(as.numeric(colData(sce)[,1])))
plot(opt_tsne, main="", col=colors[as.character(as.numeric(colData(sce)[,1]))], xlab="tSNE1", ylab="tSNE2")
mtext("SC3")

# 7) BOOTSTRAP CONSENSUS CLUSTERING
print("PERFORMING BOOTSTRAP CONSENSUS CLUSTERING")
print("**********************************************************************************************")
maxK=N_clust
reps=200
results<-ConsensusClusterPlus(as.matrix(t(opt_tsne)), maxK=maxK, reps=reps, pItem=0.8, pFeature=1, 
                              title="Consensus Clustering",clusterAlg="km", distance="euclidean", seed=123, plot="png")
plot(opt_tsne, main="", col=colors[as.character(results[[N_clust]][["consensusClass"]])], xlab="tSNE1", ylab="tSNE2")
#Consensus_Matrix<-results[[N_clust]][["consensusMatrix"]]
#dist_cc.euc<-dist(Consensus_Matrix)
#hcl_cc.euc<-hclust(dist_cc.euc, method="ward.D2")
#clusters_cc.euc<-cutree(hcl_cc.euc, N_clust)
#plot(opt_tsne, main="", col=colors[as.character(clusters_cc.euc)], xlab="tSNE1", ylab="tSNE2")
print(table(results[[N_clust]][["consensusClass"]]))
mtext("BOOTSTRAP CONSENSUS CLUSTERING")

# 8) SNN-CLIQ
print("PERFORMING SNN-CLIQ CLUSTERING")
print("**********************************************************************************************")
distan <- "euclidean"
par.k <- optPerp
par.r <- 0.5
par.m <- 0.3
source("/home/nikolay/WABI/K_Pietras/SNN.R")
SNN(data = opt_tsne, outfile = "snn-cliq.txt", k = par.k, distance = distan)
snn.res <- system(paste0("python /home/nikolay/WABI/K_Pietras/./Cliq.py -i snn-cliq.txt -o res-snn-cliq.txt -r", 
                       par.r," -m", par.m), intern=TRUE)
cat(paste(snn.res, collapse = "\n"))
snn.res <- read.table("res-snn-cliq.txt")
print(table(as.numeric(snn.res[,1])))
N_clust_snn_cliq <- length(unique(as.numeric(snn.res[,1])))
if(N_clust_snn_cliq <= 25){colors_snn_cliq <- cols25(N_clust_snn_cliq)}else{colors_snn_cliq <- rainbow(N_clust_snn_cliq)}
names(colors_snn_cliq)<-seq(from=1, by=1, to=N_clust_snn_cliq)
plot(opt_tsne, main="", col=colors_snn_cliq[as.character(as.numeric(snn.res[,1]))], xlab="tSNE1", ylab="tSNE2")
mtext("SNN-CLIQ")
system("rm snn-cliq.txt res-snn-cliq.txt")

# 9) SEURAT
print("PERFORMING SEURAT CLUSTERING")
print("**********************************************************************************************")
data <- opt_tsne + 50
rownames(data) <- paste0("CELL",c(1:dim(data)[1]))
colnames(data) <- paste0("GENE",c(1:dim(data)[2]))
#seurat_obj <- CreateSeuratObject(counts = t(data), project = "SEURAT")
seurat_obj <- CreateSeuratObject(raw.data = t(data), project = "SEURAT")
seurat_obj
#seurat_obj <- FindNeighbors(object = seurat_obj)
seurat_obj<-FindClusters(object = seurat_obj, genes.use=colnames(data), resolution=0.1, print.output=0, save.SNN=TRUE)
#seurat_obj<-FindClusters(seurat_obj, resolution=0.1, save.SNN=TRUE)
PrintFindClustersParams(object = seurat_obj)
N_clust_seurat<-length(sort(unique(as.numeric(seurat_obj@ident))))
if(N_clust_seurat <= 25){colors_seurat <- cols25(N_clust_seurat)}else{colors_seurat <- rainbow(N_clust_seurat)}
names(colors_seurat) <- sort(unique(as.numeric(seurat_obj@ident)))
print(table(as.numeric(seurat_obj@ident)))
plot(opt_tsne, main="", col=colors_seurat[as.character(as.numeric(seurat_obj@ident))], xlab="tSNE1", ylab="tSNE2")
mtext("SEURAT")

title(paste0(file_name,", N = ", N_cells,", optPC = ",optPC,", optPerp = ",optPerp,", N_clust = ", N_clust),
      line=-2,outer=TRUE,cex=2)
print("**********************************************************************************************")
print("**********************************************************************************************")
print("**********************************************************************************************")

sink()
dev.off()
}
