#CELLS ARE COLUMNS, GENES ARE ROWS, INPUT DATA IS TAB-DELIMITED
easy_scrnaseq_tsne_cluster<-function(df)
{
  library("data.table")
  library("Rtsne")
  library("matrixStats")
  library("dbscan")
  library("class")
  library("pals")
  library("foreach")
  library("doParallel")
  library("parallel")
  numCores<-detectCores()
  registerDoParallel(numCores)
  
  ################################################## READING DATA ##################################################
  file_name <- tail(unlist(strsplit(df,"/")),1)
  pdf(paste0("Easy_Output/PLOTS_",file_name,".pdf"), paper="a4r", width = 210, height = 297)
  sink(paste0("Easy_Output/LOG_",file_name),split=TRUE)
  print(paste0("START WORKING WITH FILE ", file_name))
  print("READING AND FILTERING DATA")
  expr <- suppressWarnings(as.data.frame(fread(df,sep="\t")))
  rownames(expr)<-expr$V1; expr$V1<-NULL;
  expr <- expr[grepl("ERCC_",rownames(expr))==FALSE,]
  expr <- expr[rowMeans(as.matrix(expr)) >= 1,]
  expr<-na.omit(expr)
  N_cells<-dim(expr)[2]
  print(expr[1:5,1:5])
  print(paste0("DATASET CONTAINS ",dim(expr)[1]," GENES AND ",N_cells," CELLS"))
  print(paste0("PERCENT OF MISSING VALUES = ",round(100*sum(expr==0)/(dim(expr)[1]*dim(expr)[2]),0)," %"))
  
  ######################################## SELECTING HIGHLY VARIABLE GENES ########################################
  print("CHECKING THE NUMBER OF HIGHLY VARIABLE GENES")
  mean_expr<-as.numeric(rowMeans(expr,na.rm=TRUE))
  log10_mean_expr<-log10(mean_expr)
  names(log10_mean_expr)<-rownames(expr)
  sd_expr<-rowSds(as.matrix(expr),na.rm=TRUE)
  cv_squared_expr<-(sd_expr/mean_expr)^2
  log10_cv_squared_expr<-log10(cv_squared_expr)
  names(log10_cv_squared_expr)<-rownames(expr)
  plot(log10_cv_squared_expr~log10_mean_expr,pch=20,cex=0.5,xlab="LOG10 ( MEAN EXPRESSION )",ylab="LOG10 ( CV2 )")
  
  N_bin<-100
  bin_size<-(max(log10_mean_expr)-min(log10_mean_expr))/N_bin
  mean_cv_squared_points<-vector(length=N_bin)
  sd_cv_squared_points<-vector(length=N_bin)
  mean_expr_points<-vector(length=N_bin)
  a<-min(log10_mean_expr)
  for(i in 1:N_bin)
  {
    b<-a+bin_size
    mean_cv_squared_points[i]<-mean(log10_cv_squared_expr[log10_mean_expr>=a & log10_mean_expr<=b],na.rm=TRUE)
    sd_cv_squared_points[i]<-sd(log10_cv_squared_expr[log10_mean_expr>=a & log10_mean_expr<=b],na.rm=TRUE)
    mean_expr_points[i]<-mean(log10_mean_expr[log10_mean_expr>=a & log10_mean_expr<=b],na.rm=TRUE)
    a<-b
  }
  points(mean_cv_squared_points ~ mean_expr_points,col="red",pch=20)
  boundary_cv_squared_points<-mean_cv_squared_points + 2*sd_cv_squared_points
  
  fit_mean_cv_squared_points <- loess(mean_cv_squared_points[is.finite(mean_cv_squared_points)&
                                                               is.finite(mean_expr_points)]~
                                      mean_expr_points[is.finite(mean_cv_squared_points)&
                                                         is.finite(mean_expr_points)],span=1)
  fit_boundary_cv_squared_points <- loess(boundary_cv_squared_points[is.finite(boundary_cv_squared_points)&
                                                                       is.finite(mean_expr_points)]~
                                        mean_expr_points[is.finite(boundary_cv_squared_points)
                                                         &is.finite(mean_expr_points)],span=1)
  pred_cv_squared_expr <- predict(fit_mean_cv_squared_points, log10_mean_expr)
  pred_boundary_cv_squared_expr <- predict(fit_boundary_cv_squared_points, log10_mean_expr)
  
  expr_pruned <- expr[log10_cv_squared_expr >= pred_boundary_cv_squared_expr,]
  expr_pruned <- expr_pruned[grepl("NA",rownames(expr_pruned))==FALSE,]
  text(log10_mean_expr[rownames(expr_pruned)],log10_cv_squared_expr[rownames(expr_pruned)],
       labels=rownames(expr_pruned),cex=0.5)
  lines(fit_mean_cv_squared_points$fitted~
          mean_expr_points[is.finite(mean_expr_points)],col="red",lwd=3)
  mtext(paste0("NUMBER OF HIGHLY VARIABLE GENES = ",dim(expr_pruned)[1]))
  print(paste0("AFTER SELECTEING HVG DATASET CONTAINS ",dim(expr)[1]," GENES AND ",N_cells," CELLS"))
  
  ################################### SELECTING SIGNIFICANT PRINCIPAL COMPONENTS ###################################
  print("PERFORMING PRINCIPAL COMPONENT ANALYSIS")
  PC <- prcomp(t(log10(expr + 1)), center=TRUE, scale=FALSE)
  expl_var <- PC$sdev^2/sum(PC$sdev^2)
  barplot(expl_var[1:50],ylab="EXPLAINED VARIANCE",main="VARIANCE EXPLAINED BY PRINCIPAL COMPONENTS",
          names.arg=paste0("PC",seq(1:50)),col="darkgreen")
  plot(PC$x[,1:2], col="blue", main=paste0("PCA PLOT: ",file_name))
  mtext(paste0("PERCENT OF MISSING VALUES = ",round(100*sum(expr==0)/(dim(expr)[1]*dim(expr)[2]),0)," %"))
  
  print("SELECTING SIGNIFICANT PRINCIPAL COMPONENTS")
  if(N_cells <= 90){N_perm <- 50}else{N_perm <- 10}
  expl_var_perm <- matrix(NA, ncol = length(PC$sdev), nrow = N_perm)
  for(k in 1:N_perm)
  {
    expr_perm <- apply(expr,2,sample)
    PC_perm <- prcomp(t(log10(expr_perm+1)), center=TRUE, scale=FALSE)
    expl_var_perm[k,] <- PC_perm$sdev^2/sum(PC_perm$sdev^2)
    print(paste0("FINISHED ",k," PERMUTATIONS"))
  }
  plot(expl_var[1:50]~seq(1:50),ylab="EXPLAINED VARIANCE",main="VARIANCE EXPLAINED BY PRINCIPAL COMPONENTS",
          col="darkgreen",type='o',xlab="PRINCIPAL COMPONENTS")
  lines(colMeans(expl_var_perm)[1:50]~seq(1:50),col="red")
  legend("topright",c("Explained by PCS","Explained by chance"),fill=c("darkgreen","red"),inset=0.02)
  pval <- apply(t(expl_var_perm) >= expl_var,1,sum) / N_perm
  plot(pval[1:50]~seq(1:50),col="darkred",type='o',xlab="PRINCIPAL COMPONENTS",ylab="PVALUE",
       main="SIGNIFICANCE OF PRINCIPAL COMPONENTS")
  optPC<-head(which(pval>=0.05),1)-1
  if(optPC>30){optPC<-30}else{optPC<-optPC}
  mtext(paste0("OPTIMAL NUMBER OF PRINCIPAL COMPONENTS = ",optPC))
  print(paste0("OPTIMAL NUMBER OF PRINCIPAL COMPONENTS = ",optPC,
               ", THEY TOGETHER EXPLAIN ",round(sum(expl_var[1:optPC])*100,0),"% OF VARIANCE"))
  
  ######################################### SELECTING OPTIMAL PERPLEXITY ##########################################
  print("SELECTING OPTIMAL PERPLEXITY")
  optPerp <- round(sqrt(N_cells),0)
  print(paste0("OPTIMAL PERPLEXITY = ",optPerp))
  tsne_opt_perp <- Rtsne(t(log10(expr+1)),initial_dims=optPC,verbose=FALSE,check_duplicates=FALSE,
                     perplexity=optPerp,dims=2,max_iter=10000)
  plot(tsne_opt_perp$Y,main=paste0("TSNE PLOT: ",file_name),col="blue",xlab="tSNE1",ylab="tSNE2")
  mtext(paste0("OPTIMAL PERPLEXITY = ",optPerp))
  print("SCANNING RANGE OF PERPLEXITIES")
  par(mfrow=c(3,3))
  perp_range<-vector(length=9)
  perp_range[1]<-3
  if(N_cells/3 <= 5000){perp_range[9] <- N_cells/3-1}else{perp_range[9] <- 5000}
  perp_step<-(optPerp-3)/4
  for(s in 2:8){perp_range[s]<-perp_range[s-1]+perp_step}
  print(round(perp_range,0))
  for(j in round(perp_range,0))
  {
    tsne_perp_iter<-Rtsne(t(log10(expr+1)),initial_dims=optPC,verbose=FALSE,check_duplicates=FALSE,
                     perplexity=j,dims=2,max_iter=10000)
    plot(tsne_perp_iter$Y,col="blue",xlab="tSNE1",ylab="tSNE2",cex=0.5)
    mtext(paste0("perplexity = ",j))
    print(paste0("FINISHED PERPLEXITY ",j))
  }
  title(paste0(file_name,", N = ",N_cells,", optPC = ",optPC,", optPerp = ",optPerp),line=-2,outer=TRUE,cex=2)
  par(mfrow=c(1,1))
  
  #################################### PERFORMING HDBSCAN HYPERPARAMETER TUNING #####################################
  print("PERFORMING HDBSCAN HYPERPARAMETER TUNING")
  if(N_cells <= 90){N_iter <- 120}else{N_iter <- 60}
  if(N_cells <= 200){N_pt <- 10}else{N_pt <- 50}
  score <- vector(length = N_pt-2)
  for(i in 3:N_pt)
  {
    print(paste0("MIN POINTS = ",i))
    score_iter <- foreach (1:N_iter, .combine = c) %dopar% 
    {
      tsne_iter <- Rtsne(t(log10(expr+1)),initial_dims=optPC,verbose=FALSE,check_duplicates=FALSE,
                         perplexity=optPerp,dims=2,max_iter=10000)
      res <- hdbscan(tsne_iter$Y, minPts = i)
      score_iter_temp <- sum(res$membership_prob < 0.05) / N_cells
      #score_iter_temp <- sum(res$membership_prob < 0.05) - length(unique(res$cluster))
      #if(N_cells<=500){score_iter_temp<- - 2*log(mean(res$membership_prob,na.rm=TRUE))}else{score_iter_temp<- 
      #      length(unique(res$cluster))/log(N_cells) - 2*log(mean(res$membership_prob,na.rm=TRUE))}
      return(score_iter_temp)
    }
    score[i-2] <- mean(score_iter, na.rm = TRUE)
  }
  plot(log(score + 1) ~ seq(from=3, to=N_pt, by=1), type='o', xlab="MIN SIZE OF CLUSTERS", ylab="LOG ( SCORE )")
  names(score) <- seq(from=3, to=N_pt, by=1)
  opt_score <- min(as.numeric(names(score)[score==min(score)]))
  print(paste0("OPTIMAL MIN SIZE OF CLUSTERS = ",opt_score))
  mtext(paste0("OPTIMAL MIN SIZE OF CLUSTERS = ",opt_score))
  
  ######################################## FINAL TSNE PLOT AND CLUSTERING #########################################
  print("PERFORMING FINAL TSNE DIMENSIONALITY REDUCTION AND HDBSCAN CLUSTERING")
  if(N_cells <= 200){N_tsne <- 120}else{N_tsne <- 60}
  tsne_out <- list(length = N_tsne)
  KL <- vector(length = N_tsne)
  for(k in 1:N_tsne)
  {
    tsne_out[[k]]<-Rtsne(t(log10(expr+1)),initial_dims=optPC,verbose=FALSE,check_duplicates=FALSE,
                         perplexity=optPerp,dims=2,max_iter=10000)
    KL[k]<-tail(tsne_out[[k]]$itercosts,1)
  }
  names(KL)<-c(1:N_tsne)
  opt_tsne<-tsne_out[[as.numeric(names(KL)[KL==min(KL)])]]$Y
  res_opt <- hdbscan(opt_tsne, minPts = opt_score)
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
  plot(opt_tsne,main=paste0("TSNE PLOT: ",file_name),col=colors[as.character(res_opt$cluster)],
       xlab="tSNE1",ylab="tSNE2")
  mtext(paste0("OPTIMAL NUMBER OF CLUSTERS = ",N_clust))
  write.table(data.frame(CELL=colnames(expr),CLUSTER=res_opt$cluster),
              file=paste0("Easy_Output/CLUSTER_ASSIGNMENT_",file_name),
              col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")
  sink()
  dev.off()
}
