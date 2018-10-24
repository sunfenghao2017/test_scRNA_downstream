setwd("D:/QQPCmgr/Desktop/sample_R")
args=commandArgs(T)
# set a CRAN mirror 
local({r <- getOption("repos") 
r["CRAN"] <- "http://mirrors.tuna.tsinghua.edu.cn/CRAN/"
options(repos=r)})
 
model=args[6]
if(is.na(args[7])){
   
  mycount<-read.table("input/HSMM_expr_matrix.txt",sep="\t",header=T,row.names = 1,stringsAsFactors = F)
  
}else{
  mycount<-read.table(args[7],sep="\t",header=T,row.names = 1,stringsAsFactors = F)
}
if(is.na(args[8])){
   
  mycondition<-read.table("input/HSMM_sample_sheet.txt",sep="\t",header=T,row.names = 1,stringsAsFactors = F)
  #mycondition0<-read.table("condition.txt",sep="\t",header=T,stringsAsFactors = F)
   
}else{
  mycondition<-read.table(args[8],sep="\t",header=T,row.names = 1,stringsAsFactors = F)
  #mycondition0<-read.table(args[8],sep="\t",header=T,stringsAsFactors = F)
}
if(is.na(args[9])){
  if(model=="monocle"){
    fd <-read.table("input/HSMM_gene_annotation.txt",sep="\t",header=T,row.names = 1,stringsAsFactors = F)
  }else{
  padj=0.05
  }
}else{
  if(model=="monocle"){
    fd <-read.table(args[9],sep="\t",header=T,row.names = 1,stringsAsFactors = F) 
  }else{
  padj=args[9]
  }
}
options(digits=3, width=100)
type <- factor(mycondition[,ncol(mycondition)])
#type <- factor(mycondition[1:20,ncol(mycondition)])
if(model=="-h"){
  print()
}else if(model=="Deseq"){
  suppressMessages(require(DESeq2))
  
  database <- round(as.matrix(mycount))
  coldata <- data.frame(row.names = colnames(database), type)
  dds <- DESeqDataSetFromMatrix(countData=database, colData=coldata, design=~type)
  dds <- DESeq(dds)
  res <- results(dds)
  table(res$padj <padj)
  res <- res[order(res$padj),]
  resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)
  uty<-unique(type)
   
  write.csv(resdata,file = paste("output/",uty[1],"_vs_",uty[2],"_deseq2.csv",sep=""))
}else if(model=="edger"){
  suppressMessages(require(edgeR))
  exprSet <- mycount[rowSums(mycount > 1) >= 2,]
  exprSet <- DGEList(counts = exprSet, group = type)
  exprSet <- calcNormFactors(exprSet)
  exprSet <- estimateCommonDisp(exprSet)
  exprSet <- estimateTagwiseDisp(exprSet)
  et <- exactTest(exprSet)
  tTag <- topTags(et, n=nrow(exprSet))
  tTag <- as.data.frame(tTag)
  write.csv(tTag,file = paste("output/",uty[1],"_vs_",uty[2],"_edger.csv",sep=""))
}else if(model=="scater"){
   
  suppressMessages(require("knitr"))
   
  suppressMessages(require("scater"))  
  suppressMessages(require("destiny")) 
  example_sce <- SingleCellExperiment(
    assays = list(counts = as.matrix(mycount)),
    colData = mycondition
  )
  example_sce <- normalize(example_sce)
  pdf("output/Expression_plot.pdf")
  plotExpression(example_sce, rownames(example_sce)[1:20])
  dev.off()
  example_sce <- runPCA(example_sce)
  reducedDimNames(example_sce)
  pdf("output/PCA_plot.pdf")
  #plotPCA(example_sce)
  plotReducedDim(example_sce, use_dimred = "PCA", 
                 colour_by = "State", shape_by = "Media")
  dev.off()
  pdf("output/tsne_plot.pdf")
  example_sce <- runTSNE(example_sce, perplexity=10, rand_seed=1000,
                         use_dimred="PCA")
  plotTSNE(example_sce, colour_by = row.names(example_sce)[1], size_by = row.names(example_sce)[nrow(example_sce)])
  dev.off()
  example_sce <- runDiffusionMap(example_sce)
  pdf("output/plotDif.pdf")
  plotDiffusionMap(example_sce, colour_by = row.names(example_sce)[1], size_by = row.names(example_sce)[nrow(example_sce)])
  dev.off()
}else if(model=="Seurat"){
  suppressMessages(require("Seurat"))
  suppressMessages(require("dplyr"))
  suppressMessages(require("Matrix"))
  suppressMessages(require("reshape2"))
  #mytb<-read.table("HSMM_expr_matrix.xls",sep="\t",header=T,stringsAsFactors = F)
   
  #mysp<-read.table("HSMM_sample_sheet.xls",sep="\t",header=T,stringsAsFactors = F)
  #myan<-read.table("HSMM_gene_annotation.xls",sep="\t",header=T,stringsAsFactors = F)
  mytb<-mycount
  mysp<-mycondition
  mxx <- CreateSeuratObject(raw.data = mytb, min.cells = 3, min.genes = 200, 
                             project = "10X_PBMC")
  
  mito.genes <- grep(pattern = "^MT-", x = rownames(x = mxx@data), value = TRUE)
  percent.mito <- Matrix::colSums(mxx@raw.data[mito.genes, ]) / Matrix::colSums(mxx@raw.data)
  pdf("output/QC_vlnplot.pdf")
  mxx <- AddMetaData(object = mxx, metadata = percent.mito, col.name = "percent.mito")
  VlnPlot(object = mxx, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
  dev.off()
  pdf("output/QC_geneplot.pdf")
  par(mfrow = c(1, 2))
  GenePlot(object = mxx, gene1 = "nUMI", gene2 = "percent.mito")
  GenePlot(object = mxx, gene1 = "nUMI", gene2 = "nGene")
  dev.off()
  mxx <- FilterCells(object = mxx, subset.names = c("nGene", "percent.mito"), low.thresholds = c(200, -Inf), high.thresholds = c(2500, 0.05))
  mxx <- NormalizeData(mxx)
  pdf("output/dispersion.pdf")
  mxx <- FindVariableGenes(mxx,mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
  dev.off()
  mxx <- ScaleData(mxx, vars.to.regress = c("nUMI", "percent.mito"))
  disnm<-ncol(mxx@data)-1
  if(disnm>15){
    disnm=15
  }
  mxx <- RunPCA(mxx,pc.genes = mxx@var.genes, do.print = F,pcs.compute=disnm)
  length(x = mxx@var.genes)
  pdf("output/vizpca.pdf")
  VizPCA(object = mxx, pcs.use = 1:2)
  dev.off()
  pdf("output/dot_PCA.pdf")
  PCAPlot(object = mxx, dim.1 = 1, dim.2 = 2)
  dev.off()
  pdf("output/PCheatmap.pdf")
  PCHeatmap(object = mxx, pc.use = 1:15, cells.use = disnm, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
  dev.off()
  mxx <- JackStraw(object = mxx, num.replicate = 100) 
  pdf("output/Sig_PCplot.pdf")
  JackStrawPlot(object = mxx, PCs = 1:disnm)
  dev.off()
  pdf("output/sandstone.pdf")
  PCElbowPlot(object = mxx)
  dev.off()
  mxx <- FindClusters(object = mxx, reduction.type = "pca", dims.use = 1:15, resolution = 0.6, print.output = 0, save.SNN = TRUE)
  mxx <- RunTSNE(object = mxx, reduction.use = "pca",dims.use = 1:15, do.fast = TRUE)
  pdf("output/tsneplot_Seurat.pdf")
  TSNEPlot(object = mxx)
  dev.off()
  save(mxx, file = "test.rData")
  
  mxx.markers <- FindAllMarkers(object = mxx, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
  mxx.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
  mxx.markers %>% group_by(cluster) %>% top_n(10, avg_logFC) -> top10
 pdf("output/markergenes.pdf")
  DoHeatmap(object =mxx, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE) 
  dev.off()
  pdf("output/featureplot.pdf")
  FeaturePlot(object = mxx, features.plot = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"), cols.use = c("grey", "blue"), reduction.use = "tsne")
  dev.off()
  pdf("output/vlnfeature.pdf")
  VlnPlot(object = mxx, features.plot = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"), use.raw = TRUE, y.log = TRUE)
  dev.off()
  
}else if(model=="monocle"){
  suppressMessages(require("monocle"))
  suppressMessages(require("scater", quietly = TRUE))
  suppressMessages(require("knitr"))
  options(stringsAsFactors = FALSE)
  pd <- new("AnnotatedDataFrame", data = mycondition)
  fd <- new("AnnotatedDataFrame", data = fd)
   
  HSMM<-newCellDataSet(as.matrix(mytb),phenoData =pd,featureData=fd,lowerDetectionLimit = 0.1,expressionFamily = negbinomial.size())
  HSMM
  HSMM <- estimateSizeFactors(HSMM)
  HSMM <- estimateDispersions(HSMM,method="blind",sharingMode="fit-only")
  HSMM <- detectGenes(HSMM, min_expr = 0.1)
  print(head(fData(HSMM)))
  fData(HSMM)
  expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= 10))
  length(expressed_genes)
  print(head(pData(HSMM))) 
  pData(HSMM)$Total_mRNAs <- Matrix::colSums(HSMM@assayData$exprs)
   
  HSMM <- HSMM[,pData(HSMM)$Total_mRNAs < 1e6]
  upper_bound <- 10^(mean(log10(pData(HSMM)$Total_mRNAs)) +
                       2*sd(log10(pData(HSMM)$Total_mRNAs)))
  lower_bound <- 10^(mean(log10(pData(HSMM)$Total_mRNAs)) -
                       2*sd(log10(pData(HSMM)$Total_mRNAs)))
  table(pData(HSMM)$Hours)
  pdf("output/total_mRNA.pdf")
  Hour<-factor(pData(HSMM)$Hours)
  qplot(Total_mRNAs, data = pData(HSMM), color =Hour, geom = "density") +
    geom_vline(xintercept = lower_bound) +
    geom_vline(xintercept = upper_bound)
 dev.off()
  HSMM <- HSMM[,pData(HSMM)$Total_mRNAs > lower_bound & 
                 pData(HSMM)$Total_mRNAs < upper_bound]   
  L <- log(HSMM@assayData$exprs[expressed_genes,])
  melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))
  pdf("output/melt.pdf")
  qplot(value, geom="density", data=melted_dens_df) +  stat_function(fun = dnorm, size=0.5, color='red') + 
    xlab("Standardized log(FPKM)") +
    ylab("Density")
  plot(melted_dens_df)
  dev.off()
  estimateDispersions(HSMM, modelFormulaStr = "~ 1",
                      relative_expr = TRUE, min_cells_detected = 1, remove_outliers = TRUE,
                      cores = 1)
  disp_table <- dispersionTable(HSMM)
  unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
  HSMM <- setOrderingFilter(HSMM, unsup_clustering_genes$gene_id)
  pdf("output/order.pdf")
  plot_ordering_genes(HSMM)
  dev.off()
  pdf("output/pc_var.pdf")
  plot_pc_variance_explained(HSMM, return_all = F)
  dev.off()
  HSMM <- reduceDimension(HSMM, max_components = 2, num_dim = 2,
                          reduction_method = 'tSNE', verbose = T)
  HSMM <- clusterCells(HSMM, num_clusters=2)
  pdf("output/tsne_monocle.pdf")
  plot_cell_clusters(HSMM, 1, 2,color = "CellType", markers = c("MYF5", "ANPEP"))
  dev.off()
  pdf("output/tsne_celtyp_monocle.pdf")
  plot_cell_clusters(HSMM, 1, 2, color = "Cluster") + facet_wrap(~CellType)
  dev.off()
  
  diff_test_res <- differentialGeneTest(HSMM)
  disp_table <- dispersionTable(HSMM)
  HSMM <- reduceDimension(HSMM, max_components=2, num_dim = 6, 
                          reduction_method = 'tSNE', verbose = T) 
  cth <- newCellTypeHierarchy()
   
  summary(fData(HSMM)$num_cells_expressed)
  
   
  ordering_genes <- subset(disp_table, 
                           mean_expression >= 0.5 & 
                             dispersion_empirical >= 1 * dispersion_fit)$gene_id
  HSMM <- setOrderingFilter(HSMM, ordering_genes)
  plot_ordering_genes(HSMM)
  HSMM <- reduceDimension(HSMM, max_components=2)
  HSMM <- orderCells(HSMM)
  pdf("output/trajectory_state.pdf")
  plot_cell_trajectory(HSMM, color_by="State")
  dev.off()
  pdf("output/trajectory_cluster.pdf")
  plot_cell_trajectory(HSMM, color_by="Cluster")
  dev.off()
  
  }
  
