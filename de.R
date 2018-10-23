 
args=commandArgs(T)
# set a CRAN mirror 
local({r <- getOption("repos") 
r["CRAN"] <- "http://mirrors.tuna.tsinghua.edu.cn/CRAN/"
options(repos=r)})

model=args[6]
if(is.na(args[7])){
  mycount<-read.table("input/countfile.txt",sep="\t",header=T,row.names = 1,stringsAsFactors = F)
}else{
  mycount<-read.table(args[7],sep="\t",header=T,row.names = 1,stringsAsFactors = F)
}
if(is.na(args[8])){
  mycondition<-read.table("input/condition.txt",sep="\t",header=T,row.names = 1,stringsAsFactors = F)
  #mycondition0<-read.table("condition.txt",sep="\t",header=T,stringsAsFactors = F)
}else{
  mycondition<-read.table(args[8],sep="\t",header=T,row.names = 1,stringsAsFactors = F)
  #mycondition0<-read.table(args[8],sep="\t",header=T,stringsAsFactors = F)
}
if(is.na(args[9])){
  padj=0.05
}else{
  padj=args[9]
}
options(digits=3, width=100)
type <- factor(mycondition$Treatment)
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
   
  write.csv(resdata,file = paste(uty[1],"_vs_",uty[2],"_deseq2.csv",sep=""))
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
  write.csv(tTag,file = paste(uty[1],"_vs_",uty[2],"_edger.csv",sep=""))
}else if(model=="scater"){
   
  suppressMessages(require("knitr"))
   
  suppressMessages(require("scater"))  
  suppressMessages(require("destiny")) 
  example_sce <- SingleCellExperiment(
    assays = list(counts = as.matrix(mycount)),
    colData = mycondition
  )
  example_sce <- normalize(example_sce)
  pdf("Expression_plot.pdf")
  plotExpression(example_sce, rownames(example_sce)[1:6],
                 colour_by = "Cell_Cycle", shape_by = "Mutation_Status",
                 x = "Treatment", exprs_values = "logcounts")
  dev.off()
  example_sce <- runPCA(example_sce)
  reducedDimNames(example_sce)
  pdf("PCA_plot.pdf")
  #plotPCA(example_sce)
  plotReducedDim(example_sce, use_dimred = "PCA", 
                 colour_by = "Treatment", shape_by = "Mutation_Status")
  dev.off()
  pdf("tsne_plot.pdf")
  example_sce <- runTSNE(example_sce, perplexity=10, rand_seed=1000,
                         use_dimred="PCA")
  plotTSNE(example_sce, colour_by = "Gene_0001", size_by = "Gene_1000")
  dev.off()
  example_sce <- runDiffusionMap(example_sce)
  pdf("plotDif.pdf")
  plotDiffusionMap(example_sce, colour_by = "Gene_0001", size_by = "Gene_1000")
  dev.off()
  }
  
