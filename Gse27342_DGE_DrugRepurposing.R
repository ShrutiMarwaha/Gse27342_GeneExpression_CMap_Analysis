library(Biobase)
source("http://bioconductor.org/biocLite.R")
biocLite("GEOquery")
library(GEOquery)
biocLite("limma")
library(limma)

#########################################################################################
# GSE27342 analysis
#########################################################################################
gset27342 <- getGEO("GSE27342", GSEMatrix =TRUE)
# save and load the stored gset file in future
#load("/Users/shruti/Dropbox/SHRUTIM/Microarray/Microarray_R_Scripts/final/gse27342/gset27342.rda")

GPLid <- levels((gset27342[[1]])$platform_id)
if (length(gset27342) > 1) 
{
  #idx <- grep("GPL5175", attr(gset27342, "names")) 
  idx <- grep(GPLid, attr(gset27342, "names"))
} else 
{ 
  idx <- 1
}
gset27342 <- gset27342[[idx]]

#########################################################################################
# data muging
#########################################################################################
## look at the pattern in which the sample types are arranged. 
# pData(gset27342)[1:2,]

# log2 transform
ex27342 <- exprs(gset27342)
qx <- as.numeric(quantile(ex27342, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex27342[which(ex27342 <= 0)] <- NaN
            exprs(gset27342) <- log2(ex27342) }

fvarLabels(gset27342) <- make.names(fvarLabels(gset27342))
SampleInformation <- pData(gset27342)$title
## remove anything that is not a digit
PatientID <- as.numeric(sub(pattern="\\D+", replacement="", x=SampleInformation,ignore.case=T))
TissueType <- sub(pattern=" from gastric cancer patient [[:digit:]]*", replacement="", x=SampleInformation,ignore.case=T)
TissueType <- make.names(TissueType)

#########################################################################################
# running limma
#########################################################################################

## convert the vector into factor groups
Block <- factor(PatientID)
Treatment <- factor(TissueType,levels=c("control","tumor.tissue")) ## order of levels is important
design27342 <- model.matrix(~0+Block+Treatment)
## fit linear model for each gene given a series of arrays
fit <- lmFit(gset27342, design27342)
## tell which levels should be compared
cont.matrix <- makeContrasts(Treatmenttumor.tissue, levels=design27342)
## given a linear model fit to microarray data, compute estimated coefficients & standard errors for given set of contratsts
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)

## to get the no.of genes in each array
nrow27342 <- nrow(ex27342)
tT <- topTable(fit2, adjust="fdr", sort.by="logFC",number=nrow27342)
tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","logFC","t","B","GB_LIST","gene_assignment"))
DEG_gse27342 <- subset(tT,adj.P.Val<0.05,)

## function to add annotation to limma output
Annotation <- function(limma_output)
{
  ## add annotation 
  GeneAssignment <- as.character(limma_output$gene_assignment)
  Annotation <- data.frame(matrix(data=NA,nrow=length(GeneAssignment),ncol=5,dimnames=list(c(),c("ProbeID","GenBankId","GeneSymbol","GeneName","GeneID"))),stringsAsFactors=F)
  for(i in seq_along(GeneAssignment))
  {
    Annotation$ProbeID[i] <- limma_output$ID[i]
    ## extract information about gene annotation 
    GeneAssignmentSplit <- strsplit(GeneAssignment[i]," [//]+ ")
    Annotation$GenBankId[i] <- GeneAssignmentSplit[[1]][1]
    Annotation$GeneSymbol[i] <- GeneAssignmentSplit[[1]][2]
    Annotation$GeneName[i] <- GeneAssignmentSplit[[1]][3]
    Annotation$GeneID[i] <- GeneAssignmentSplit[[1]][5]
  }
  GseFinal <- merge(x=Annotation,y=limma_output[, -grep("GB_LIST|gene_assignment",names(limma_output))],by.x="ProbeID",by.y="ID")
  GseFinal <- GseFinal[order(-abs(GseFinal$logFC)),]
  return(GseFinal)
}
DEG_gse27342 <- Annotation(DEG_gse27342)

#function to remove rows which are not associated with any gene 
RemoveRows <- function(MatrixName,ColumnName){
  complete <- complete.cases(MatrixName)
  MatrixName <- MatrixName[complete,]
  MatrixName <- subset(MatrixName,ColumnName!="",)
  return(MatrixName)
}

DEG_gse27342 <- RemoveRows(DEG_gse27342,"GeneID")
gse27342_up <- unique(subset(DEG_gse27342,logFC>0,GeneID))
gse27342_down <- unique(subset(DEG_gse27342,logFC<0,GeneID))

#########################################################################################
# Connectivity Map analysis
#########################################################################################
# generate iterative input lists for cmap input and find compounds which are common in among cmap output for all the input lists
#load("/Users/shruti/Dropbox/SHRUTIM/Microarray/Microarray_R_Scripts/GastricCancer/GSE27342/final/2014/Gse27342DEGFinal.rda") # Gse27342DEGFinal is same as DEG_gse27342
setwd("/Users/shruti/Dropbox/SHRUTIM/Microarray/Microarray_R_Scripts/GastricCancer/GSE27342/rough/2014/cMAP/Iteration")

Gse27342CmapInput <- DEG_gse27342[,c("GeneID","logFC","adj.P.Val")]
Gse27342CmapInput <-  Gse27342CmapInput[order(-abs(Gse27342CmapInput$logFC)),]
grep("/",Gse27342CmapInput$GeneID)
# if GeneID contains more than 1 id, remove additional ones
#GeneIds <- sapply(Gse27342CmapInput$GeneID,function(x) {
#  sub("[ ]/.+","",x)
#})
#Gse27342CmapInput$GeneID <- GeneIds
write.table(Gse27342CmapInput,file="./Gse27342CmapInput.txt",sep="\t",quote=F,col.names=T,row.names=F) 
#######################################################################################

#######################################################################################
# convert entrez ids to hgu1331a afffy probeset ids.
#######################################################################################
## use biomaRt package to convert entrez ids to hgu1331a afffy probeset ids. Cross check by converting entrezIds using http://central.biomart.org/converter/#!/ID_converter/gene_ensembl_config_2
#source("http://bioconductor.org/biocLite.R")
#biocLite("biomaRt")
library("biomaRt")
listMarts()
## select biomart database and dataset
database <- useMart("ensembl")
grep("sapiens",listDatasets(database)$description,)
listDatasets(database)[32,]
#data <- useDataset("hsapiens_gene_ensembl",mart=database)
## connect to a specified biomart database
data <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
## list the filters available
filters <- listFilters(data)
filters
test <- grep("hg.?u.?133.?a",filters$description,ignore.case=T)
filters[test,]
## attribites are values that you are interested in to retrieve
attributes <- listAttributes(data)
#attributes
grep("entrez",attributes$description,ignore.case=T)
grep("hg.?u.?133.?a",attributes$description,ignore.case=T)
#attributes[c(55,56,100,101),]
attributes[c(59,60,108,109),]

GSE27342EntrezGeneID <- as.character(unique(Gse27342CmapInput$GeneID))
GSE27342probesetIDs <- getBM(attributes=c('entrezgene','affy_hg_u133a'), filters = 'entrezgene', values = GSE27342EntrezGeneID, mart = data)
GSE27342probesetIDs <- subset(GSE27342probesetIDs,affy_hg_u133a!="")
GSE27342probesetIDs_FC <- merge(x=Gse27342CmapInput[,c("GeneID","logFC")],y=GSE27342probesetIDs,by.x="GeneID",by.y="entrezgene",sort=F,all.y=F)
GSE27342probesetIDs_FC <- GSE27342probesetIDs_FC[!duplicated(GSE27342probesetIDs_FC$affy_hg_u133a),]
GSE27342probesetIDs_FC <- GSE27342probesetIDs_FC[order(-abs(GSE27342probesetIDs_FC$logFC)),]
#######################################################################################

GSE27342cMapInputUp <- subset(GSE27342probesetIDs_FC,logFC>0,affy_hg_u133a)
GSE27342cMapInputDown <- subset(GSE27342probesetIDs_FC,logFC<0,affy_hg_u133a)
## count no.of up and Down DEGs to get idea about iterations you will like to run 
dim(GSE27342cMapInputUp)
dim(GSE27342cMapInputDown)
setwd("./Iteration")

## save probe ids for Up DEGs
for(i in seq(from=1,by=100,to=600))
{
  cMapInputUp <- GSE27342cMapInputUp[1:(i+99),1]
  write.table(cMapInputUp,file=sprintf("GSE27342CmapInputUpTop%i.txt.grp",(i+99)),sep="\t",quote=F,col.names=F,row.names=F)
}

## save probe ids for Down DEGs
for(i in seq(from=1,by=100,to=400))
{
  cMapInputDown <- GSE27342cMapInputDown[1:(i+99),1]
  write.table(cMapInputDown,file=sprintf("GSE27342CmapInputDownTop%i.txt.grp",(i+99)),sep="\t",quote=F,col.names=F,row.names=F)
}
## check the last files where there can be NAs if no.of DEgs is not exactly equal to iteration interval
##########################################################

## run iterations in cmap: 100up, 100 dn; 200up, 200 dn; 300up, 300 dn; 400up, 400 dn; 600up, 400 dn;
# provide the iteration interval on which cmap results have been named
Iteration <- c(200,400,600,800,1000)
## following list will store compound names from each iteration
CompoundList <- list()

library("xlsx")
for(i in seq_along(Iteration))  
{
  ## each cmap result file has been named as "gse27342CmapResult" followed by the number of input genes
  CmapResults <- read.xlsx2(file=sprintf("gse27342CmapResultTop%i.xls",Iteration[i]),sheetIndex=1,colClasses=c(c("numeric","character"),rep("numeric",4)))
  ## select compounds which have p.value < 0.05 and negative enrichment score
  CompoundList[[i]] <- subset(CmapResults,p<0.05 & enrichment<0,cmap.name,drop=T)
}
## find common compounds between hits from each iteration
Reduce(intersect,list(CompoundList[[1]],CompoundList[[2]],CompoundList[[3]],CompoundList[[4]],CompoundList[[5]]))
CommonCompoundsGse27342 <- Reduce(intersect,list(CompoundList[[2]],CompoundList[[3]],CompoundList[[4]],CompoundList[[5]]))
save(CommonCompoundsGse27342,file="./CommonCompoundsGse27342.rda")

