# introduction to seurat - set Rstudio, libpaths, install packages
# introduction to seurat and rds object, read, load data
# integration - SCtransform: Normalization, scaling, find variable features
# DR - pca, umap, tsne
# clustering - elbow plot,resolution, k.param, marker genes, silhouttee width, clustree
# visualizations - exercise : dotplots, heatmaps, umap's, feature plot, viloin plots

#######################

### CHECK LIB PATH ###

########################
.libPaths()

#######################

### LOAD LIBRARIES ###

########################
library(Seurat)
library(dplyr)
library(tidyr)
library(Matrix)
library(ggplot2)
library(patchwork)
library(DropletUtils)
library(ggpol)
library(gghalves)
library(tibble)


#############################

### DEFINE COLOR PALLETE  ###

############################
custom_colors <- list()

colors_dutch <- c(
  '#FFC312','#C4E538','#12CBC4','#FDA7DF','#ED4C67',
  '#F79F1F','#A3CB38','#1289A7','#D980FA','#B53471',
  '#EE5A24','#009432','#0652DD','#9980FA','#833471',
  '#EA2027','#006266','#1B1464','#5758BB','#6F1E51'
)

colors_spanish <- c(
  '#40407a','#706fd3','#f7f1e3','#34ace0','#33d9b2',
  '#2c2c54','#474787','#aaa69d','#227093','#218c74',
  '#ff5252','#ff793f','#d1ccc0','#ffb142','#ffda79',
  '#b33939','#cd6133','#84817a','#cc8e35','#ccae62'
)

custom_colors$discrete <- c(colors_dutch, colors_spanish)

##############################

### READ IN TO DIRECTORIES ###

#############################

data_dir <- c(
  FS30X1_1_comb = "1-1_FS30X_comb",
  FS30X2_1_comb = "2-1_FS30X_comb")
library_id<- c("FS30X1_1_comb", "FS30X2_1_comb")
scRNA_data <- Read10X(data.dir = data_dir)

#count data 
seurat_object = CreateSeuratObject(counts = scRNA_data)
count <-  as.matrix(GetAssayData(object = seurat_object, layer = 'counts'))
keep <- rowSums(count) > 0
count <- count[keep,]


# Feature data - no of gene detected/cell and mito removal
feature <- data.frame(ID = rownames(count))
rownames(feature) <- feature$ID
con <- gzfile(file.path(data_dir[1], "features.tsv.gz"))
mitoGenes <- read.table("mitogene_id.csv",header = T)
feature$Mitochondrial <- feature$EnsemblID %in% mitoGenes$x
table(feature$Mitochondrial)

#barcode data 
barcode <-data.frame(barcode = colnames(seurat_object))
barcode$SampleID <- gsub("([^_]+).+", "\\1", barcode$barcode, perl = TRUE)
barcode$BarBak <- pDat$barcode
pDat <- pDat %>% separate(BarBak, c("Sam","Loupe"))
pDat <- pDat[,-3] # remove Sam column
for (i in seq_along(library_id)){ 
  barcode$Loupe <- ifelse(barcode$SampleID == library_id[i], paste0(barcode$Loupe, paste0("-",i)), barcode$Loupe)
}
barcode$GenesDetected <- colSums(count!=0)

##plot 
ggplot(barcode, aes(x=SampleID, y=GenesDetected, fill= SampleID)) +
  geom_violin(draw_quantiles=1)+
  ylab("Total number of genes detected") +RotatedAxis()
pDat$UmiSums<- colSums(count)
ggplot(pDat, aes(x=SampleID,y=UmiSums, fill=SampleID)) +
  geom_violin(draw_quantiles=0.5)+
  ylab("Total number of molecules")
mMito <- count[feature$Mitochondrial,]
barcode$prcntMito <- colSums(mMito)/colSums(count)
barcode$DuplicatedBarcodes <- duplicated(rownames(barcode)) | duplicated(rownames(barcode), fromLast = TRUE)
table(barcode$DuplicatedBarcodes)

ggplot(barcode, aes(x=prcntMito,y=..density..)) + geom_histogram(fill="white",color="black",bins=100) +
  scale_x_continuous(breaks =seq(0, .5, .05), lim = c(0, .5)) + ylim(0,4) + facet_wrap(~SampleID) +
  geom_vline(aes(xintercept=1),color="red",lty="longdash") + 
  RotatedAxis()
ggplot(barcode, aes(x=GenesDetected,y=..density..)) + geom_histogram(fill="white",color="black",bins=500) +
  scale_x_continuous(breaks = seq(0, 2000, 250), lim = c(0, 2000)) + RotatedAxis() +
  geom_vline(aes(xintercept=200),color="red",lty="longdash") + 
  facet_wrap(~SampleID)

ggplot(barcode, aes(x=UmiSums,y=..density..)) + geom_histogram(fill="white",color="black",bins=500) +
  scale_x_continuous(breaks = seq(0, 5000, 250), lim = c(0, 5000)) + RotatedAxis() +
  geom_vline(aes(xintercept=400),color="red",lty="longdash") + 
  facet_wrap(~SampleID)

#pass
barcode <- mutate(barcode, PassViability=prcntMito < 0.05, PassGenesDet=GenesDetected > 200, PassLibSize=UmiSums > 500,
               PassBarcodeFreq=DuplicatedBarcodes==FALSE, PassAll= PassViability & PassGenesDet & PassLibSize & PassBarcodeFreq)

#check before proceeding
feature <- feature[match(rownames(count), fDat$ID), ]
rownames(barcode) <- barcode$barcode
stopifnot(identical(as.character(rownames(barcode)),colnames(count)))
stopifnot(identical(as.character(feature$ID),rownames(count)))
out <- list()
out[["counts"]] <- count
out[["barData"]] <- barcode
out[["featureData"]] <- feature

saveRDS(out,file=file.path("Practice_QC.rds")) # this saves all of our information before filtering out low quality cells
#filter process
count <- count[,barcode$PassAll]
barcode <- barcode[barcode$PassAll,]

All <- CreateSeuratObject(counts = count, meta.data = barcode) # create Seurat object of counts & barcode data
saveRDS(All,file=file.path("file_to_use.rds"))
