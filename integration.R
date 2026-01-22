#######################

### READ DATA ###

########################
All<- readRDS("file_to_use.rds")
head(All)
head(All@meta.data)

#set idents 
Idents(All) <- "SampleID"

#check numbers for each sample
FS30X_1 <- subset(All, ident = "FS30X1"); dim(FS30X_1) # 27265, 1168
FS30X_2 <- subset(All, ident = "FS30X2"); dim(FS30X_2) # 27265  1247

#split the object
All.list <- SplitObject(All, split.by = "SampleID") # split by sample IDs
All.list <- All.list[c("FS30X1", "FS30X2")]


#######################

### INTEGRATE DATA ###

########################
#Note that this single command replaces NormalizeData, ScaleData, and FindVariableFeatures
#Transformed data will be available in the SCT assay, which is set as the default after running sctransform
for (i in 1:length(All.list)) { # normalize data using SCTransform method, returns 3,000 variable feature
  All.list[[i]] <- SCTransform(All.list[[i]], method = "glmGamPoi", variable.features.n = 3000, verbose = FALSE)
}
All.features <- SelectIntegrationFeatures(All.list, verbose = TRUE)
All.list <- PrepSCTIntegration(All.list, anchor.features = All.features, verbose = TRUE)
All.anchors <- FindIntegrationAnchors(All.list, normalization.method = "SCT", anchor.features = All.features, dims = 1:30)
All.integrated <- IntegrateData(All.anchors, normalization.method = "SCT", dims = 1:20)
All.integrated <- RunPCA(All.integrated, npcs = 100, verbose = TRUE)
ElbowPlot(All.integrated, ndims = 100)

#SET NUMVER OF PC'S
PCdims <- 1:11

#OTHER DR 
#Create UMAP dims:
All.integrated <- RunUMAP(All.integrated, 
                          dims = PCdims,
                          n.components = 3, # calculate 3 plot dimensions in case we want to try 3D plotting later
                          assay = "SCT") # create UMAP

#Create t-SNE dimensions:
All.integrated <- RunTSNE(All.integrated, 
                          dims = PCdims, # use our calculated number of PCs
                          dim_embed = 3, # calculate 3 plot dimensions in case we want to try 3D plotting later
                          assay = "SCT")
