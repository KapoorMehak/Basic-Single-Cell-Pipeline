#######################

### CLUSTER DATA ###

########################

All.integrated <- FindNeighbors(All.integrated, dims = PCdims, verbose = TRUE)
#1. RESOLUTION 
All.integrated <- FindClusters(All.integrated, 
                               dims = PCdims, 
                               resolution = c(.5, 1, 1.5,1.85, 2, 2.5, 3, 3.5, 4, 4.5, 5), graph.name = "integrated_snn",
                               verbose = FALSE)
DefaultAssay(All.integrated) <- "RNA"
All.integrated <- NormalizeData(All.integrated,  normalization.method = "LogNormalize", scale.factor = 10000, assay = "RNA")
All.integrated <- ScaleData(All.integrated, assay = "RNA")

#VISUALIZE THE LAST RESOLUTION
DimPlot(All.integrated,label=TRUE,raster = FALSE)



##############################

### VALIDATION FOR CLUSTERS ###

################################

#2. CLUSTREE
library(clustree)
clustree(All.integrated, prefix = "integrated_snn_res.", layout = "sugiyama")
#clustree(All.integrated, prefix = "integrated_snn_res.")

DefaultAssay(All.integrated) <- "SCT"


#CHECK FOR MARKER GENES 
clustree(All.integrated, 
         prefix = "integrated_snn_res.",
         node_colour = "CD3E", 
         node_colour_aggr = "median")
clustree(All.integrated, 
         prefix = "integrated_snn_res.",
         node_colour = "CD4", 
         node_colour_aggr = "median")
DefaultAssay(All.integrated) <- "SCT"
FeaturePlot(All.integrated,
            features = c("CD4", "CD3E"))


#CHECK HOW UMAP CHANGES
DimPlot(All.integrated, 
        label = TRUE,
        reduction = 'umap',
        group.by = "integrated_snn_res.0.5")
DimPlot(All.integrated, 
        label = TRUE,
        reduction = 'umap',
        group.by = "integrated_snn_res.1.5")
DimPlot(All.integrated, 
        label = TRUE,
        reduction = 'umap',
        group.by = "integrated_snn_res.1.85")
DimPlot(All.integrated, 
        label = TRUE,
        reduction = 'umap',
        group.by = "integrated_snn_res.2")
DimPlot(All.integrated, 
        label = TRUE,
        reduction = 'umap',
        group.by = "integrated_snn_res.3")
DimPlot(All.integrated, 
        label = TRUE,
        reduction = 'umap',
        group.by = "integrated_snn_res.4")

#VALIDATE INTEGRATION
DimPlot(All.integrated, 
        label = TRUE,
        reduction = 'umap',  
        group.by = "SampleID")


#CHOOSE A RESOLUTION

# CLUSTER TREE
Idents(All.integrated)<- "integrated_snn_res.1.85"
All.integrated <- BuildClusterTree(
  All.integrated,
  dims = 1:10,  # Replace with the appropriate dimensions for your data
  assay = "PCA"  # Replace with the correct assay name or type
)

PlotClusterTree(All.integrated,
                edge.width = 5)
data.tree <- Tool(object = All.integrated,
                  slot = "BuildClusterTree") # pull the tree
ape::plot.phylo(x = data.tree,
                direction = "downwards", # plot the tree without node labels
                edge.width = 1.5)
plot(data.tree, direction = 'downwards', edge.width = 1.5, font = 1)

#SET LEVELS ACCORDING TO YOUR TREE
levels(All.integrated) <- c('18','13','14','7','17','15','0','3','9','11','2','10','5','16','1','8','12','4','6')
All.integrated$phyloorder_1.85_re <- Idents(All.integrated)


#CHECK TO SEE EXPRESSION OF MARKERS IN CLUSTERS
features = c('CD86', 'ENSSSCG00000028461', 'CD14', 'CD163', 'CSF1R', 'NLRP3', 'TLR4', 
             'FLT3', 'HLA-DRA', 'ENSSSCG00000001455', 'FCER1A', 
             'TCF4', 'XBP1', 'CLEC12A', 'CD93', 'IRF8',
             'CD79A', 'CD19', 'PAX5', 'ENSSSCG00000028674', 'IRF4', 'PRDM1', 
             'CD3E', 'CD2', 'CD4', 
             'CD8A', 'CD5', 'CD6', 'PRF1', 'KLRK1', 'KLRB1', 'TYROBP', 
             'HBM', 'AHSP')
DotPlot(All.integrated, features = features, cols = c('yellow', 'red')) + RotatedAxis()
#Can you try to make more plots ?
# Plot a legend to map colors to expression levels
FeaturePlot(All.integrated ,features = "PAX5" )
RidgePlot(All.integrated , features = "PAX5")
VlnPlot(All.integrated , features = "PAX5" )
DoHeatmap(subset(All.integrated, downsample = 100), features = features, size = 3)


#SAVE YOUR FILE
saveRDS(All.integrated,file=file.path("Practice_int.rds"))

## STATS to check
table_samples_by_clusters <- All.integrated@meta.data %>%
  group_by(SampleID, integrated_snn_res.1.85) %>%
  summarize(count = n()) %>%
  spread(integrated_snn_res.1.85, count, fill = 0) %>%
  ungroup() %>%
  mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
  dplyr::select(c('SampleID', 'total_cell_count', everything())) %>%
  arrange(factor(SampleID, levels = levels(All.integrated@meta.data$SampleID)))

knitr::kable(table_samples_by_clusters)


table_clusters_by_samples <- All.integrated@meta.data %>%
  group_by(integrated_snn_res.1.85, SampleID) %>%
  summarize(count = n()) %>%
  spread(SampleID, count, fill = 0) %>%
  ungroup() %>%
  mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
  dplyr::select(c('integrated_snn_res.1.85', 'total_cell_count', everything())) %>%
  arrange(factor(integrated_snn_res.1.85, levels = levels(All.integrated@meta.data$integrated_snn_res.1.85)))

knitr::kable(table_clusters_by_samples)



################ WHEN YOU ARE DONE, YOU CAN TRY SOME MORE TECHNIQUES #######################
#SOME EXTRA PLOTS

##1. distribution of number of transcripts and expressed genes per cell by sample
temp_labels <- All.integrated@meta.data %>%
  group_by(SampleID) %>%
  tally()

p1 <- ggplot() +
  geom_half_violin(
    data = All.integrated@meta.data, aes(SampleID, nCount_RNA, fill = SampleID),
    side = 'l', show.legend = FALSE, trim = FALSE
  ) +
  geom_half_boxplot(
    data = All.integrated@meta.data, aes(SampleID, nCount_RNA, fill = SampleID),
    side = 'r', outlier.color = NA, width = 0.4, show.legend = FALSE
  ) +
  geom_text(
    data = temp_labels,
    aes(x = SampleID, y = -Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1),
    color = 'black', size = 2.8
  ) +
  scale_color_manual(values = custom_colors$discrete) +
  scale_fill_manual(values = custom_colors$discrete) +
  scale_y_continuous(labels = scales::comma, expand = c(0.08,0)) +
  theme_bw() +
  labs(x = '', y = 'Number of transcripts') +
  theme(
    panel.grid.major.x = element_blank(),
    axis.title.x = element_blank()
  )
