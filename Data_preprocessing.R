lapply(c("dplyr","Seurat","HGNChelper","patchwork","ggplot2","harmony","viridis","tidyr","scater","openxlsx","miloR","dittoSeq"), library, character.only = T)

#read in (filtered) samples
data.dirs <- c(read.table("CellPlex/samples.txt"), sep=" ")
#remove whitespace element at the end
data.dirs <- data.dirs[1:24]

#Create individual seurat objects
allobjs <- sapply(data.dirs,function(x){
  
  data <- Read10X(data.dir = x)
  CreateSeuratObject(counts = data$`Gene Expression`,
  project = gsub(".*Day_20_|/outs.*",'',x))
  
}, simplify=F)

#merge objects
day20 <- merge(allobjs[[1]],
               y=allobjs[2:24],
               add.cell.ids=sapply(data.dirs, function(x){gsub(".*_outs/|/count.*",'',x)}),
               project="day20")

cell_lines <- lapply(sub("(.*?_.*?)_.*", "\\1", rownames(day20@meta.data)),`[[`,1)
day20@meta.data[["Cell_lines"]] <- unlist(cell_lines)

day20[["percent.mt"]] <- PercentageFeatureSet(day20, pattern = "^MT-")
head(day20@meta.data, 5)

#filter plots
plot1 <- FeatureScatter(day20, feature1 = "nFeature_RNA", feature2 = "percent.mt") + scale_y_continuous(breaks=seq(0,60,5)) +
  theme(text = element_text(size = 30), axis.text = element_text(size = 30))+
  ggtitle(NULL)+xlab("Number of features")+ylab("Percentage of mitochondrial genes")
plot2 <- FeatureScatter(day20, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
    theme(text = element_text(size = 30), axis.text = element_text(size = 30))+
  ggtitle(NULL)+xlab("Number of detected (RNA) molecules")+ylab("Number of features")
plot1 + plot2
ggsave("CellPlex/feature_scatter_day20.pdf",width=20,height=12)

cells <- day20@meta.data[order(day20$nFeature_RNA),]
filter_plt <- ggplot(cells, mapping=aes(x=seq(1:nrow(cells)),y=nFeature_RNA, color=Cell_lines))+
  geom_point()+
  geom_hline(yintercept = c(2000,10000),linetype="dashed",col="red")+
  xlab("Cells sorted by nFeature_RNA")+
  ylab("Number of expressed genes per cell")+
  facet_wrap(~orig.ident)+ theme_bw() +
  theme(text = element_text(size = 30), axis.text = element_text(size = 30), legend.position = "none")
ggsave("CellPlex/Num_expressed_genes_day20.pdf",width=20,height=12)

#calculating number of cells before and after filtering (and ratio) for each cell line, per pool
counts <- as.data.frame(table(day20@meta.data[,c(1,4)]))
counts <- counts[counts$Freq!=0,]
day20_filtered <- subset(day20, subset = nFeature_RNA > 2000 & nFeature_RNA < 8000 & percent.mt < 10)
after <- as.data.frame(table(day20_filtered@meta.data[,c(1,4)]))
counts$after <- after$Freq[after$Freq!=0]
counts$percentage <- c(counts$after/counts$Freq*100)
#percentage after filtering looks okay to proceed

#filtering & pre-processing
day20 <- subset(day20, subset = nFeature_RNA > 2000 & nFeature_RNA < 8000 & percent.mt < 10)
day20 <- NormalizeData(day20)
day20 <- FindVariableFeatures(day20, selection.method = 'vst', nfeatures = 3000)
top10 <- head(VariableFeatures(day20), 10)
plot1 <- VariableFeaturePlot(day20)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
ggsave("variable_features_day20.png",plot2,width=10,height=6,bg="white")

#scaling all genes, not just variable
all.genes <- rownames(day20)
day20 <- ScaleData(day20, features = all.genes)
#better not to regress mt or cell cycle for exploratory analysis

#PCA
day20 <- RunPCA(day20, features = VariableFeatures(object = day20))

#checking number of relevant PCs
print(day20[['pca']], dims = 1:5, nfeatures = 5)
VizDimLoadings(day20, dims = 1:2, reduction = 'pca')
ggsave("CellPlex/VizDimLoadings_day20_3k.png",width=10,height=6)

DimPlot(day20, reduction = 'pca', group.by='Cell_lines', split.by='orig.ident')
ggsave("CellPlex/DimPlot_day20_3k.png",width=10,height=6)

ElbowPlot(day20)+
  xlab("Number of principal components")+
  theme(text = element_text(size = 25), axis.text = element_text(size = 25))+
  ggtitle("Elbow plot before cell cycle regression")
ggsave("CellPlex/elbowplot_day20_3k.pdf",width=10,height=6,bg="white")

day20 <- JackStraw(day20, num.replicate = 100, dims=50)
day20 <- ScoreJackStraw(day20, dims = 1:50)
JackStrawPlot(day20, dims = 1:50)
ggsave("CellPlex/JackStrawPlot_day20_50PCs_3k.png",width=5,height=6,bg="white")

#harmony for batch correction
day20 <- RunHarmony(day20, "orig.ident")

#addding additional metadata
day20@meta.data[["Cell_lines"]] <- paste0(day20@meta.data$orig.ident,"_",day20@meta.data$Cell_lines)
Clones <- sapply(day20@meta.data$Cell_lines, function(x){
  gsub(".*Pool1_|Pool2_|_.*",'',x)
}, simplify=F)
#arbitrarily assigning all 47-2 to 11-4 and 61-2 to 61-1 (clones)
Clones <- gsub("47-2.*","11-4",Clones)
Clones <- gsub("61-2.*","61-1",Clones)
day20@meta.data[["Clones"]] <- unlist(Clones)
Tech_reps <- sapply(day20@meta.data$Cell_lines, function(x){
  gsub("(.*?_.*?)_.*", "\\1", x)
}, simplify=F)
day20@meta.data[["Tech_reps"]] <- unlist(Tech_reps)

DimPlot(day20, reduction = 'harmony', group.by='Cell_lines')
ggsave("CellPlex/Dimplot_after_harmony.png",width=10,height=6,bg="white")

#clustering
day20 <- FindNeighbors(day20, dims = 1:15)
#try different resolutions
day20 <- FindClusters(day20, resolution = 0.3)
day20 <- RunUMAP(day20, reduction="harmony", dims = 1:15)
DimPlot(day20, reduction = 'umap')+
  theme(text = element_text(size = 30), axis.text = element_text(size = 30))+
  ggtitle("UMAP based on resolution=0.3")
ggsave("CellPlex/UMAP_0.3.pdf",width=14,height=8,bg="white")

day20 <- FindClusters(day20, resolution = 0.8)
day20 <- RunUMAP(day20, reduction="harmony", dims = 1:15)
DimPlot(day20, reduction = 'umap')+
  theme(text = element_text(size = 30), axis.text = element_text(size = 30))+
  ggtitle("UMAP based on resolution=0.8")
ggsave("CellPlex/UMAP_0.8.pdf",width=14,height=8,bg="white")

#some meta.data changes
day20$pool <- day20$orig.ident
day20$orig.ident <- NULL
day20$donor <- day20$Clones
day20$Clones <- NULL
Samples <- sapply(day20@meta.data$Tech_reps, function(x){
  gsub("Pool*._", "", x)
}, simplify=F)
day20@meta.data[["Samples"]] <- unlist(Samples)

#color by alternate metadata cols
DimPlot(day20, reduction = 'umap', group.by='Samples', split.by='pool')+
  ggtitle(NULL)
ggsave("CellPlex/UMAP_after_harmony_clones.pdf",width=14,height=8,bg="white")
DimPlot(day20, reduction = 'umap', group.by='pool')+ ggtitle(NULL)
ggsave("CellPlex/UMAP_after_harmony_pool.pdf",width=14,height=8,bg="white")

#compute cell cycle scores and plot UMAP colored by cell cycle
day20 <- CellCycleScoring(day20, s.features = cc.genes.updated.2019$s.genes,
                          g2m.features = cc.genes.updated.2019$g2m.genes,
                          set.ident = FALSE)
DimPlot(day20, reduction = 'umap', group.by='Phase')
ggsave("CellPlex/UMAP_cell_cycle_after_harmony.png",width=10,height=6,bg="white")

#save object before cell cycle regression
saveRDS(day20,file="CellPlex/day20")

#regress out  difference between the G2M and S phase scores
#https://satijalab.org/seurat/articles/cell_cycle_vignette.html
DimPlot(day20, reduction = 'pca', group.by='Phase')+
  theme(text = element_text(size = 30), axis.text = element_text(size = 30))+
  ggtitle("PCA plot without cell cyle regression")
ggsave("CellPlex/PCA_before_cc_reg.pdf",width=14,height=8,bg="white")

#performing regression
day20$CC.Difference <- day20$S.Score - day20$G2M.Score
day20_cc <- ScaleData(day20, vars.to.regress = "CC.Difference", features = rownames(day20))
#Plot PCs after regression
day20_cc <- RunPCA(day20_cc, features = VariableFeatures(day20_cc), nfeatures.print = 5)
DimPlot(day20_cc, reduction = 'pca', group.by='Phase')+
  theme(text = element_text(size = 30), axis.text = element_text(size = 30))+
  ggtitle("PCA plot after cell cyle regression")
ggsave("CellPlex/PCA_cc_reg.pdf",width=14,height=8,bg="white")

ElbowPlot(day20_cc)+
  xlab("Number of principal components")+
  ggtitle("Elbow plot after cell cycle regression")+
  theme(text = element_text(size = 25), axis.text = element_text(size = 25))
ggsave("CellPlex/elbowplot_day20_cc_reg.pdf",width=10,height=6,bg="white")

#try clustering with the cell cycle regressed data
day20 <- FindNeighbors(day20, dims = 1:15)
day20 <- FindClusters(day20, resolution = 0.8)
day20 <- RunUMAP(day20, reduction = "harmony", dims = 1:15)
DimPlot(day20, reduction = 'umap', group.by='seurat_clusters')+
  theme(text = element_text(size = 30), axis.text = element_text(size = 30))+
  ggtitle("PCA plot without cell cyle regression")
ggsave("CellPlex/PCA_before_cc_reg.pdf",width=14,height=8,bg="white")
#saveRDS(day20,file="CellPlex/day20_CC_regress")

#compare clusters assigned per cell ID in unregressed vs regressed object
all(rownames(day20)==rownames(day20_cc))
cluster_comparison <- table(day20$seurat_clusters, day20_cc$seurat_clusters)
#seen that the clusters change quite a bit after regression - annotate both ways and decide
day20.markers <- FindAllMarkers(day20, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(day20.markers %>% group_by(cluster) %>% slice_max(n = 10, order_by = avg_log2FC),"CellPlex/FindAllMarkers_top10_0.3.csv")
day20.cc.markers <- FindAllMarkers(day20_cc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(day20.cc.markers %>% group_by(cluster) %>% slice_max(n = 10, order_by = avg_log2FC),"CellPlex/FindAllMarkers_top10_cc_0.3.csv")