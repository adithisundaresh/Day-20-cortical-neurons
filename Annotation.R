lapply(c("dplyr","Seurat","HGNChelper","patchwork","ggplot2","harmony","viridis","tidyr","scater","openxlsx","miloR","dittoSeq"), library, character.only = T)

#loading unregressed object
day20 <- readRDS(file="CellPlex/day20")

#for sc-type web tool for annotation
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
#using our annotations
gs_list = gene_sets_prepare("CellPlex/sc_type_new_AS.edits.xlsx", "Brain")
# get cell-type by cell matrix
es.max = sctype_score(scRNAseqData = day20[["RNA"]]@scale.data, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(day20@meta.data$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(day20@meta.data[day20@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(day20@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)
# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])

day20@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  day20@meta.data$customclassif[day20@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

#final annotation plot
DimPlot(day20, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif')+
  theme(text = element_text(size = 30), axis.text = element_text(size = 30))+
  ggtitle("Annotation via scType")
ggsave("CellPlex/plots/UMAP_scType_0.8_final.pdf",width=14,height=8,bg="white")

#testing GABAergic & mesodermal cell markers
mesoderm <- c("ACTA2","CNN1","MYH3","MYH8","MYL1","MYLPF","BGN","FN1","IGFBP7","COL3A1","LUM")
DotPlot(object=day20, features=mesoderm)+
  theme(axis.text.x = element_text(size = 15, angle=45, vjust = 0.8, hjust=1))
ggsave("CellPlex/DotPlot_mesoderm.png",width=8,height=6,bg="white")

DotPlot(object=day20, features=c("GSX1","GSX2","ASCL1","DLX2","DLX1","DLX5","DLX6","GAD1","GAD2","NKX2-1","LHX6"))+
  theme(axis.text.x = element_text(size = 15, angle=45, vjust = 0.8, hjust=1))
ggsave("CellPlex/DotPlot_GABAergic.png",width=8,height=6,bg="white")

#wnt-shh signaling effects
signaling <- c("PTCH1", "GLI1", "GAS1", "AXIN2", "TNFRSF19")
DotPlot(object=day20, features=signaling)
ggsave("CellPlex/DotPlot_0.8_signaling.markers.png",width=12,height=6,bg="white")

#feature plot
FeaturePlot(day20, features = c("PAX6", "EOMES", "TBR1", "TUBB3", "MKI67", "FOXG1", "IGFBP7", "DLX1", "NEUROG1"))

#average gene expression per cell line 
qPCR <- c("PAX6","NES","TBR1","NEUROG1","NEUROG2","OTX1","VIM","ZBTB16","MKI67","EOMES","BCL11B","SOX1","FOXG1","TUBB3")
avgexpr <- AverageExpression(subset(x = day20, subset = orig.ident == "Pool1"),features=qPCR,group.by="Bio_reps")
write.csv(avgexpr,"CellPlex/avgexpr_qPCRgenes.csv")

#fetal annotation
#reference data
load("CellPlex/fetal/raw_counts_mat.rdata")
meta <- read.csv("CellPlex/fetal/cell_metadata.csv")
rownames(meta) <- meta[,1]
meta <- meta[,-1]
ref <- CreateSeuratObject(raw_counts_mat, project = "GeschwindFetal", assay = "RNA", meta.data = meta)
ref <- ref[,!is.na(ref$Number_UMI)]

#pre-processing reference
ref <- NormalizeData(ref)
ref <- FindVariableFeatures(ref, selection.method = 'vst', nfeatures = 3000)
all.genes <- rownames(ref)
ref <- ScaleData(ref, features = all.genes, vars.to.regress = c('Number_UMI','Library','Donor'))
ref <- RunPCA(ref, features = VariableFeatures(object = ref), npcs=40)

#mapping
ElbowPlot(ref)
anchors <- FindTransferAnchors(reference = ref, query = day20, dims = 1:20, reference.reduction = "pca")
predictions <- TransferData(anchorset = anchors, refdata = ref$Cluster, dims = 1:20)
day20 <- AddMetaData(day20, metadata = predictions)
saveRDS(day20,file="CellPlex/day20_ref.mapped")

#adding full-forms
long <- c("ExDp1"="Excitatory deep layer 1",
"ExDp2" ="Excitatory deep layer 2",
"ExM-U"="Maturing excitatory upper enriched",
"InCGE" = "Interneuron caudal ganglionic eminence",
"InMGE" = "Interneuron medial ganglionic eminence",
"vRG" = "Ventricular radial glia",
"IP" = "Intermediate progenitors",
"oRG" = "Outer radial glia",
"Per" = "Pericytes",
"PgG2M" = "G2/M phase cycling progenitors",
"PgS" = "S phase cycling progenitors",
"ExN" = "Newborn migrating excitatory")

predicted.id.long <- c()
for(i in day20$predicted.id){
  predicted.id.long <- c(predicted.id.long,long[i])
}

day20@meta.data[["predicted.id.long"]] <- predicted.id.long

p1 <- DimPlot(day20, reduction = 'umap', group.by='predicted.id', label=T, repel=T) + 
  scale_fill_discrete(labels = day20$predicted.id.long, breaks = day20$predicted.id) +
  ggtitle("Annotation from fetal reference")+
  theme(text = element_text(size = 30), axis.text = element_text(size = 30))
ggsave("CellPlex/UMAP_fetal.map_20dims.pdf",width=12,height=8,bg="white")

#cell type proportions
tt <- as.data.frame(table(ref$Cluster))
colnames(tt) <- c("Cell_type","Percentage")
tt$Percentage <- tt$Percentage/sum(tt$Percentage)*100

bar_data <- as.data.frame(table(day20$predicted.id))#,day20$Samples
colnames(bar_data) <- c("Cell_type","Percentage")#"Sample",
bar_data$Percentage <- as.numeric(bar_data$Percentage)
bar_data$Cell_type <- as.character(bar_data$Cell_type)
bar_data$Percentage <- bar_data$Percentage*100/(sum(bar_data$Percentage))
for(i in tt$Cell_type){
  if(!(i %in% bar_data$Cell_type)){
    bar_data <-rbind(bar_data,c(i,0))
  }
}

new <- rbind(tt,bar_data)
new$group <- rep(c("Fetal","In vitro"),times=c(nrow(tt),nrow(bar_data)))
new$Percentage <- as.numeric(new$Percentage,2)

ggplot(new, aes(fill=group, y=Percentage, x=Cell_type))+ 
  geom_bar(position="dodge", stat="identity")+
  scale_x_discrete(guide = guide_axis(angle = 90))+
  labs(title="Comparison of cell type proportions", fill="Dataset")+
  xlab("Cell type")+
  ylab("Percentage of cells per dataset")+theme_bw()+
  theme(text = element_text(size = 30), axis.text = element_text(size = 30))
ggsave("CellPlex/fetal_vs_day20_bar.pdf", width=20,height=12,bg="white")