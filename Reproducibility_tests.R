lapply(c("dplyr","Seurat","HGNChelper","patchwork","ggplot2","harmony","viridis","tidyr","scater","openxlsx","miloR","dittoSeq"), library, character.only = T)

#loading reference mapped object
readRDS(file="CellPlex/day20_ref.mapped")

#miloR for batch reproducibility (after harmony correction)
day20.sce <- as.SingleCellExperiment(day20)
day20_milo <- Milo(day20.sce)
day20_milo <- buildGraph(day20_milo, k = 40, d = 40, reduced.dim = "HARMONY")
day20_milo <- makeNhoods(day20_milo, prop = 0.1, k = 40, d=40, refined = TRUE, reduced_dims = "HARMONY")

plotNhoodSizeHist(day20_milo)
ggsave("CellPlex/plots/neighbor_dist_noharm.png",width=8,height=6,bg="white")
day20_milo <- countCells(day20_milo, meta.data = as.data.frame(colData(day20_milo)), sample="Cell_lines")
head(nhoodCounts(day20_milo))

day20_design <- data.frame(colData(day20_milo))[,c("Samples", "pool", "Cell_lines")]
day20_design$Cell_lines <- as.factor(day20_design$Cell_lines)
day20_design$pool <- as.factor(day20_design$pool)
day20_design <- distinct(day20_design)
rownames(day20_design) <- day20_design$Cell_lines
day20_design
day20_milo <- calcNhoodDistance(day20_milo, d=40, reduced.dim = "HARMONY")
da_results <- testNhoods(day20_milo, design = ~ Samples + pool, design.df = day20_design)

table(da_results$SpatialFDR<0.05)
head(da_results)
da_results %>%
  arrange(SpatialFDR) %>%
  head()
#p-value distributions
ggplot(da_results, aes(PValue)) + geom_histogram(bins=50)
ggsave("CellPlex/miloR_pval.hist.png",width=8,height=6,bg="white")
ggplot(da_results, aes(logFC, -log10(SpatialFDR))) + 
  geom_point() +
  geom_hline(yintercept = 1)
ggsave("CellPlex/miloR_volcano.png",width=8,height=6,bg="white")

day20_milo <- buildNhoodGraph(day20_milo)

## Plot single-cell UMAP
umap_pl <- plotReducedDim(day20_milo, dimred = "UMAP", colour_by="pool", text_by = "Samples", 
                          text_size = 3, point_size=0.5) +
  guides(fill="none")

## Plot neighbourhood graph
nh_graph_pl <- plotNhoodGraphDA(day20_milo, da_results, layout="UMAP",alpha=0.1) 
  
umap_pl + nh_graph_pl +
  plot_layout(guides="collect")

ggsave("CellPlex/miloR_UMAP.png",width=8,height=6,bg="white")

#repeat with custom annotation
da_results <- annotateNhoods(day20_milo, da_results, coldata_col = 'predicted.id.long')
ggplot(da_results, aes('predicted.id.long_fraction')) + geom_histogram(bins=50)
#ggsave("CellPlex/miloR_celltype_dist.png",width=8,height=6,bg="white")
da_results$predicted.id.long <- ifelse(da_results$predicted.id.long_fraction < 0.7, "Mixed", da_results$predicted.id.long)

#plotting beeswarm, repeat with fetal annotation
col <- c()
for(i in 1:nrow(da_results)){
  if(da_results[i,'PValue']>0.05){
    col[i] <- 'slategrey'
  } else if(da_results[i,'PValue']<0.05 && abs(da_results[i,'logFC'])>2){
    col[i] <- 'brown'
  } else{
    col[i] <- 'blue'
  }
}
da_results['color'] <- col
ggplot(da_results,aes(x=logFC,y=predicted.id.long,color=color))+geom_jitter(width=0.35,height=0.35)+
  ggtitle("Differential abundance of cell types between batches")+
  scale_color_manual(values=c('blue','grey','brown'), labels = c("p-value<0.05,|logFC|<2","p-value>0.05","p-value<0.05,|logFC|>2"))+
  ylab("Cell types based on fetal annotation")+ theme_bw()+
  theme(text = element_text(size = 30), axis.text = element_text(size = 20))
ggsave("CellPlex/miloR_scatter_ref.pdf",width=20,height=10,bg="white")

#clone reproducibility
data <- read.csv("CellPlex/Clone_reps.csv")
data$Clone.1 <- as.factor(data$Clone.1)
cor_qPCR=cor.test(data$qPCR_Clone1,data$qPCR_Clone2, method='pearson')
cor_sc=cor.test(data$scRNA.seq_Clone1,data$scRNA.seq_Clone2, method='pearson')

ggplot(data)+geom_point(aes(x=qPCR_Clone1,y=qPCR_Clone2,color=Clone.1,size=2))+
  annotate("text", x = min(data$qPCR_Clone1), y = max(data$qPCR_Clone2),hjust=0.2,vjust=1,size=10,
    label=paste0("r=",signif(cor_qPCR$estimate,3),"\n","  p=",signif(cor_qPCR$p.val,3)))+
  geom_smooth(method=lm,aes(x=qPCR_Clone1,y=qPCR_Clone2))+
  theme(panel.background = element_blank(),legend.title=element_blank(),axis.line = element_line(color = 'black'))+
  xlab("Gene expression in first clone (Average ΔCt)")+ylab("Gene expression in second clone (Average ΔCt)")+
  ggtitle("Comparison of clones from qPCR")+
  theme(text = element_text(size = 30), axis.text = element_text(size = 20))
ggsave("CellPlex/plots/qPCR/qPCR_clones.png",width=14,height=12,bg="white")

ggplot(data)+geom_point(aes(x=scRNA.seq_Clone1,y=scRNA.seq_Clone2,color=Clone.1,size=2))+
  geom_smooth(method=lm,aes(x=scRNA.seq_Clone1,y=scRNA.seq_Clone2))+
  annotate("text", x = min(data$scRNA.seq_Clone1), y = max(data$scRNA.seq_Clone2),hjust=0.2,vjust=1,size=10,
    label=paste0("r=",signif(cor_sc$estimate,3),"\n","  p=",signif(cor_sc$p.val,3)))+
  theme(panel.background = element_blank(),legend.title=element_blank(),axis.line = element_line(color = 'black'))+
  xlab("Log normalized gene expression in first clone")+ylab("Log normalized gene expression in second clone")+
  ggtitle("Comparison of clones from scRNA-seq")+
  theme(text = element_text(size = 30), axis.text = element_text(size = 20))
ggsave("CellPlex/plots/qPCR/sc_clones.png",width=14,height=12,bg="white")

#clone reproducibility per cell type
#filtering for genes not expressed by at least 1% of cells
selected_f <- names(which(rowSums(day20@assays$RNA[]>0)>0.001*dim(day20)[1]))
day20_sub <- subset(day20, features = selected_f)
day20_sub
all.genes <- rownames(day20_sub)
df <- data.frame()
clone11 <- list()
clone61 <- list()
for(i in unique(day20_sub$customclassif)){

  data <- AverageExpression(day20_sub[,day20_sub$customclassif==i],features=all.genes,group.by="Samples",slot='scale.data')
  data <- as.data.frame(data)

  if('RNA.11.4' %in% colnames(data) & 'RNA.47.2' %in% colnames(data)){

      cor=cor.test(data$RNA.11.4,data$RNA.47.2, method='pearson')
      df <- rbind(df,c('11.4','47.2',i,cor$estimate,cor$p.value))     
      plt <- ggplot(data, aes(x=RNA.11.4 ,y=RNA.47.2))+geom_point()+
        annotate("text", x = min(data$RNA.11.4), y = max(data$RNA.47.2),hjust=0.2,vjust=1,
          label=paste0("r=", signif(cor$estimate,3)),size=7)+
        geom_smooth(method=lm)+
        labs(x="11.4",y="47.2",title=paste0('Scaled gene expression of all genes in ',i))+
        theme(text = element_text(size = 20), axis.text = element_text(size = 25))+
        theme_bw()
      ggsave(paste0("CellPlex/plots/clone_cor/clone11_cor_",i,".png"),width=8,height=6,bg="white")
      #clone11 <- c(clone11,list(plt))

  }
   if('RNA.61.1' %in% colnames(data) & 'RNA.61.2' %in% colnames(data)){

      cor=cor.test(data$RNA.61.1,data$RNA.61.2, method='pearson')
      df <- rbind(df,c('61.1','61.2',i,cor$estimate,cor$p.value))
      plt <- ggplot(data, aes(x=RNA.61.1 ,y=RNA.61.2))+geom_point()+
        annotate("text", x = min(data$RNA.61.1), y = max(data$RNA.61.2),hjust=0.2,vjust=1,
          label=paste0("r=", signif(cor$estimate,3)),size=7)+
        geom_smooth(method=lm)+
        labs(x="61.1",y="61.2",title=paste0('Scaled gene expression of all genes in ',i))+
        theme(text = element_text(size = 20), axis.text = element_text(size = 25))+
        theme_bw()
      ggsave(paste0("CellPlex/plots/clone_cor/clone61_cor_",i,".png"),width=8,height=6,bg="white")
      #clone61 <- c(clone61,list(plt))

  }

}
colnames(df) <- c('Clone1','Clone2','Cell_type','correlation','p.val')
df$adj.p.val <- p.adjust(df$p.val,"BH")

#62.4 differences
customclassif.table <- table(day20$Samples,day20$customclassif)
customclassif.table <- customclassif.table/colSums(customclassif.table)*100

#violin plots of specific qPCR genes
day20 <- SetIdent(day20, value = day20@meta.data$Samples)
VlnPlot(subset(day20, customclassif=='Intermediate progenitor cells'), features = 'NES')+geom_boxplot(width=0.2)+
  theme(text = element_text(size = 30), axis.text = element_text(size = 30))+
  ggtitle("NES expression in IPCs")+xlab(NULL)
ggsave(paste0("CellPlex/plots/VlnPlot_NES_IPCs.pdf"),width=16,height=18,bg="white")

#heatmap for qPCR genes
pdf(file="CellPlex/heatmap_final.pdf", width=8,height=6)
dittoHeatmap(day20, genes=qPCR, annot.by = c("customclassif","Samples","pool"),show_rownames = FALSE,
cluster_rows=F, cluster_cols=F, scale='row', breaks = seq(-5, 5, length.out = 51), slot='scale.data')
dev.off()