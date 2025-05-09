
if (T) {
  dir.create("scripts")
  dir.create("results")
  dir.create("files")
  dir.create("figures")
  dir.create("origin_datas/GEO",recursive = T)
  dir.create("origin_datas/TCGA")
}

library(stringr)
library(tidydr)
library(openxlsx)
library(data.table)
library(reshape2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(clusterProfiler)
library(pheatmap)
library(ComplexHeatmap)
library(GSVA)
library(GSEABase)
library(fgsea)
library(corrplot)
library(colorspace)
library(survival)
library(survminer)
library(maftools)
library(vegan)
library(forcats)
library(ggpubr)
library(ggsci)
library(ggplot2)
library(rstatix)
library(ggstatsplot)
library(ggcor)
library(ggstance)
options(stringsAsFactors = F)



#GSE131182#########
GSE131182.data=read.xlsx('origin_datas/GEO/GSE131182_CircRNA_Expression_Profiling.xlsx')
grep('hsa_circ_0000831',GSE131182.data$circBaseID)
colnames(GSE131182.data)
GSE131182.exp=GSE131182.data[,c(2:14)]
GSE131182.exp=GSE131182.exp %>% drop_na(circBaseID)
rownames(GSE131182.exp)=GSE131182.exp$circBaseID
GSE131182.exp=GSE131182.exp[,-13]

GSE131182.samples=data.frame(Samples=colnames(GSE131182.exp),
                             type=rep(c('healthy','OSCC'),c(6,6)))
GSE131182.samples

#############
GSE131182.Degs=read.xlsx('results/01.DEcircRNA/GSE131182_DEcircRNA.xlsx')
head(GSE131182.Degs)
pdf('results/01.DEcircRNA/GSE131182_volcano_plot.pdf',height = 10,width = 12)
ggplot(GSE131182.Degs,aes(x=Log2FC,y=-log10(`p-value`),color=Change))+
  geom_point()+
  scale_color_manual(values=c("blue","grey","red"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
    # legend.title = element_text(family = 'Times',size=15),
    text = element_text(size=20,face = 'bold'),
    plot.margin  = unit(c(2,1,1,1), 'cm'),
    legend.position = 'right')+
  labs(title = 'GSE131182 OSCC vs ctrl')+
  ylab('-log10(pvalue)')+xlab('log2(Fold change)')+
  geom_vline(xintercept=c(-1,1),linetype='dashed',col="grey",lwd=0.5) +
  geom_hline(yintercept = -log10(0.05),linetype='dashed',col="grey",lwd=0.5)
dev.off()

###############
table(GSE131182.Degs$Change)
GSE131182.heatmap=GSE131182.exp[GSE131182.Degs$circBase.ID[GSE131182.Degs$Change!='stable'],]

pdf('results/01.DEcircRNA/GSE131182_pheatmap_plot.pdf',height = 8,width = 8)
pheatmap(GSE131182.heatmap[,GSE131182.samples$Samples],
         scale = 'row',name = 'Expression',
         color = colorRampPalette(c(rev(brewer.pal(3,"RdYlBu"))))(100),
         breaks = seq(-1, 1, length.out = 100),
         display_numbers = F, 
         annotation_col = GSE131182.samples[,'type',drop=F],
         annotation_colors = list(type=setNames(c("#BEBADA", "#FDB462"),
                                                c("healthy", "OSCC"))),
         cluster_cols = T, 
         cluster_rows = T,
         show_rownames = F, 
         show_colnames = T)
dev.off()




################
circrna.data=read.xlsx('origin_datas/1_circRNA_expression.xlsx',rowNames = T)
colnames(circrna.data)
circrna.exp=circrna.data[,c(10:15)]

circrna.samples=data.frame(Samples=colnames(circrna.exp),
                             type=rep(c('Ctrl','OSCC'),c(3,3)))
circrna.samples



circrna.Degs=read.xlsx('results/01.DEcircRNA/1_OSCCVSCtrl_circRNA_differential_expression.xlsx')
###############
table(circrna.Degs$significant)
circrna.heatmap=circrna.exp[circrna.Degs$Accession[circrna.Degs$significant=='yes'],]

library(tinyarray)
pdf('results/01.DEcircRNA/circRNA_pheatmap_plot.pdf',height = 8,width = 8)
pheatmap(circrna.heatmap[,circrna.samples$Samples],
         scale = 'row',name = 'Expression',
         color = colorRampPalette(c(rev(brewer.pal(3,"RdYlBu"))))(100),
         breaks = seq(-1, 1, length.out = 100),
         display_numbers = F, 
         annotation_col = circrna.samples[,'type',drop=F],
         annotation_colors = list(type=setNames(c("#BEBADA", "#FDB462"),
                                                c("Ctrl", "OSCC"))),
         cluster_cols = T, 
         cluster_rows = T,
         show_rownames = T,
         show_colnames = T)
dev.off()


##############
head(GSE131182.Degs)
GSE131182.anno=read.xlsx('origin_datas/GEO/GSE1131182_circRNA_anno.xlsx')
colnames(GSE131182.anno)
GSE131182.up.DEcirCRNA=GSE131182.anno[GSE131182.anno$circBaseID %in% GSE131182.Degs$circBase.ID,'GeneName']
GSE131182.up.GeneName=unique(GSE131182.up.DEcirCRNA)

library(clusterProfiler)
library(org.Hs.eg.db)
DEG.entrez_id = mapIds(x = org.Hs.eg.db,
                       keys = GSE131182.up.GeneName,
                       keytype = "SYMBOL",
                       column = "ENTREZID")
DEG.entrez_id = na.omit(DEG.entrez_id)

up.DEGs.enrichKEGG=enrichKEGG(DEG.entrez_id,
                              organism = "hsa",
                              keyType = "kegg",
                              pvalueCutoff = 0.05,
                              pAdjustMethod = "BH",
                              minGSSize = 10,
                              maxGSSize = 500)
enrichplot::dotplot(up.DEGs.enrichKEGG)
up.DEGs.enrichKEGG.res=up.DEGs.enrichKEGG@result
write.xlsx(up.DEGs.enrichKEGG.res,'results/01.DEcircRNA/GSE131182_KEGG_enrichment.xlsx',overwrite = T)

head(up.DEGs.enrichKEGG.res)
table(up.DEGs.enrichKEGG.res$p.adjust<0.05)
up.DEGs.enrichKEGG.res=up.DEGs.enrichKEGG.res[up.DEGs.enrichKEGG.res$p.adjust<0.05,]
up.DEGs.enrichKEGG.res=up.DEGs.enrichKEGG.res  %>% slice_min(n =25, order_by = p.adjust)
pdf('results/01.DEcircRNA/GSE131182_KEGG_enrichment.pdf',height = 9,width = 10)
ggplot(data=up.DEGs.enrichKEGG.res,aes(x=Count,y=reorder(Description ,Count),
                                       color = -log10(p.adjust))) +
  geom_point(aes(size=Count),show.legend = T) +
  scale_color_continuous(type = "viridis")+
  labs(x='gene number',y='pathway')+theme_bw()+
  theme(text = element_text(family = 'Times',size=20))
dev.off()

DEGs.enrichGO=enrichGO(DEG.entrez_id,
                       OrgDb = "org.Hs.eg.db", 
                       ont="ALL", pvalueCutoff=0.05)

DEGs.enrichGO.res=DEGs.enrichGO@result
write.xlsx(DEGs.enrichGO.res,'results/01.DEcircRNA/GSE131182_GO_enrichment.xlsx',overwrite = T)
DEGs.enrichGO.res=DEGs.enrichGO.res %>% group_by(ONTOLOGY) %>% slice_max(n =15, order_by = Count)
DEGs.enrichGO.res$Description <- factor(DEGs.enrichGO.res$Description, 
                                        levels = DEGs.enrichGO.res$Description[order(DEGs.enrichGO.res$ONTOLOGY,decreasing = T)])

pdf('results/01.DEcircRNA/GSE131182_GO_enrichment.pdf',height = 12,width = 18)
ggplot(DEGs.enrichGO.res, aes(x = Description, y = Count, 
                              fill = ONTOLOGY, group = Description)) +
  geom_bar(stat="identity", position="dodge", colour=aes(ONTOLOGY))+
  scale_fill_manual(values =c("#80B1D3","#FDB462","#91C133"))+
  xlab('')+ggtitle('GO enrichment')+theme_bw()+
  # scale_x_discrete(labels=function(y)str_wrap(y,width = 25))+
  theme(text = element_text(family = 'Times',size = 18,face = 'bold',color='black'),
        axis.text.x = element_text(size = 18,angle = 60,hjust = 1),
        legend.position = 'top')
dev.off()




circrna.enrichGO.res=read.xlsx('results/01.DEcircRNA/OSCCVSCtrl_GO_enrichment.xlsx')
head(circrna.enrichGO.res)
circrna.enrichGO.res=circrna.enrichGO.res %>% drop_na(GO_function)
circrna.enrichGO.res=circrna.enrichGO.res %>% group_by(GO_function) %>% slice_max(n =15, order_by = S.gene.number)
circrna.enrichGO.res$GO_Term <- factor(circrna.enrichGO.res$GO_Term, 
                                        levels = circrna.enrichGO.res$GO_Term[order(circrna.enrichGO.res$GO_function,decreasing = T)])

pdf('results/01.DEcircRNA/circrna_GO_enrichment.pdf',height = 12,width = 18)
ggplot(circrna.enrichGO.res, aes(x = GO_Term, y = S.gene.number, 
                              fill = GO_function, group = GO_Term)) +
  geom_bar(stat="identity", position="dodge", colour=aes(GO_function))+
  scale_fill_manual(values =c("#80B1D3","#FDB462","#91C133"))+
  xlab('')+ggtitle('GO enrichment')+theme_bw()+ylab('Count')+
  # scale_x_discrete(labels=function(y)str_wrap(y,width = 25))+
  theme(text = element_text(family = 'Times',size = 18,face = 'bold',color='black'),
        axis.text.x = element_text(size = 18,angle = 60,hjust = 1),
        legend.position = 'top')
dev.off()



circrna.enrichKEGG.res=read.xlsx('results/01.DEcircRNA/OSCCVSCtrl_KEGG_enrichment.xlsx')
head(circrna.enrichKEGG.res)
circrna.enrichKEGG.res=circrna.enrichKEGG.res %>% slice_min(n =25, order_by = pvalue)

pdf('results/01.DEcircRNA/circRNA_KEGG_enrichment.pdf',height = 9,width = 10)
ggplot(data=circrna.enrichKEGG.res,aes(x=S.gene.number,y=reorder(pathway_name,S.gene.number),
                                       color = -log10(pvalue))) +
  geom_point(aes(size=S.gene.number),show.legend = T) +
  scale_color_continuous(type = "viridis")+
  labs(x='gene number',y='pathway')+theme_bw()+
  theme(text = element_text(family = 'Times',size=20))
dev.off()

#hsa_circ_0000831ROC#########
dir.create('results/02.hsa_circ_0000831')
GSE131182.expr.df=cbind.data.frame(hsa_circ_0000831=as.numeric(GSE131182.exp['hsa_circ_0000831',GSE131182.samples$Samples]),
                                   type=GSE131182.samples$type)

GSE131182.expr.df

pdf('results/02.hsa_circ_0000831/GSE131182_Hubgenes_expr.pdf',height = 5,width = 5)
ggboxplot(GSE131182.expr.df, x = "type", y = "hsa_circ_0000831",color="type")+
  ggpubr::stat_compare_means(aes(group=type), label = "p.format", method = 'wilcox.test')+
  scale_color_manual(values = group.cols)+ylab('Normalized expression level
')+xlab('')+
  theme(text = element_text(family = 'Times',size=15),legend.position = 'none')+
  ggtitle('hsa_circ_0000831')
dev.off()

library(pROC)
pdf('results/02.hsa_circ_0000831/GSE131182_Hubgenes_ROC.pdf',height = 5,width = 5)
GSE131182.roc1 <- plot.roc(GSE131182.expr.df$type,
                           as.numeric(GSE131182.expr.df$hsa_circ_0000831),
                           percent=TRUE, col="2")
title(main = paste0("hsa_circ_0000831 AUC=",round(GSE131182.roc1$auc[1],2),'%'))
dev.off()

circrna.expr.df=cbind.data.frame(hsa_circ_0000831=as.numeric(log2(circrna.exp['hsa_circ_0000831',circrna.samples$Samples]+1)),
                                 type=circrna.samples$type)

circrna.expr.df

pdf('results/02.hsa_circ_0000831/circrna_Hubgenes_expr.pdf',height = 5,width = 5)
ggboxplot(circrna.expr.df, x = "type", y = "hsa_circ_0000831",color="type")+
  ggpubr::stat_compare_means(aes(group=type), label = "p.format", method = 'wilcox.test')+
  scale_color_manual(values = group.cols)+ylab('Normalized expression level
')+xlab('')+
  theme(text = element_text(family = 'Times',size=15),legend.position = 'none')+
  ggtitle('hsa_circ_0000831')
dev.off()

pdf('results/02.hsa_circ_0000831/circrna_Hubgenes_ROC.pdf',height = 5,width = 5)
circrna.roc1 <- plot.roc(circrna.expr.df$type,
                 as.numeric(circrna.expr.df$hsa_circ_0000831),
                 percent=TRUE, col="2")
title(main = paste0("hsa_circ_0000831 AUC=",round(circrna.roc1$auc[1],2),'%'))
dev.off()


