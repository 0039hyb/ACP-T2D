load("step1.Rdata")
library(clusterProfiler)
##############GSEA
library(ReactomePA)
library(tidyverse)
library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
library(biomaRt)
library(enrichplot)


DEG <- degs_cpg1$deg   #
DEG <- DEG[DEG$change!="stable",]  # 
DEG <- DEG[,c(1,8)] 
gene_map <- select(org.Hs.eg.db, keys=DEG$symbol, keytype="SYMBOL", columns=c("ENTREZID")) 
colnames(DEG)[2] <- "SYMBOL"   # 
DEG<-inner_join(gene_map,DEG,by = "SYMBOL")    # 
DEG <- DEG[order(DEG$logFC,decreasing = T),]  # 
genename <- as.character(DEG[,2])  # 
geneList = DEG[,3]     # 
names(geneList) <- genename     # 

geneList  #查看
Go_gseresult <- gseGO(geneList, 'org.Hs.eg.db', keyType = "ENTREZID", ont="all", nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=1)
KEGG_gseresult <- gseKEGG(geneList, organism = "hsa",keyType = "kegg",nPerm = 1000, minGSSize = 5, maxGSSize = 1000, pvalueCutoff=0.5)
Go_Reactomeresult <- gsePathway(geneList, nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=1)
gseaplot2(Go_gseresult,"GO:0030178", pvalue_table = TRUE)  
gseaplot2(Go_gseresult,"GO:0098793", pvalue_table = TRUE)
gseaplot2(KEGG_gseresult, "hsa05200", pvalue_table = TRUE)
gseaplot2(Go_Reactomeresult, "R-HSA-6809371", pvalue_table = TRUE)


library(tinyarray)

geneset = rio::import("cellMarker.csv")     
geneset = split(geneset$Metagene,geneset$`Cell type`)    
lapply(geneset[1:3], head)   

library(GSVA)
re <- gsva(cpg1, geneset, method="ssgsea",
           mx.diff=FALSE, verbose=FALSE)
library(stringr)
Group = group_cpg1   #将group_cpg1向量赋值给Group，即样本分组信息
library(tinyarray)
draw_boxplot(re,Group,color = c("#1d4a9d","#e5171b"),width = 0.6) #

re1 <- gsva(dia1, geneset, method="ssgsea",  
           mx.diff=FALSE, verbose=FALSE)
library(stringr)
Group = group_dia1
library(tinyarray)
draw_boxplot(re1,Group,color = c("#1d4a9d","#e5171b"),width = 0.6)



library(Hmisc)
identical(colnames(re),colnames(cpg1))  
nc = t(rbind(re,cpg1[choose_gene_1se,])) #
m = rcorr(nc)$r[1:nrow(re),(ncol(nc)-length(choose_gene_1se)+1):ncol(nc)] #

p = rcorr(nc)$P[1:nrow(re),(ncol(nc)-length(choose_gene_1se)+1):ncol(nc)] #
p[1:4,1:4]  #
library(dplyr)
tmp = matrix(case_when(p<0.01~"**",  
                       p<0.05~"*",
                       T~""),nrow = nrow(p))  
library(pheatmap)
p1 <- pheatmap(t(m),
               display_numbers =t(tmp),
               angle_col =45,
               color = colorRampPalette(c("#92b7d1", "white", "#d71e22"))(100),
               border_color = "white",
               cellwidth = 20, 
               cellheight = 20,
               width = 7, 
               height=9.1,
               treeheight_col = 0,
               treeheight_row = 0)
identical(colnames(re1),colnames(dia1))
nc = t(rbind(re1,dia1[choose_gene_1se,]))
m = rcorr(nc)$r[1:nrow(re1),(ncol(nc)-length(choose_gene_1se)+1):ncol(nc)]

p = rcorr(nc)$P[1:nrow(re1),(ncol(nc)-length(choose_gene_1se)+1):ncol(nc)]
p[1:4,1:4]
library(dplyr)
tmp = matrix(case_when(p<0.01~"**",
                       p<0.05~"*",
                       T~""),nrow = nrow(p))
library(pheatmap)
p1 <- pheatmap(t(m),
               display_numbers =t(tmp),
               angle_col =45,
               color = colorRampPalette(c("#92b7d1", "white", "#d71e22"))(100),
               border_color = "white",
               cellwidth = 20, 
               cellheight = 20,
               width = 7, 
               height=9.1,
               treeheight_col = 0,
               treeheight_row = 0)



###############################################################
library(ggpubr)
expa_choose <- cpg2[choose_gene_1se,]   
lable2 <- ifelse(group_cpg2=="normal","Normal","ACP") 
WNT5B <- as.numeric(expa_choose[1,])  
box_a_WNT5B <- as.data.frame(WNT5B)   
box_a_WNT5B$TYPE <- lable2     
ggplot(box_a_WNT5B,aes(x=TYPE,y=WNT5B,color=TYPE),size=3.5)+
  geom_boxplot()+
  theme_bw()+#改变绘图主题
  stat_compare_means(aes(label = ..p.signif..),comparisons = list(c("Normal","ACP")),method="wilcox.test")+#添加检验
  xlab("")#修改横坐标

KRT14 <- as.numeric(expa_choose[2,])
box_a_KRT14 <- as.data.frame(KRT14)
box_a_KRT14$TYPE <- lable2
ggplot(box_a_KRT14,aes(x=TYPE,y=KRT14,color=TYPE),size=3.5)+
  geom_boxplot()+
  theme_bw()+#改变绘图主题
  stat_compare_means(aes(label = ..p.signif..),comparisons = list(c("Normal","ACP")),method="wilcox.test")+#添加检验
  xlab("")#修改横坐标
FGF19 <- as.numeric(expa_choose[3,])
box_a_FGF19 <- as.data.frame(FGF19)
box_a_FGF19$TYPE <- lable2
ggplot(box_a_FGF19,aes(x=TYPE,y=FGF19,color=TYPE),size=3.5)+
  geom_boxplot()+
  theme_bw()+#改变绘图主题
  stat_compare_means(aes(label = ..p.signif..),comparisons = list(c("Normal","ACP")),method="wilcox.test")+#添加检验
  xlab("")#修改横坐标
DKK1 <- as.numeric(expa_choose[4,])
box_a_DKK1 <- as.data.frame(DKK1)
box_a_DKK1$TYPE <- lable2
ggplot(box_a_DKK1,aes(x=TYPE,y=DKK1,color=TYPE),size=3.5)+
  geom_boxplot()+
  theme_bw()+#改变绘图主题
  stat_compare_means(aes(label = ..p.signif..),comparisons = list(c("Normal","ACP")),method="wilcox.test")+#添加检验
  xlab("")#修改横坐标
MMP12 <- as.numeric(expa_choose[5,])
box_a_MMP12 <- as.data.frame(MMP12)
box_a_MMP12$TYPE <- lable2
ggplot(box_a_MMP12,aes(x=TYPE,y=MMP12,color=TYPE),size=3.5)+
  geom_boxplot()+
  theme_bw()+#改变绘图主题
  stat_compare_means(aes(label = ..p.signif..),comparisons = list(c("Normal","ACP")),method="wilcox.test")+#添加检验
  xlab("")#修改横坐标
IKBKB <- as.numeric(expa_choose[6,])
box_a_IKBKB <- as.data.frame(IKBKB)
box_a_IKBKB$TYPE <- lable2
ggplot(box_a_IKBKB,aes(x=TYPE,y=IKBKB,color=TYPE),size=3.5)+
  geom_boxplot()+
  theme_bw()+#改变绘图主题
  stat_compare_means(aes(label = ..p.signif..),comparisons = list(c("Normal","ACP")),method="wilcox.test")+#添加检验
  xlab("")#修改横坐标

PLAU <- as.numeric(expa_choose[7,])
box_a_PLAU <- as.data.frame(PLAU)
box_a_PLAU$TYPE <- lable2
ggplot(box_a_PLAU,aes(x=TYPE,y=PLAU,color=TYPE),size=3.5)+
  geom_boxplot()+
  theme_bw()+#改变绘图主题
  stat_compare_means(aes(label = ..p.signif..),comparisons = list(c("Normal","ACP")),method="wilcox.test")+#添加检验
  xlab("")#修改横坐标

library(ggpubr)

expa_choose <- dia2[choose_gene_1se,]
lable2 <- ifelse(group_dia1=="normal","Normal","DIA")
WNT5B <- as.numeric(expa_choose[1,])
box_a_WNT5B <- as.data.frame(WNT5B)
box_a_WNT5B$TYPE <- lable2
ggplot(box_a_WNT5B,aes(x=TYPE,y=WNT5B,color=TYPE),size=3.5)+
  geom_boxplot()+
  theme_bw()+#改变绘图主题
  stat_compare_means(aes(label = ..p.signif..),comparisons = list(c("Normal","DIA")),method="t.test")+#添加检验
  xlab("")#修改横坐标
KRT14 <- as.numeric(expa_choose[2,])
box_a_KRT14 <- as.data.frame(KRT14)
box_a_KRT14$TYPE <- lable2
ggplot(box_a_KRT14,aes(x=TYPE,y=KRT14,color=TYPE),size=3.5)+
  geom_boxplot()+
  theme_bw()+#改变绘图主题
  stat_compare_means(aes(label = ..p.signif..),comparisons = list(c("Normal","DIA")),method="t.test")+#添加检验
  xlab("")#修改横坐标
FGF19 <- as.numeric(expa_choose[3,])
box_a_FGF19 <- as.data.frame(FGF19)
box_a_FGF19$TYPE <- lable2
ggplot(box_a_FGF19,aes(x=TYPE,y=FGF19,color=TYPE),size=3.5)+
  geom_boxplot()+
  theme_bw()+#改变绘图主题
  stat_compare_means(aes(label = ..p.signif..),comparisons = list(c("Normal","DIA")),method="t.test")+#添加检验
  xlab("")#修改横坐标
DKK1 <- as.numeric(expa_choose[4,])
box_a_DKK1 <- as.data.frame(DKK1)
box_a_DKK1$TYPE <- lable2
ggplot(box_a_DKK1,aes(x=TYPE,y=DKK1,color=TYPE),size=3.5)+
  geom_boxplot()+
  theme_bw()+#改变绘图主题
  stat_compare_means(aes(label = ..p.signif..),comparisons = list(c("Normal","DIA")),method="t.test")+#添加检验
  xlab("")#修改横坐标
MMP12 <- as.numeric(expa_choose[5,])
box_a_MMP12 <- as.data.frame(MMP12)
box_a_MMP12$TYPE <- lable2
ggplot(box_a_MMP12,aes(x=TYPE,y=MMP12,color=TYPE),size=3.5)+
  geom_boxplot()+
  theme_bw()+#改变绘图主题
  stat_compare_means(aes(label = ..p.signif..),comparisons = list(c("Normal","DIA")),method="t.test")+#添加检验
  xlab("")#修改横坐标
IKBKB <- as.numeric(expa_choose[6,])
box_a_IKBKB <- as.data.frame(IKBKB)
box_a_IKBKB$TYPE <- lable2
ggplot(box_a_IKBKB,aes(x=TYPE,y=IKBKB,color=TYPE),size=3.5)+
  geom_boxplot()+
  theme_bw()+#改变绘图主题
  stat_compare_means(aes(label = ..p.signif..),comparisons = list(c("Normal","DIA")),method="t.test")+#添加检验
  xlab("")#修改横坐标

PLAU <- as.numeric(expa_choose[7,])
box_a_PLAU <- as.data.frame(PLAU)
box_a_PLAU$TYPE <- lable2
ggplot(box_a_PLAU,aes(x=TYPE,y=PLAU,color=TYPE),size=3.5)+
  geom_boxplot()+
  theme_bw()+#改变绘图主题
  stat_compare_means(aes(label = ..p.signif..),comparisons = list(c("Normal","DIA")),method="t.test")+#添加检验
  xlab("")#修改横坐标


M<-cor(t(dia1[choose_gene_1se,]))  
library(paletteer)
library(corrplot)
my_color = rev(paletteer_d("RColorBrewer::RdYlBu"))
my_color = colorRampPalette(my_color)(10)
png("Fig3C.png", height = 800,width = 800)
corrplot(M, type="upper",
         method="pie",
         order="hclust",
         col=my_color,
         tl.col="black",
         #tl.pos = "d",
         tl.srt=45)
dev.off()






DEG=intersect(degs_cpg1$cgs$diff$diffprobes,degs_dia1$cgs$diff$diffprobes)

gene <- bitr(DEG, fromType="SYMBOL",
             toType="ENTREZID", 
             OrgDb="org.Hs.eg.db") # ??symbolת??Ϊentrz
kk <- enrichKEGG(gene$ENTREZID, 
                 organism="hsa",
                 keyType = "kegg",
                 pvalueCutoff=0.1,
                 pAdjustMethod="BH",
                 qvalueCutoff=1)
barplot(kk)   #KEGG可视化
dotplot(kk)

ego <- enrichGO(gene=gene$ENTREZID,
                OrgDb='org.Hs.eg.db',
                ont= "ALL",
                pAdjustMethod="BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
barplot(ego, split = "ONTOLOGY", font.size = 10, 
        showCategory = 5) + 
  facet_grid(ONTOLOGY ~ ., space = "free_y",scales = "free_y") + 
  scale_y_discrete(labels = function(x) str_wrap(x, width = 45))
dotplot(ego, split = "ONTOLOGY", font.size = 10, 
        showCategory = 5) + 
  facet_grid(ONTOLOGY ~ ., space = "free_y",scales = "free_y") + 
  scale_y_discrete(labels = function(x) str_wrap(x, width = 45))

save(geneList,gene,file = "1.Rdata")
