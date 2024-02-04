library(tinyarray)
library(data.table)
cpg1 <- fread("GSE94349.csv")  #读取手动下载的数据集
cpg1_pd <- fread("GSE94349_PD.csv")
#cpg1 <- geo_download("GSE94349")
cpg2 <- geo_download("GSE68015")   #用代码下载GEO数据集
dia1 <- geo_download("GSE38642")
dia2 <- geo_download("GSE41762")

cpg1 <- as.data.frame(cpg1)      
rownames(cpg1) <- cpg1$ID_REF    
cpg1 <- cpg1[,-1]                
cpg1_pd <- cpg1_pd[,-1]    
cpg1_pd <- t(cpg1_pd)            #

cpg2_pd <- cpg2$pd               #
cpg2 <- cpg2$exp                 #

dia1_pd <- dia1$pd
dia1 <- dia1$exp

dia2_pd <- dia2$pd
dia2 <- dia2$exp

find_anno("gpl570")                      #
library(hgu133plus2.db)
ids_570 <- toTable(hgu133plus2SYMBOL)    #
cpg1 <- trans_array(cpg1,ids_570)        #
cpg2 <- trans_array(cpg2,ids_570)

find_anno("gpl6244")
library(hugene10sttranscriptcluster.db)
ids_6244 <- toTable(hugene10sttranscriptclusterSYMBOL)
dia1 <- trans_array(dia1,ids_6244)
dia2 <- trans_array(dia2,ids_6244)

gene <- intersect(rownames(cpg1),rownames(cpg2))   #
gene <- intersect(gene,rownames(dia1))
gene <- intersect(gene,rownames(dia2))

cpg1 <- cpg1[gene,]             #
cpg2 <- cpg2[gene,]
dia1 <- dia1[gene,]
dia2 <- dia2[gene,]

library(stringr)
#
cpg1_pd <- cpg1_pd[str_detect(cpg1_pd,"normal")|str_detect(cpg1_pd,"adamantinomatous craniopharyngioma"),]
cpg1_pd <- data.frame(names(cpg1_pd),cpg1_pd)   #
cpg1 <- cpg1[,rownames(cpg1_pd)]                #

cpg2_pd <- cpg2_pd[str_detect(cpg2_pd$title,"ACP")|str_detect(cpg2_pd$title,"NORMAL"),]
cpg2 <- cpg2[,rownames(cpg2_pd)]

#g
group_cpg1 <- factor(ifelse(cpg1_pd$cpg1_pd=="disease state: normal","normal","CPG"),levels = c("normal","CPG")) 
group_cpg2 <- factor(ifelse(str_detect(cpg2_pd$title,"ACP"),"CPG","normal"),levels = c("normal","CPG"))
group_dia1 <- factor(ifelse(dia1_pd$characteristics_ch1.5=="status: T2D donors","DIA","normal"),levels = c("normal","DIA"))
group_dia2 <- factor(ifelse(dia2_pd$`status:ch1`=="Diabetic donor","DIA","normal"),levels = c("normal","DIA"))

m <- cbind(cpg1,cpg2,dia1,dia2)   
boxplot(m)
m <- limma::normalizeBetweenArrays(m)    #
boxplot(m)
cpg1 <- m[,colnames(cpg1)]      #
cpg2 <- m[,colnames(cpg2)]
dia1 <- m[,colnames(dia1)]
dia2 <- m[,colnames(dia2)]

group_cpg1 <- ifelse(group_cpg1=="CPG","ACP","Normal")    
group_cpg1 <- factor(group_cpg1,levels = c("Normal","ACP"))  #
fids <- data.frame(rownames(cpg1),rownames(cpg1))        #
degs_cpg1 <- get_deg_all(                  
  cpg1,
  group_cpg1,
  fids,
  logFC_cutoff = 1,
  pvalue_cutoff = 0.05,
  adjust = T,
  entriz = F,
  show_rownames = FALSE,
  pkg = 4,
  color_volcano = c("#2874C5", "grey", "#f87669")
)
degs_cpg1$plots        #

degs_cpg2 <- get_deg_all(
  cpg2,
  group_cpg2,
  fids,
  logFC_cutoff = 1,
  pvalue_cutoff = 0.05,
  adjust = T,
  entriz = F,
  show_rownames = FALSE,
  pkg = 4,
  color_volcano = c("#2874C5", "grey", "#f87669")
)
degs_cpg2$plots

degs_dia2 <- get_deg_all(
  dia2,
  group_dia2,
  fids,
  logFC_cutoff = 0.2,
  pvalue_cutoff = 0.05,
  adjust = F,
  entriz = F,
  show_rownames = FALSE,
  pkg = 4,
  color_volcano = c("#2874C5", "grey", "#f87669")
)
degs_dia2$plots

group_dia1_2 <- as.data.frame(group_dia1)
rownames(group_dia1_2) <- colnames(dia1)
group_dia1_2$V2 <- rep(1,times=63)
group_dia1_2 <- group_dia1_2[order(group_dia1_2$group_dia1),]
dia1 <- dia1[,rownames(group_dia1_2)]
group_dia1 <- factor(group_dia1_2$group_dia1,levels = c("normal","DIA"))
degs_dia1 <- get_deg_all(
  dia1,
  group_dia1,
  fids,
  logFC_cutoff = 0.2,
  pvalue_cutoff = 0.05,
  adjust = F,
  entriz = F,
  show_rownames = FALSE,
  pkg = 4,
  color_volcano = c("#2874C5", "grey", "#f87669")
)
degs_dia1$plots

save(cpg1,cpg2,dia1,dia2,cpg1_pd,cpg2_pd,dia1_pd,dia2_pd,degs_cpg1,degs_cpg2,degs_dia1,
     degs_dia2,group_cpg1,group_cpg2,group_dia1,group_dia2,file = "step1.Rdata")
