load("step1.Rdata")

library(WGCNA)


fpkm <- cpg1  #读取基因表达数据
##########################################################################
##########################################################################
#
WGCNA_matrix = t(fpkm[order(apply(fpkm,1,mad), decreasing = T)[1:10000],]) #
##??????ȫgene
#WGCNA_matrix = t(fpkm)



datExpr <- WGCNA_matrix      #
sampleNames = rownames(datExpr)       #将样本名保存到名为 sampleNames 的变量中，以备后续使用
##??????????
sampleTree = hclust(dist(datExpr), method = "average") #
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", 
     xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)#
abline(h = 180, col = "red")
clust = cutreeStatic(sampleTree, cutHeight = 180, minSize = 10) 
table(clust)

keepSamples = (clust==1)
datExpr = datExpr[keepSamples, ]
sampleTree = hclust(dist(datExpr), method = "average") # 
sizeGrWindow(12,9)  #绘图
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="",
     xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
############################################################
############################################################
powers = c(c(1:10), seq(from = 12, to=20, by=2))   #
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

#绘图
par(mfrow = c(1,2));#
cex1 = 0.85;### 
# 
# 绘图
png(file = "SoftThreshold.png",width = 800, height = 600)
par(mfrow = c(1,2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");

abline(h=0.90,col="red") #

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

dev.off()
power = sft$powerEstimate   # 将 WGCNA 分析中自动选择的最优软阈值幂次赋值给 power 变量。该变量的值代表了网络中基因表达数据的相关性和差异性的平衡点
power
net = blockwiseModules(
  datExpr,
  power = 7,
  maxBlockSize = 6000,
  TOMType = "unsigned", minModuleSize = 100,
  reassignThreshold = 0, mergeCutHeight = 0.25,
  numericLabels = TRUE, pamRespectsDendro = FALSE,
  saveTOMs =F,
  saveTOMFileBase = "AS-green-FPKM-TOM",
  verbose = 3
)
table(net$colors)#  

#绘制基因树状图和每个基因模块的颜色信息
png(file = "moduleCluster.png", width = 1200, height = 800)

mergedColors = labels2colors(net$colors)

plotDendroAndColors(dendro = net$dendrograms[[1]], 
                    colors = mergedColors[net$blockGenes[[1]]],
                    groupLabels = "Module colors",
                    dendroLabels = FALSE, 
                    hang = 0.03,
                    addGuide = TRUE, 
                    guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

MEs = moduleEigengenes(datExpr, mergedColors)$eigengenes
MEs = orderMEs(MEs)


####################################################################
####################################################################
datTraits <- data.frame(ifelse(group_cpg1=="Normal",0,1))
rownames(datTraits) <- rownames(cpg1_pd) #
colnames(datTraits) <- "ACP"


moduleTraitCor=cor(MEs, datTraits, use="p") #
write.table(file="Step04-modPhysiological.cor.xls",moduleTraitCor,sep="\t",quote=F) #
moduleTraitPvalue=corPvalueStudent(moduleTraitCor, nSamples)  
write.table(file="Step04-modPhysiological.p.xls",moduleTraitPvalue,sep="\t",quote=F)


png(file="Module_trait_relationships.png",width=800,height=900)
textMatrix=paste(signif(moduleTraitCor,2),"\n(",signif(moduleTraitPvalue,1),")",sep="")
dim(textMatrix)=dim(moduleTraitCor)  #

labeledHeatmap(Matrix=moduleTraitCor,
               xLabels=colnames(datTraits),
               yLabels=names(MEs),
               ySymbols=names(MEs),
               colorLabels=FALSE,
               colors=blueWhiteRed(50),
               textMatrix=textMatrix,
               setStdMargins=TRUE,
               cex.text=1.2,
               cex.lab=0.9,
               zlim=c(-1,1),
               main=paste("Module-trait relationships"))

dev.off()

cpg_stage = datTraits 

modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("ACP", modNames, sep="");
names(MMPvalue) = paste("p.ACP", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, cpg_stage, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));


names(geneTraitSignificance) = paste("GS.", names(cpg_stage), sep="")
names(GSPvalue) = paste("p.GS.", names(cpg_stage), sep="")

module = "blue"

column = match(module, modNames); 

moduleColors=mergedColors    

moduleGenes = moduleColors==module;



png(file="Module_membership_vs_gene_significance.png")
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for ACP",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()

blue_gene <- colnames(datExpr)[moduleGenes]
save(blue_gene,file = "blue_gene.Rdata")







inn_gene <- intersect(degs_dia1$cgs$diff$diffgenes,degs_cpg1$cgs$diff$diffgenes)
inn <- intersect(inn_gene,blue_gene)
write(inn,file = "inn.txt")

library(data.table)
g11 <- fread("g11.csv")   
g11 <- g11$`shared name`

save(cpg1,cpg2,dia1,dia2,group_cpg1,group_cpg2,group_dia1,group_dia2,
     cpg1_pd,cpg2_pd,dia1_pd,dia2_pd,g11,file="step2.Rdata")
