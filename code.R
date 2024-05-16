###Load the required libraries #################################################
# 我的R版本
# R version 4.2.1 (2022-06-23 ucrt)
library(CellChat)
library(patchwork)
library(data.table)
library(tidyverse)
library(ggpubr)
library(Seurat)

setwd("D:\\Data\\TRM\\3")
###Create a directory to save figures ##########################################
# data.dir <- 'D:\\Data\\TRM\5\\comparison"'
# dir.create(data.dir)
# setwd(data.dir)



###Load CellChat object of each dataset and then merge together ################
# Please do updateCellChat if you have CellChat objects 
# that are obtained using the earlier version (< 1.6.0).
cellchatT<- readRDS("chatTumor.rds")
cellchatN <- readRDS("chatNormal.rds")
AAA=c(levels(cellchatT@idents))
levels(cellchatN@idents)=factor(cellchatN@idents,ordered =T,levels =AAA )
cellchat@idents[["Normal"]]=factor(cellchat@idents[["Normal"]],ordered =T,levels = c(levels(cellchatT@idents)))
object.list <- list(Normal = cellchatN, Tumor = cellchatT)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

cellchat

###Part I: Predict general principles of cell-cell communication ###############
#1. Compare the total number of interactions and interaction strength
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

levels(cellchat@idents) #展示以下现在的细胞分组
#2. Compare the number of interactions and interaction strength among different cell populations
#red  (or blue) colored edges represent increased (or decreased) signaling in the second dataset compared to the first one
#pdf("上皮和TRM互作.pdf",width = 15,height = 15)

#par(mfrow = c(2,2), xpd=TRUE)
P1=netVisual_diffInteraction(cellchat, weight.scale = T,sources.use = c(3),
                             targets.use = c(2,5))
P2=netVisual_diffInteraction(cellchat, weight.scale = T,sources.use = c(2,5),
                             targets.use = c(3))
P3=netVisual_diffInteraction(cellchat, weight.scale = T,sources.use = c(3),
                             targets.use = c(2,5), measure = "weight")
P4=netVisual_diffInteraction(cellchat, weight.scale = T,sources.use = c(2,5),
                             targets.use = c(3), measure = "weight")

netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

# in a greater details using a heatmap
# The top colored bar plot represents the sum of column of values displayed in the heatmap (incoming signaling). 
# The right colored bar plot represents the sum of row of values (outgoing signaling)
gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
gg1 + gg2


# 下面的图：不是在图中展示谁比谁多，而是在图中展示每个数据集的实际互作数量
# To better control the node size and edge weights of the inferred networks across different datasets, 
# we compute the maximum number of cells per cell group and the maximum number of interactions (or interaction weights) across all datasets
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
# object.list[[1]]@idents %>% table()
# object.list[[2]]@idents %>% table()
# object.list[[1]]@net$count %>% max()
# object.list[[2]]@net$count %>% max()

par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}



###Part II: Identify the conserved and context-specific signaling pathways
#1. Identify signaling networks with larger (or less) difference as well as signaling groups based on their functional/structure similarity
# Functional similarity: High degree of functional similarity indicates major senders and receivers are similar, 
# and it can be interpreted as the two signaling pathways or two ligand-receptor pairs exhibit similar and/or redundant roles. 
# NB: Functional similarity analysis is not applicable to multiple datsets with different cell type composition.

#1.1 Identify signaling groups based on their functional similarity
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)

#1.2 Identify signaling groups based on structure similarity
#1.3 Compute and visualize the pathway distance in the learned joint manifold
# identify the signaling networks with larger (or less) difference based on their Euclidean distance in the shared two-dimensions space.
# Those signaling pathways that are only identified in one dataset are not considered here.
rankSimilarity(cellchat, type = "functional")


#2. Identify and visualize the conserved and context-specific signaling pathways
#2.1 Compare the overall information flow of each signaling pathway
# We can identify the conserved and context-specific signaling pathways by simply comparing the information flow for each signaling pathway, 
# which is defined by the sum of communication probability among all pairs of cell groups in the inferred network (i.e., the total weights in the network).
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2
dev.off()
#2.2 Compare outgoing (or incoming) signaling associated with each cell population
library(ComplexHeatmap)
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
pdf("outgoing.pdf",width = 10,height = 20)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 20)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 20)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
#dev.off()
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 26, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 26, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))



###Part III: Identify the upgulated and down-regulated signaling ligand-receptor pairs
# Identify dysfunctional signaling by comparing the communication probabities
gg1 <-netVisual_bubble(cellchat, sources.use = c(2,3), targets.use = c(1),  comparison = c(1, 2), angle.x = 45)
gg2 <-netVisual_bubble(cellchat, sources.use = c(1), targets.use = c(2,3),  comparison = c(1, 2), angle.x = 45)
gg1 + gg2
ggsave("细胞通讯气泡图.pdf", plot = gg1 + gg2, width = 10, height = 6.5) 

gg1 <- netVisual_bubble(cellchat, sources.use = c(1), targets.use = c(2,3),  comparison = c(1, 2),  title.name = "Increased signaling in Day5", angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = c(2,3), targets.use = c(1),  comparison = c(1, 2),  title.name = "Decreased signaling in Day5", angle.x = 45, remove.isolate = T)
gg1 + gg2


tmpp2=netVisual_bubble(cellchat, sources.use = c(1), targets.use = c(2,3), comparison = c(1, 2),remove.isolate = FALSE)
tmpp2=tmpp2$data
###Part IV: Visually compare cell-cell communication using Hierarchy plot, Circle plot or Chord diagram
group.cellType <- A[c(2,3,5)] # grouping cell clusters into fibroblast, DC and TC cells
names(group.cellType) <-levels(object.list[[1]]@idents)
#########################
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
################################
df.net <- subsetCommunication(cellchat)
write.csv(df.net$Tumor, "cell-cell_communications.Tumor.csv")
write.csv(df.net$Normal, "cell-cell_communications.Normal.csv")

pathways.show <- c("MHC-I")

pdf(paste0(pathways.show,".pdf"),width = 16,height = 8)

par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  # i=1
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord",
                      #edge.weight.max = weight.max[1], edge.width.max = 10,
                      group = "CXCL13+CD8+TRM",
                      signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()
####################################
pathways.show <- c("TIGIT")
par(mfrow = c(1, 2), xpd=TRUE)
# compare all the interactions sending from Inflam.FIB to DC cells
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], sources.use = c(2,3), targets.use = c(1), lab.cex = 0.5, title.name = paste0("", names(object.list)[i]))
}
pathways.show <- c("TIGIT") 
par(mfrow=c(1,2))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
#和弦图 
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")

dev.off()
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
#ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))

"hierarchy", "circle", "chord"

par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[[2]],
                     sources.use = 4,targets.use = c(5:11),
                     slot.name = 'net',net = net.up, 
                     lab.cex = 0.8,small.gap = 3.5,
                     title.name = paste0("Up-regulated signaling in ",names(object.list)[2]))
netVisual_chord_gene(object.list[[1]],
                     sources.use = 4, targets.use = c(5:11),
                     slot.name = 'net', net = net.down,
                     lab.cex = 0.8, small.gap = 3.5,
                     title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
###Save the merged CellChat object #############################################
library(RColorBrewer)
#library(cellchat)

# 设置自定义配色方案
my_colors <- brewer.pal(n = 8, name = "Set1")
my_colors=c("#E78AC3","#8DA0CB")
cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("Normal", "Tumor")) # set factor level
plotGeneExpression(cellchat, signaling = "TIGIT", split.by = "datasets",# type ="bar",
                   color.use = my_colors,colors.ggplot = T)



saveRDS(cellchat, file = "cellchat_comparisonAnalysis.rds")



### 参考 #######################################################################
# https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/Comparison_analysis_of_multiple_datasets.html
###Part IV: Visually compare cell-cell communication using Hierarchy plot, Circle plot or Chord diagram
pathways.show <- c("CXCL")
weight.max <- getMaxWeight(cellchatT, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, 
                      layout = "hierarchy",  
                      signaling.name = paste(pathways.show, names(object.list)[i]))
}
par(mfrow = c(1, 2), xpd=TRUE)
# compare all the interactions sending from Inflam.FIB to DC cells
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], sources.use = c(2,3), targets.use = c(1), lab.cex = 0.5,  title.name = paste0("Signaling from Inflam.FIB - ", names(object.list)[i]))
}



group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4)) # grouping cell clusters into fibroblast, DC and TC cells
names(group.cellType) <- levels(object.list[[1]]@idents)



pathways.show <- c("CXCL") 
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))





### R version 4.1.0 (2021-05-18) ###############################################
library(tidyverse)
library(Seurat)
library(CellChat)
library(RColorBrewer)
library(scales)

setwd("./Day5_CC_population/")



#导入cellchat对象###############################################################
#If USERS load previously calculated CellChat object (version < 1.6.0), 
#please update the object via updateCellChat
cellchat=readRDS("./../Day5.rds")
cellchat=updateCellChat(cellchat)
levels(cellchat@idents) # show factor levels of the cell labels
table(cellchat@idents) # number of cells in each cell group

#简单看看大概有啥信号 
#cellchat@netP$pathways

#设置参考数据库#################################################################
showDatabaseCategory(CellChatDB.mouse)
cellchat@DB <- CellChatDB.mouse

cellchat=cellchatT

#表达数据预处理#################################################################
##有点填补的意思（借助相关性比较高的基因，基于PPI网络），降低dropout的影响
##对于多亚基的情况比较有用（不然，如果一个受体亚基没有测到，整个受体的表达量都为0）
cellchat <- subsetData(cellchat) ##提取信号基因的表达数据，加快运行速度
##dim(cellchat@data)
##dim(cellchat@data.signaling)

#future::plan("multiprocess", workers = 6) #并行
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)



#计算通讯可能性，推断细胞互作网络###############################################
#cellchat <- computeCommunProb(cellchat) 
#默认population.size = FALSE; raw.use = TRUE
cellchat <- computeCommunProb(cellchat,population.size = T,raw.use = FALSE)
##cellchat <- filterCommunication(cellchat, min.cells = 10) 
cellchat <- computeCommunProbPathway(cellchat)
#整合细胞通讯网络###############################################################
##两种形式：counting the number of links or summarizing the communication probability
cellchat <- aggregateNet(cellchat)
##网络图########################################################################
##颜色
color.cluster=c(brewer.pal(9, "Set1")[-c(1,6)],"#8dd2c6","#fcb461","#bb80bc","#d9d9d9","#fbcde4")
names(color.cluster)=c("a_fib_0", "a_fib_1", "a_fib_2", "a_fib_3", "a_fib_4", "a_fib_5", "a_fib_6", "Bcell", "granulocytes", "mye", "NKcell", "Tcell")
##要保证不同样本、不同亚群有可比性
groupSize <- as.numeric(table(cellchat@idents))
cat(paste("max(cellchat@net$count): ",max(cellchat@net$count),"\n",sep = ""),file="for.max.txt",append=T)
cat(paste("max(cellchat@net$weight): ",max(cellchat@net$weight),"\n",sep = ""),file="for.max.txt",append=T)
cat(paste("max(groupSize): ",max(groupSize),"\n",sep = ""),file="for.max.txt",append=T)
##到这，保存运行环境后退出，比较不同样本的max，挑一个最大的
pdf('Day5.interaction.count.pdf')
netVisual_circle(cellchat@net$count, color.use = color.cluster, vertex.weight = groupSize, vertex.weight.max=8055,weight.scale = T, edge.weight.max = 81,label.edge= F, title.name = "Number of interactions")
dev.off()
pdf('Day5.interaction.strength.pdf')
netVisual_circle(cellchat@net$weight, color.use = color.cluster, vertex.weight = groupSize, vertex.weight.max=8055,weight.scale = T, edge.weight.max = 0.9,label.edge= F, title.name = "Interaction weights/strength")
dev.off()
#8055是所有样本中最大的"max(groupSize)"
#81是所有样本中最大的"max(cellchat@net$count)"
#0.9是所有样本中最大的"max(cellchat@net$weight)"



##展开图 #######################################################################
pdf('Day5.interaction.count.2.pdf',width = 12,height = 12)
mat <- cellchat@net$count
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, color.use = color.cluster, vertex.weight = groupSize, vertex.weight.max=8055,weight.scale = T, edge.weight.max = 81, title.name = rownames(mat)[i])
}
dev.off()

pdf('Day5.interaction.strength.2.pdf',width = 12,height = 12)
mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, color.use = color.cluster, vertex.weight = groupSize, vertex.weight.max=8055,weight.scale = T, edge.weight.max = 0.9, title.name = rownames(mat)[i])
}
dev.off()


setwd("E:\\data2\\TRM\\5\\comparison2")
#信号通路水平的可视化###########################################################
cellchat@netP$pathways #显示显著的通路
for (psi in c("CXCL")) {
  vertex.receiver = c(2,4) # a numeric vector
  pdf(paste('',psi,".circle.pdf",sep = ""))
  netVisual_aggregate(cellchat, #color.use = color.cluster, 
                      signaling = psi,
                      vertex.receiver = vertex.receiver,
                      layout = "circle")
  #layout: "hierarchy", "circle", "chord"
  dev.off()
}



#对于每一个信号通路，每一个受配体对的贡献#######################################
#单个受配体对的可视化
cellchat=cellchatN
for (psi in c("MHC-I")) {
  netAnalysis_contribution(cellchat, signaling = psi)
  tmpp=netAnalysis_contribution(cellchat, signaling = psi)
  tmpp=tmpp$data
  tmpp=tmpp%>%filter(contribution > 0)
  tmpp=tmpp%>%arrange(contribution)
  tmpp$name=factor(tmpp$name,levels = as.character(tmpp$name))
  tmpp%>%ggplot(aes(x=name,y=contribution))+geom_bar(stat = "identity",aes(fill=contribution))+
    scale_x_discrete("",expand = c(0,0.5))+
    scale_y_continuous("Relative contribution",expand = c(0.01,0))+
    scale_fill_gradient(low = "#c7e9c0",high = "#238b45")+
    theme_bw()+
    coord_flip()
  ggsave(paste('Day5.',psi,".contribution.pdf",sep = ""),width = 5,height = 0.5*dim(tmpp)[1])
  
  ##把最突出最显著的pair画出来
  pairLR <- extractEnrichedLR(cellchat, signaling = psi, geneLR.return = FALSE)
  LR.show <- intersect(rev(rownames(tmpp)),pairLR[,1])[1]
  
  vertex.receiver = seq(1:3)
  pdf(paste('Day5.',psi,".",LR.show,".hierarchy.pdf",sep = ""))
  netVisual_individual(cellchat,  signaling = psi,  pairLR.use = LR.show, vertex.receiver = vertex.receiver,layout = "circle")
  dev.off()
}



#多种受配体对的可视化，气泡图###################################################
netVisual_bubble(cellchat, sources.use = c(1,6), targets.use = c(10:12), remove.isolate = FALSE)
ggsave("Day5_bubble.pdf",width = 9,height = 30,units = "cm")

###
tmpp=netVisual_bubble(cellchat, sources.use = c(1,6), targets.use = c(10:12), remove.isolate = FALSE)
tmpp=tmpp$data
tmpp=tmpp[,c("prob", "pval", "interaction_name_2", "source.target")]
my.pair=read.table("./../pair.txt",header = F,sep = "\t",stringsAsFactors = F)
tmpp=tmpp[tmpp$interaction_name_2 %in% my.pair$V1,]
my.pair=setdiff(my.pair$V1,unique(tmpp$interaction_name_2))
v1=c()
v2=c()
v3=c()
v4=c()
for (sti in unique(tmpp$source.target)) {
  for (inj in my.pair) {
    v1=append(v1,NA)
    v2=append(v2,NA)
    v3=append(v3,inj)
    v4=append(v4,sti)
  }
}
my.pair=data.frame(
  prob=v1,pval=v2,interaction_name_2=v3,source.target=v4
)
tmpp=rbind(tmpp,my.pair)
tmpp$interaction_name_2=factor(
  as.character(tmpp$interaction_name_2),
  levels = rev(sort(unique(as.character(tmpp$interaction_name_2)))))

source("./../plothtm.R")
plothtm(tmpp)
ggsave(paste("Day5_bubble.2.pdf",sep = ""),width = 9,height = 15,units = "cm")


#在特定的信号通路中，每个细胞亚群扮演的角色#####################################
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
for (psi in c("CXCL")) {
  pdf(paste('Day5.',psi,".Centrality.pdf",sep = ""))
  netAnalysis_signalingRole_network(cellchat, signaling = psi, width = 8, height = 2.5, font.size = 10,color.use = color.cluster)
  dev.off()
}



#保存###########################################################################
saveRDS(cellchat, file = "Day5.res.rds")
pairLR.CXCL <- extractEnrichedLR(cellchatN, signaling = "MHC-I", geneLR.return = FALSE)
LR.show <- pairLR.CXCL[1,] # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
#> [[1]]
# Circle plot
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
###########数据读取#############
library(Seurat)
library(data.table)
library(tidyverse)
library(ggpubr)
library(fs)
files=list.files(pattern = 'txt')
filelent <- length(files)

# 设置你的根目录路径
root_dir <- getwd()

# 获取所有子文件夹路径
subdirs <- dir_ls(root_dir, type = "directory")

# 初始化一个空的列表来存储Seurat对象
scList <- list()

# 循环遍历每个子文件夹
for (subdir in subdirs) {
  # 构建filtered_feature_bc_matrix文件夹的路径
  matrix_dir <- file.path(subdir, "filtered_feature_bc_matrix")
  
  # 检查filtered_feature_bc_matrix文件夹是否存在
  if (dir_exists(matrix_dir)) {
    # 读取数据
    data <- Read10X(data.dir = matrix_dir)
    
    # 提取样本名（子文件夹名称）
    sample_name <- basename(subdir)
    
    # 创建Seurat对象，使用子文件夹名称作为project参数
    scList <- CreateSeuratObject(counts = data, project = sample_name)
    
    # 将Seurat对象存储到列表中
    seurat_objects[[sample_name]] <- scList
  } else {
    message("目录 ", matrix_dir, " 不存在。")
  }
}

#############去除双细胞#############
library(Seurat)
library(dplyr)
library(DoubletFinder) 

Find_doublet <- function(data){
  data <- NormalizeData(object = data) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData() %>%
    RunPCA()
  
  ElbowPlot(data, ndims=50)
  dim.usage <- 20
  data <- FindNeighbors(data, dims = 1:dim.usage) %>%
    FindClusters(resolution = 1) %>%
    RunUMAP(dims = 1:dim.usage) %>% 
    RunTSNE(dims = 1:dim.usage)
  
  sweep.res.list <- paramSweep_v3(data, PCs = 1:dim.usage, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  
  mpK <- as.numeric(as.vector(bcmvn[bcmvn$BCmetric==max(bcmvn$BCmetric),]$pK))
  DoubletRate = ncol(data)*8*1e-6
  nExp_poi <- round(DoubletRate*ncol(data))
  
  homotypic.prop <- modelHomotypic(data@meta.data$seurat_clusters)
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  data <- doubletFinder_v3(data, PCs = 1:dim.usage, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  colnames(data@meta.data)[ncol(data@meta.data)] = "doublet_low"
  data <- doubletFinder_v3(data, PCs = 1:dim.usage, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
  colnames(data@meta.data)[ncol(data@meta.data)] = "doublet_high"
  
  return(data)
}

# 使用for循环遍历scList
for(i in 1:length(scList)) {
  scList[[i]] <- Find_doublet(scList[[i]])
  scList[[i]] <- subset(scList[[i]], doublet_high=="Singlet")
}
# 确保 scList 是一个已命名的列表
if (is.null(names(scList))) {
  stop("scList should be a named list")
}

# 获取 scList 中的名称
add.cell.ids <- names(scList)

scRNA <- merge(scList[[1]],
               y = scList[-1], 
               add.cell.ids = add.cell.ids,
               project = "CESC")

#################质控#####################
library(tidyverse)   
library(patchwork)
dir.create('QC')
### 计算质控指标
# 计算细胞中线粒体基因比例
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
# 计算细胞中核糖体基因比例
scRNA[["percent.rb"]] <- PercentageFeatureSet(scRNA, pattern = "^RP[LS]")
# 计算红细胞比例
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1",
              "HBG1","HBG2","HBM","HBQ1","HBZ")
HB.genes <- CaseMatch(HB.genes, rownames(scRNA))
scRNA[["percent.HB"]]<-PercentageFeatureSet(scRNA, features=HB.genes) 
theme.set2 = theme(axis.title.x=element_blank())
plot.featrures = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb", "percent.HB")
group = "orig.ident"
# 质控前小提琴图
plots = list()
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(scRNA, group.by=group, pt.size = 0,
                       features = plot.featrures[i]) + theme.set2 + NoLegend()}
violin <- wrap_plots(plots = plots, nrow=2)    
# 设置质控指标
quantile(scRNA$nFeature_RNA, seq(0.01, 0.1, 0.01))
quantile(scRNA$nFeature_RNA, seq(0.9, 1, 0.01))
plots[[1]] + geom_hline(yintercept = 500) + geom_hline(yintercept = 8000)
quantile(scRNA$nCount_RNA, seq(0.9, 1, 0.01))
plots[[2]] + geom_hline(yintercept = 70000)
quantile(scRNA$percent.mt, seq(0.9, 1, 0.01))
plots[[3]] + geom_hline(yintercept = 20)
quantile(scRNA$percent.HB, seq(0.9, 1, 0.01))
plots[[5]] + geom_hline(yintercept = 1)
# 设置质控标准
minGene=500
maxGene=8000
maxUMI=70000
pctMT=20
pctHB=1
# 数据质控并绘制小提琴图
scRNA <- subset(scRNA, subset = nFeature_RNA > minGene & 
                  nFeature_RNA < maxGene & nCount_RNA < maxUMI &
                  percent.mt < pctMT & percent.HB < pctHB )
plots = list()
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(scRNA, group.by=group, pt.size = 0.01,
                       features = plot.featrures[i]) + theme.set2 + NoLegend()}
violin <- wrap_plots(plots = plots, nrow=2)    
ggsave("QC/CCvlnplot_after_qc.pdf", plot = violin, width = 9, height = 8) 

###### 细胞周期评分
scRNA <- NormalizeData(scRNA)  # 解决每个细胞测序深度不同的问题
g2m_genes <- cc.genes$g2m.genes
g2m_genes <- CaseMatch(search=g2m_genes, match=rownames(scRNA))
s_genes <- cc.genes$s.genes    
s_genes <- CaseMatch(search=s_genes, match=rownames(scRNA))
scRNA <- CellCycleScoring(scRNA, g2m.features=g2m_genes, s.features=s_genes)
library(pheatmap)
library(future)
library(tidyverse)
library(patchwork)
###### 数据标准化
scRNA <- NormalizeData(scRNA) %>% FindVariableFeatures(nfeatures = 2000)
scRNA <- ScaleData(scRNA, vars.to.regress = c("percent.mt","S.Score","G2M.Score"))
###### 降维聚类
scRNA <- RunPCA(scRNA, verbose = F)
########多样本整合
library(harmony)
scRNA <- RunHarmony(scRNA, group.by.vars="orig.ident", assay.use="RNA", max.iter.harmony = 30)
pcs = 1:20
library(clustree)
pbmc <- FindNeighbors(scRNA, reduction = "harmony", dims = pcs) 
pbmc <- FindClusters(pbmc,resolution=c(seq(.1,1.0,.1)))

scRNA <- FindNeighbors(scRNA, reduction = "harmony", dims = pcs) %>%FindClusters(resolution=0.7)
scRNA <- RunUMAP(scRNA, reduction = "harmony", dims = pcs) %>% RunTSNE(reduction = "harmony", dims = pcs)
p1 <- DimPlot(scRNA, reduction = "tsne", group.by="seurat_clusters",label = T)
p2 <- DimPlot(scRNA, reduction = "umap",  group.by="seurat_clusters",label = T)
p1 = p1 + ggsci::scale_color_igv()
p2 = p2 + ggsci::scale_color_igv()
p <- p1|p2
p1 <- DimPlot(scRNA, group.by = "Phase", label = T,reduction = "tsne")
p1 <- DimPlot(scRNA, group.by = "orig.ident", label = T,reduction = "tsne")
p2 <- DimPlot(scRNA, reduction = "tsne", label = T)
p <- p1|p2
p
############群体鉴定############
genes=c("EPCAM","CDH1",#Epithelial
        "PTPRC","ITGAM","LYZ","ITGAX","CSF1R",#Myeloid
        "CD68","CD14",#Macrophage/Mono
        "CD2","CD3D","CD3E","CD3G","CD4","CD8A",#T
        "CD19","MS4A1","CD79A","CD79B","BLK",#B
        "NCAM1","FCGR3A","KLRB1","KLRC1",#NK
        "COL1A1","COL1A2","COL3A1",
        "COL6A3","LUM","DCN","ACTA2","PDPN","FAP","SPARC","VIM","PDGFRB",#CAF
        "PECAM1","VWF","CDH5","ENG","CD34","FLT4","ICAM1","MCAM","SELE","VCAM1",#Endothelial
        "SDC1","SLAMF7","IGKC","IGLC2","IGHA1","IGHG3","MZB1")#Plasma
VlnPlot(scRNA, features =panmakerT, stack = T, flip = T,group.by = "seurat_clusters") + 
  NoLegend()
scRNA=subset(scRNA,celltype%in%"Tcells")
############figS1 A#########
DimPlot(scRNA, group.by = "group",reduction = "tsne", label = T,label.size = 6)
DimPlot(scRNA, group.by = "group",reduction = "tsne", label = T,label.size = 6)
############figS1 B#########
#小提琴图1-----------------------------------------
library(RColorBrewer)
library(knitr)
#remotes::install_github("lyc-1995/MySeuratWrappers") 
library(MySeuratWrappers)
library(ggplot2)
markers <- c("EPCAM","CDH1","CD2","CD3D","BLK",
             "MS4A1","IGKC","MZB1","LYZ","CD68",
             "COL1A1","PDGFRB","PECAM1","VWF")
my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', 
               '#D6E7A3', '#57C3F3', '#476D87',
               '#E95C59', '#E59CC4', '#AB3282', 
               '#BD956A',  '#585658','#F3B1A0', 
               '#9FA3A8', '#E0D4CA', '#5F3D69', 
               '#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#E39A35', 
               '#6778AE', '#91D0BE', '#B53E2B',
               '#712820', '#DCC1DD', '#CCE0F5',  
               '#CCC9E6', '#625D9E', '#68A180', 
               '#3A6963','#968175','#23452F','#C1E6F3')
VlnPlot(scRNA, features = markers,
        stacked=T,pt.size=0,
        cols = my36colors,
        direction = "vertical",
        x.lab = '', y.lab = '')+
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())
###### TRMs鉴定###########
library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)
panmakerT <- c("PTPRC", "CD2",  "CD3D", "CD3E","CD4","CD8A")
TRM=c("ITGAE","CD69")
FeaturePlot(scRNA, features = panmakerT, ncol = 3,
            cols = c("gray","red"),reduction = "tsne")
exp <- data.frame(FetchData(object = scRNA, vars = TRM))
exp$type=ifelse(exp$ITGAE>0.5,"Positive","Negative")
scRNA@meta.data$type=exp$type
scRNA@meta.data$type=ifelse(scRNA@meta.data$type%in%c("TRM"),"TRM","non-TRM")
#########fig1b############
DimPlot(scRNA, group.by = "group",reduction = "tsne", label = T,label.size = 6)
FeaturePlot(scRNA, features = TRM, ncol = 3,cols = c("#FEF1EB","#FF0019"),reduction = "tsne")
#########fig1c############
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db) ##加载人类
library(tidyverse)
library(ggrepel)
library(knitr)
library(RColorBrewer)

Idents(scRNA)="group"
deg<- FindMarkers(scRNA, ident.1 = "Tumor", min.pct = 0.1)
deg <- deg[which(deg$p_val<0.05),]
deg <- deg[!grepl("^RP[SL]", rownames(deg), ignore.case = F),]
deg <- deg[!grepl("^MT-", rownames(deg), ignore.case = F),]
colnames(deg)
deg$difference <- deg$pct.1 - deg$pct.2
deg_sig <- deg[which(deg$p_val<0.05 & abs(deg$avg_log2FC) >0.25),]
deg_sig$label <- rownames(deg_sig)
top20=deg_sig %>% top_n(30,abs(avg_log2FC))
colors <- rep(colors, each = 2)[1:20]
p1=ggplot(deg_sig, aes(x=difference, y=avg_log2FC)) + 
  geom_point(size=2, color="grey60") + 
  geom_text_repel(data = top20, aes(label=label),max.overlaps = 200000,
                  color="black",fontface="italic")+
  geom_point(data=deg_sig[which(deg_sig$p_val<0.05 & deg_sig$avg_log2FC>0.25),],
             aes(x=difference, y=avg_log2FC),
             size=2, color="#E41A1C",alpha=0.5)+
  geom_point(data=deg_sig[which(deg_sig$p_val<0.05 & deg_sig$avg_log2FC< -0.25),],
             aes(x=difference, y=avg_log2FC),
             size=2, color="#377EB8",alpha=0.5)+
  theme_classic()+
  theme(axis.text.x = element_text(colour = 'black',size = 12),
        axis.text.y = element_text(colour = 'black',size = 12),
        axis.title = element_text(colour = 'black',size = 15),
        axis.line = element_line(color = 'black', size = 1))+
  geom_hline(yintercept = 0,lty=2,lwd = 1)+
  geom_vline(xintercept = 0,lty=2,lwd = 1)+
  ylab("Log-fold Change")+
  xlab("Delta Percent")
##############fig1d##########
singlecell_gene_test <- function(SerautObj, 
                                 genes.use, 
                                 group.by=NULL, 
                                 assay = "RNA", 
                                 comp = NULL, 
                                 alpha_start = .05, 
                                 Bonferroni = T,
                                 only_postive =F) {
  p_val.out <- c()
  stat.out <- c()
  condition.out <- c()
  gene.out <- c()
  if (only_postive == F){
    for (gene in genes.use){
      group1_cellname = rownames(SerautObj@meta.data[SerautObj@meta.data[[group.by]] == comp[1],])
      group1_exp = SerautObj@assays[[assay]]@data[gene, group1_cellname] 
      
      group2_cellname = rownames(SerautObj@meta.data[SerautObj@meta.data[[group.by]] == comp[2],])
      group2_exp = SerautObj@assays[[assay]]@data[gene, group2_cellname]
      t_out = t.test(group1_exp, group2_exp)
      cond = paste(comp[1], comp[2], sep = "_")
      condition.out <- c(condition.out, cond)
      stat.out <- c(stat.out, t_out[["statistic"]])
      p_val.out <- c(p_val.out, t_out[["p.value"]])
      gene.out <- c(gene.out, gene)
    }
  }
  else{
    for (gene in genes.use){
      group1_cellname = rownames(SerautObj@meta.data[SerautObj@meta.data[[group.by]] == comp[1],])
      group1_exp = SerautObj@assays[[assay]]@data[gene, group1_cellname]
      group1_exp <- group1_exp[which(group1_exp>0)] 
      
      
      group2_cellname = rownames(SerautObj@meta.data[SerautObj@meta.data[[group.by]] == comp[2],])
      group2_exp = SerautObj@assays[[assay]]@data[gene, group2_cellname]
      group2_exp <- group2_exp[which(group2_exp>0)] 
      
      t_out = t.test(group1_exp, group2_exp)
      cond = paste(comp[1], comp[2], sep = "_")
      condition.out <- c(condition.out, cond)
      stat.out <- c(stat.out, t_out[["statistic"]])
      p_val.out <- c(p_val.out, t_out[["p.value"]])
      gene.out <- c(gene.out, gene)
    }
    
  }
  
  if (Bonferroni == T){
    new_alpha = alpha_start/(2*length(genes.use))
    cat(paste("\n", "P-value for significance: p <", new_alpha, "\n"))
    sig_out = p_val.out < new_alpha
    dfOUT<- data.frame(gene=gene.out, condition = condition.out, p_val = p_val.out, statistic = stat.out, significant = sig_out)
    
    dfOUT$sig = ifelse(dfOUT$p_val > 0.05, "ns",
                       ifelse(dfOUT$p_val > 0.01, '*',
                              ifelse(dfOUT$p_val > 0.001, "**", "****")))
    
  }
  
  else{
    dfOUT<- data.frame(gene=gene.out, condition = condition.out, p_val = p_val.out, statistic = stat.out)
    dfOUT$sig = ifelse(dfOUT$p_val > 0.05, "ns",
                       ifelse(dfOUT$p_val > 0.01, '*',
                              ifelse(dfOUT$p_val > 0.001, "**", "****")))
  }
  
  return(dfOUT)
}
gene= c('TGFB1','ITGB7','GZMB','PRF1','CXCL13','CCL5')
A <- singlecell_gene_test(scRNA, 
                          genes.use =gene,
                          group.by = 'group', 
                          comp = c("TRMs", "non−TRMs"))
anno_pvalue <- format(A$p_val, scientific = T,digits = 3) 
anno_sig <- A$sig

plots_violins <- VlnPlot(scRNA, 
                         cols = c("#E24840", "##1F78B4"),
                         pt.size = 0,
                         group.by = "group",
                         features = gene, 
                         ncol = 2, 
                         log = FALSE,
                         combine = FALSE)

for(i in 1:length(plots_violins)) {
  data <- plots_violins[[i]]$data
  colnames(data)[1] <- 'gene'
  plots_violins[[i]] <- plots_violins[[i]] + 
    theme_classic() + 
    theme(axis.text.x = element_text(size = 10,color="black"),
          axis.text.y = element_text(size = 10,color="black"),
          axis.title.y= element_text(size=12,color="black"),
          axis.title.x = element_blank(),
          legend.position='none')+
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+
    scale_x_discrete(labels = c("TRMs", "non−TRMs"))+
    geom_signif(annotations = anno_sig[i],
                y_position = max(data$gene)+0.5,
                xmin = 1,
                xmax = 2,
                tip_length = 0)
}
CombinePlots(plots_violins)
#############fig1e###################
library(Seurat)
library(patchwork)
library(clusterProfiler)
library(org.Hs.eg.db) 
library(tidyverse)
library(ggplot2)
library(stringr)
Idents(scRNA)<-"group"
Positive <- FindMarkers(scRNA, ident.1 = "TRMs", min.pct = 0.1)
sig_dge.all <- subset(Positive, p_val<0.05&abs(avg_log2FC)>0.25)
sig_dge.up <- subset(Positive, p_val<0.05&avg_log2FC>0.25)
sig_dge.up <- sig_dge.up[order(sig_dge.up$avg_log2FC,decreasing = T),]
up_go_BP <- enrichGO(gene          = row.names(sig_dge.up),
                     OrgDb         = 'org.Hs.eg.db',
                     keyType       = 'SYMBOL',
                     ont           = "ALL",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05) 
up_go_BP <- data.frame(up_go_BP)
all2<- up_go_BP[order(up_go_BP$Count,decreasing = T),]
all2$ratio<- all2$Count/525
num=c("0042110","0019221","0007159","0009615","0002253")
select=which(all2$ID%in%paste0("GO:",num))
all2=all2[select,]
all2<- all2[order(all2$ratio,decreasing = T),]
all2$Description<- factor(all2$Description,levels = rev(all2$Description))
ggplot(data = all2, aes(x = ratio, y = Description, fill = pvalue)) +
  scale_fill_gradientn(colors = c("#E83D7E", "#EF84AD","#EFABC5")) +
  geom_bar(stat = "identity", width = 0.9,color = "gray",size = 0) +
  theme_classic2() +
  labs(x = "Number of Gene", y = "", title = "GO TRMs vs Non-TRM")
##############fig1f################
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(Seurat)
library(patchwork)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
Idents(scRNA)<-"group"
deg<- FindMarkers(scRNA, ident.1 = "TRMs", min.pct = 0.05)
deg <- subset(deg, p_val<0.05 & abs(avg_log2FC)>0.25)
nrDEG=deg[,c(2,1)]
colnames(nrDEG)=c('log2FoldChange','pvalue')
head(nrDEG) 
gene <- bitr(rownames(nrDEG),    
             fromType = "SYMBOL",     
             toType =  "ENTREZID",   
             OrgDb = org.Hs.eg.db)   
gene$logFC <- nrDEG$log2FoldChange[match(gene$SYMBOL,rownames(nrDEG))] 
geneList=gene$logFC
names(geneList)=gene$ENTREZID 
geneList=sort(geneList,decreasing = T)
head(geneList)
go_gse <- gseGO(geneList     = geneList,
                OrgDb        = "org.Hs.eg.db",
                ont          = "BP",
                nPerm        = 1000,  
                minGSSize    = 10,
                maxGSSize    = 500,
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH",
                verbose      = FALSE)  
# Extract results and set gene IDs readable
go <- DOSE::setReadable(go_gse, OrgDb = 'org.Hs.eg.db', keyType = 'ENTREZID')
# Sort results by enrichment score
sortgo <- go[order(go$enrichmentScore, decreasing = TRUE),]
library(enrichplot)
library(GseaVis)
A1="GO:"
gseaNb(object = go_gse,
       geneSetID = A1,
       newGsea = T,
       addPval = T,
       rmHt = T,
       pvalX = 0.9,
       pvalY = 0.5,
       newCurveCol = c("#377EB8","white","#E41A1C"))
#############fig1g&h################
Singlecellratio_plotstat <- function (seu, by = "type",meta.include = NULL, 
                                      group_by = NULL, shape_by = NULL,
                                      custom_fill_colors = NULL, group_by.point = NULL, color_by = NULL, 
                                      pb = FALSE, comparisons = my_comparisons, 
                                      ncol = NULL, label = c("p.format","p.signif"), 
                                      label.x = NA, pt.size = NA) 
{
  by <- match.arg(by)  
  if (is.null(group_by)){ 
    group_by <- "null.group" 
  } 
  shapes <- NULL 
  if (!is.null(shape_by)) {
    shapes <- c(16,15,3,7,8,18,5,6,2,4,1,17)
  }
  fq <- prop.table(table(seu@meta.data$type, seu@meta.data[,"orig.ident"]), margin=2) *100
  df <- reshape2::melt(fq, value.name = "freq", varnames = c("type", "orig.ident"))  
  uniques <- apply(seu@meta.data, 2, function(x) length(unique(x))) 
  ei <- unique(seu@meta.data[, names(uniques[uniques<=100])])
  ei <- unique(ei[,colnames(ei) %in% meta.include])
  df <- merge(df, ei, by = "orig.ident")
  df <- cbind(df, null.group = paste("1"))
  df$orig.ident <- as.factor(df$orig.ident)
  if (is.null(x = ncol)) {
    ncol <- 3
    if (length(unique(df$celltype)) > 9) {
      ncol <- 4
    }
    if (length(unique(df$celltype)) > 20) {
      ncol <- 5
    }
  }
  custom_fill_colors = c(RColorBrewer::brewer.pal(9, "Oranges")[2], 
                         RColorBrewer::brewer.pal(9, "Reds")[6], 
                         RColorBrewer::brewer.pal(9, "Oranges")[5], 
                         RColorBrewer::brewer.pal(9, "Blues")[4:9])
  p <- ggplot(df, aes_string(y = "freq", x = group_by)) + labs(x = NULL,y = "Proportion (%)") + theme_bw() + theme(panel.grid.minor = element_blank(), 
                                                                                                                   panel.grid.major = element_blank(), strip.background = element_rect(fill = NA,color = NA), 
                                                                                                                   strip.text = element_text(face = "bold", size = 14), 
                                                                                                                   axis.ticks.x = element_blank(), axis.text = element_text(color = "black"), 
                                                                                                                   axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,color = 'black', size = 12),
                                                                                                                   axis.text.y = element_text(color = 'black', hjust = 1, vjust = 0.5, size = 12),
                                                                                                                   axis.title.y = element_text(color = 'black', size = 14))
  if(by=="type" && color_by=="type") {
    p + facet_wrap(group_by, scales = "free_x") + 
      geom_bar(aes_string(x = "orig.ident", fill = "factor(type)"), position = "fill", stat = "identity") + 
      scale_fill_manual("type", values = c("#FB8072", "#1965B0", "#7BAFDE", "#882E72","#B17BA6", 
                                           "#FF7F00", "#FDB462", "#E7298A", "#E78AC3","#33A02C", 
                                           "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D","#E6AB02")) + 
      scale_y_continuous(expand = c(0, 0), labels = seq(0, 100, 25)) + theme(panel.border = element_blank())
  }
  else {
    switch(by, type = p + facet_wrap("type", scales = "free_y", 
                                     ncol = ncol) + guides(fill = FALSE) + geom_boxplot(aes_string(x = group_by), 
                                                                                        alpha = 0.25, outlier.color = NA) + geom_point(size = 4, position = position_jitter(width = 0.25), 
                                                                                                                                       aes_string(x = group_by, y = "freq", color = color_by, 
                                                                                                                                                  shape = shape_by)) + scale_shape_manual(values = shapes) + 
             theme(panel.grid.major = element_line(color = "grey", 
                                                   size = 0.25)) + scale_color_manual(values = custom_fill_colors) + scale_fill_manual(values = custom_fill_colors)) + 
      scale_y_continuous(expand = expand_scale(mult = c(0.05, 0.2)))+
      ggpubr::stat_compare_means(mapping = aes_string(group_by), comparisons = comparisons, label = label,method = "t.test")
  }
}
my_comparisons <- list(c("Normal", "Tumor"))
Singlecellratio_plotstat(AT, group_by = "group",
                         meta.include = c("group","orig.ident"),
                         ncol = 2,
                         color_by = 'type')
Singlecellratio_plotstat(AT, group_by = "group",
                         meta.include = c("group","orig.ident"),
                         comparisons = my_comparisons, color_by = 'group',
                         group_by.point = "orig.ident",label.x = 1, pt.size = 3,
                         label = 'p.format', ncol =2)
############fig2a##############
DimPlot(scRNA, group.by = "group",reduction = "tsne", label = T,label.size = 6)
DimPlot(scRNA, split.by = "group",reduction = "tsne", label = T,label.size = 6)
############fig2b##############
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(Seurat)
library(patchwork)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(EnhancedVolcano)

Idents(scRNA)<-"group"
deg<- FindMarkers(scRNA, ident.1 = "Tumor", min.pct = 0.1)
colors=brewer.pal(11, 'Paired')
colors
group<-ifelse(
  deg$avg_log2FC<(-0.25)&deg$p_val_adj<0.05,'#3E7C17',
  ifelse(deg$avg_log2FC>(0.25)&deg$p_val_adj<0.05,'#DA1212',
         '#b5b5b5'))
group[is.na(group)]<-'#b5b5b5'
names(group)[group=='#DA1212']<-'Up'
names(group)[group=='#b5b5b5']<-'Nodiff'
names(group)[group=='#3E7C17']<-'Down'
p3=EnhancedVolcano(deg,
                   x="avg_log2FC",
                   y ="p_val_adj",
                   lab=rownames(deg),
                   FCcutoff=0.25,
                   pointSize=5,
                   labSize=7,
                   xlim=c(-4.5, 4.5),
                   ylim=c(0,300),
                   #colAlpha=0.4,
                   selectLab=c('CXCL13','HLA-DRB5','CD74','HLA-DRA','IFITM1',
                               "GZMB","HSPA6","XCL2","ANXA1","DNAJB1","FOS","XCL1"
                   ),
                   xlab=bquote(~Log[2]~'fold change'),
                   labCol='black',
                   labFace='bold',
                   #labSize=5,
                   boxedLabels=TRUE,
                   drawConnectors=T,
                   widthConnectors=0.8,
                   endsConnectors="last",
                   colConnectors='black',
                   colCustom=group,
                   colAlpha=0.55,
                   cutoffLineType='longdash',
                   cutoffLineCol='lightgrey',
                   cutoffLineWidth=0.7,
                   title="Volcano Plot",
                   subtitle="Differential expression")
############fig2c##############
library(pheatmap)
gene<- c("CCRT","TCF7","LEF1","SELL","PRF1", "GNLY", "NKG7", "GZMB", "GZMA", "GZMH", 
         "MKI67", "STMN1", "LAG3", "TIGIT", "PDCD1", "CTLA4", "HAVCR2")
Idents(scRNA)="orig.ident"
mat <- AverageExpression(scRNA, assays = "RNA", slot = "data", features = gene)[[1]]
anno <- data.frame(celltype=colnames(mat), row.names = colnames(mat))
pheatmap(mat, show_colnames = F, cluster_rows = F, cluster_cols = F, border_color = NA,
         color = colorRampPalette(c("#344CB7", "white", "#CD1818"))(50),
         annotation_col = anno, scale = "row",  fontsize = 2)
#############fig2d##############
library(Seurat)
library(patchwork)
library(clusterProfiler)
library(org.Hs.eg.db) 
library(tidyverse)
library(ggplot2)
library(stringr)
Idents(scRNA)<-"type"
Positive <- FindMarkers(scRNA, ident.1 = "Tumor", min.pct = 0.1)
sig_dge.all <- subset(Positive, p_val<0.05&abs(avg_log2FC)>0.25)
sig_dge.up <- subset(Positive, p_val<0.05&avg_log2FC>0.25)
sig_dge.up <- sig_dge.up[order(sig_dge.up$avg_log2FC,decreasing = T),]
up_go_BP <- enrichGO(gene          = row.names(sig_dge.up),
                     OrgDb         = 'org.Hs.eg.db',
                     keyType       = 'SYMBOL',
                     ont           = "ALL",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05) 
up_go_BP <- data.frame(up_go_BP)
all2<- up_go_BP[order(up_go_BP$Count,decreasing = T),]
all2$Count2<- all2$Count
num=c("00")
select=which(all2$ID%in%paste0("GO:",num))
all2=all2[select,]
all2<- all2[order(all2$Count2,decreasing = T),]
all2$Description<- factor(all2$Description,levels = rev(all2$Description))
ggplot(data = all2, aes(x = Count2, y = Description, fill = pvalue)) +
  scale_fill_gradientn(colors = c("#E83D7E", "#EF84AD","#EFABC5")) +
  geom_bar(stat = "identity", width = 0.9,color = "gray",size = 0) +
  theme_classic2() +
  labs(x = "Number of Gene", y = "", title = "GO Tumor vs Normal")
#############fig2e##############
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(Seurat)
library(patchwork)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
Idents(scRNA)<-"type"
deg<- FindMarkers(scRNA, ident.1 = "Tumor", min.pct = 0.05)
deg <- subset(deg, p_val<0.05 & abs(avg_log2FC)>0.25)
nrDEG=deg[,c(2,1)]
colnames(nrDEG)=c('log2FoldChange','pvalue')
head(nrDEG) 
gene <- bitr(rownames(nrDEG),    
             fromType = "SYMBOL",     
             toType =  "ENTREZID",   
             OrgDb = org.Hs.eg.db)   
gene$logFC <- nrDEG$log2FoldChange[match(gene$SYMBOL,rownames(nrDEG))] 
geneList=gene$logFC
names(geneList)=gene$ENTREZID 
geneList=sort(geneList,decreasing = T)
head(geneList)
go_gse <- gseGO(geneList     = geneList,
                OrgDb        = "org.Hs.eg.db",
                ont          = "BP",
                nPerm        = 1000,  
                minGSSize    = 10,
                maxGSSize    = 500,
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH",
                verbose      = FALSE)  
go <- DOSE::setReadable(go_gse, OrgDb = 'org.Hs.eg.db', keyType = 'ENTREZID')
sortgo <- go[order(go$enrichmentScore, decreasing = TRUE),]
library(enrichplot)
library(GseaVis)
A1="GO:"
gseaNb(object = go_gse,
       geneSetID = A1,
       newGsea = T,
       addPval = T,
       rmHt = T,
       pvalX = 0.9,
       pvalY = 0.5,
       newCurveCol = c("#377EB8","white","#E41A1C"))
############figs2#################
library(UCell)
library(Seurat)
library(GSEABase)
library(BiocParallel)
library(ggplot2)
library(remotes)
library(ggpubr)
library(gridExtra)

run_UCell <- function(rds,signatures,chunk.size=1000,assay=NULL,slot="data",seed=1,ncores=1){
  set.seed(seed)
  geneSets <- getGmt(signatures)
  markers <- geneIds(geneSets)
  rds1 <- AddModuleScore_UCell(rds, features = markers,assay=assay,slot=slot,ncores=ncores,chunk.size=chunk.size)
  metadata <- rds1@meta.data
  
  pathway <- c()
  for (i in names(markers)){
    pathway <- c(pathway,colnames(metadata)[grep(i,colnames(metadata))])
  }
  meta.merge <- FetchData(rds1,vars=pathway)
  colnames(meta.merge) <- gsub("_UCell","",colnames(meta.merge))
  rds <- Seurat::AddMetaData(rds, as.data.frame(meta.merge))
  return(rds)    
}

a="type"
pathways <- c("Activation of immune response", 
              "T cell proliferation", 
              "Cytokine mediated signaling pathway", 
              "T cell activation", 
              "Regulation of immune effector process", 
              "Positive regulation of cell activation")
plot_list <- list()
for (i in 1:length(pathways)) {
  Set <- pathways[i]
  signatures <- paste0(Set,".gmt")
  rds1 <- run_UCell(rds = scRNA, signatures = signatures, seed = 1, assay = "RNA", slot = "data")
  metadata <- rds1@meta.data
  rds1$type <- factor(rds1$type, ordered = TRUE, levels = c("Tumor", "Normal"))
  my_comparisons <- list(c("Tumor", "Normal"))
  p <- ggviolin(rds1@meta.data, x = a, y = Set, color = a, add = 'mean_sd', fill = a,
                add.params = list(color = "black")) + 
    stat_compare_means(comparisons = my_comparisons, label = "p.signif") + 
    scale_color_manual(values = my9color) + 
    scale_fill_manual(values = my9color) +
    theme(axis.text.x.bottom = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
    NoLegend() + labs(x = '')
  plot_list[[i]] <- p
}

grid.arrange(grobs = plot_list, ncol = 3)
############fig3a##############
DimPlot(scRNA, reduction = "tsne", label = T,label.size = 6)
markers=c("CD8A","CD4","CXCL13","PLAC8","ITM2C","STMN1","FCER1G","XIST")
FeaturePlot(scRNA, features = markers, ncol = 3,cols = c("#FEF1EB","#FF0019"),reduction = "tsne")
############fig3b&c##############
library(pheatmap)
gene1<- c("CXCL13","KRT86","TRBV2","FGFBP2","FCGR3A", "PLAC8", "GZMK", "ITM2C", "CCL4L2", "STMN1", 
          "TUBA1B", "TUBB", "TYROBP", "FCER1G", "GNLY", "MALAT1", "XIST","S100A9","TNFRSF4","LTB","FOXP3")
gene2<- c("CCRT","TCF7","LEF1","SELL","CD28", "ICOS", "TNFRSF9", "TNFRSF14", "GZMA", "GZMB", 
          "GZMK", "GNLY", "NKG7", "PRF1", "LAG3", "TIGIT", "PDCD1", "HAVCR2", "CTLA4")
Idents(scRNA)="celltype"
mat <- AverageExpression(scRNA, assays = "RNA", slot = "data", features = gene1)[[1]]
anno <- data.frame(celltype=colnames(mat), row.names = colnames(mat))
pheatmap(mat, show_colnames = F, cluster_rows = F, cluster_cols = F, border_color = NA,
         color = colorRampPalette(c("#344CB7", "white", "#CD1818"))(50),
         annotation_col = anno, scale = "row",  fontsize = 2)
mat <- AverageExpression(scRNA, assays = "RNA", slot = "data", features = gene2)[[1]]
anno <- data.frame(celltype=colnames(mat), row.names = colnames(mat))
pheatmap(mat, show_colnames = F, cluster_rows = F, cluster_cols = F, border_color = NA,
         color = colorRampPalette(c("#344CB7", "white", "#CD1818"))(50),
         annotation_col = anno, scale = "row",  fontsize = 2)
############fig3d################
Singlecellratio_plotstat <- function (seu, by = "cell.type",meta.include = NULL, 
                                      group_by = NULL, shape_by = NULL,
                                      custom_fill_colors = NULL, group_by.point = NULL, color_by = NULL, 
                                      pb = FALSE, comparisons = my_comparisons, 
                                      ncol = NULL, label = c("p.format","p.signif"), 
                                      label.x = NA, pt.size = NA) 
{
  by <- match.arg(by)  
  if (is.null(group_by)){ 
    group_by <- "null.group" 
  } 
  shapes <- NULL 
  if (!is.null(shape_by)) {
    shapes <- c(16,15,3,7,8,18,5,6,2,4,1,17)
  }
  fq <- prop.table(table(seu@meta.data$celltype, seu@meta.data[,"orig.ident"]), margin=2) *100
  df <- reshape2::melt(fq, value.name = "freq", varnames = c("cell.type", "orig.ident"))  
  uniques <- apply(seu@meta.data, 2, function(x) length(unique(x))) 
  ei <- unique(seu@meta.data[, names(uniques[uniques<=100])])
  ei <- unique(ei[,colnames(ei) %in% meta.include])
  df <- merge(df, ei, by = "orig.ident")
  df <- cbind(df, null.group = paste("1"))
  df$orig.ident <- as.factor(df$orig.ident)
  if (is.null(x = ncol)) {
    ncol <- 3
    if (length(unique(df$celltype)) > 9) {
      ncol <- 4
    }
    if (length(unique(df$celltype)) > 20) {
      ncol <- 5
    }
  }
  custom_fill_colors = c(RColorBrewer::brewer.pal(9, "Oranges")[2], 
                         RColorBrewer::brewer.pal(9, "Reds")[6], 
                         RColorBrewer::brewer.pal(9, "Oranges")[5], 
                         RColorBrewer::brewer.pal(9, "Blues")[4:9])
  p <- ggplot(df, aes_string(y = "freq", x = group_by)) + labs(x = NULL,y = "Proportion (%)") + theme_bw() + theme(panel.grid.minor = element_blank(), 
                                                                                                                   panel.grid.major = element_blank(), strip.background = element_rect(fill = NA,color = NA), 
                                                                                                                   strip.text = element_text(face = "bold", size = 14), 
                                                                                                                   axis.ticks.x = element_blank(), axis.text = element_text(color = "black"), 
                                                                                                                   axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,color = 'black', size = 12),
                                                                                                                   axis.text.y = element_text(color = 'black', hjust = 1, vjust = 0.5, size = 12),
                                                                                                                   axis.title.y = element_text(color = 'black', size = 14))
  if(by=="cell.type" && color_by=="cell.type") {
    p + facet_wrap(group_by, scales = "free_x") + 
      geom_bar(aes_string(x = "orig.ident", fill = "factor(cell.type)"), position = "fill", stat = "identity") + 
      scale_fill_manual("cell.type", values = c("#FB8072", "#1965B0", "#7BAFDE", "#882E72","#B17BA6", 
                                                "#FF7F00", "#FDB462", "#E7298A", "#E78AC3","#33A02C", 
                                                "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D","#E6AB02")) + 
      scale_y_continuous(expand = c(0, 0), labels = seq(0, 100, 25)) + theme(panel.border = element_blank())
  }
  else {
    switch(by, cell.type = p + facet_wrap("cell.type", scales = "free_y", 
                                          ncol = ncol) + guides(fill = FALSE) + geom_boxplot(aes_string(x = group_by), 
                                                                                             alpha = 0.25, outlier.color = NA) + geom_point(size = 4, position = position_jitter(width = 0.25), 
                                                                                                                                            aes_string(x = group_by, y = "freq", color = color_by, 
                                                                                                                                                       shape = shape_by)) + scale_shape_manual(values = shapes) + 
             theme(panel.grid.major = element_line(color = "grey", 
                                                   size = 0.25)) + scale_color_manual(values = custom_fill_colors) + scale_fill_manual(values = custom_fill_colors)) + 
      scale_y_continuous(expand = expand_scale(mult = c(0.05, 0.2)))+
      ggpubr::stat_compare_means(mapping = aes_string(group_by), comparisons = comparisons, label = label,method = "wilcox.test")
  }
}
my_comparisons <- list(c("Normal", "Tumor"))
Singlecellratio_plotstat(scRNA, group_by = "group",
                         meta.include = c("group","orig.ident"),
                         ncol = 2,
                         color_by = 'cell.type')
Singlecellratio_plotstat(scRNA, group_by = "group",
                         meta.include = c("group","orig.ident"),
                         comparisons = my_comparisons, color_by = 'group',
                         group_by.point = "orig.ident",label.x = 1, pt.size = 3,
                         label = 'p.format', ncol =3)
############figs3a##############
library(scRNAtoolVis)
library(RColorBrewer)
library(knitr)
Idents(scRNA)="celltype"
ClusterMarker <- FindAllMarkers(scRNA, assay = "RNA",
                                slot = "data",
                                logfc.threshold = 0.25,
                                only.pos = T)
ClusterMarker1 <- ClusterMarker1[!grepl("^RP[SL]", ClusterMarker1$gene, ignore.case = F),]
ClusterMarker1 <- ClusterMarker1[!grepl("^MT-", ClusterMarker1$gene, ignore.case = F),]
colour=brewer.pal(11, 'Paired')[c(1:11)]
jjVolcano(diffData = ClusterMarker1,
          log2FC.cutoff = 0.25,
          fontface = 'italic',
          topGeneN=10,
          tile.col = colour)
############figs3b##############
Idents(scRNA)="celltype"
signature<- c("Naive", 
              "Inhibitory", 
              "Cytoxicity", 
              "T cell activation", 
              "Costimulatory")
a="celltype"
plot_list <- list()
for (i in 1:length(signature)) {
  Set <- signature[i]
  signatures <- paste0(Set,".gmt")
  rds1 <- run_UCell(rds = scRNA, signatures = signatures, seed = 1, assay = "RNA", slot = "data")
  metadata <- rds1@meta.data
  p <-VlnPlot(scRNA,features=Set,pt.size=0)
  plot_list[[i]] <- p
}
grid.arrange(grobs = plot_list, ncol = 2)
############fig4a##############
Idents(scRNA)="CXCL13TRM"
DimPlot(scRNA, reduction = "tsne", label = T,label.size = 6)
############fig4b##############
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(ggrepel)
library(knitr)
library(RColorBrewer)
deg<- FindMarkers(scRNA, ident.1 = "CXCL13TRM", min.pct = 0.1)
deg <- deg[which(deg$p_val<0.05),]
deg <- deg[!grepl("^RP[SL]", rownames(deg), ignore.case = F),]
deg <- deg[!grepl("^MT-", rownames(deg), ignore.case = F),]
deg$difference <- deg$pct.1 - deg$pct.2
deg_sig <- deg[which(deg$p_val<0.05 & abs(deg$avg_log2FC) >0.25),]
deg_sig$label <- rownames(deg_sig)
top20=deg_sig %>% top_n(30,abs(avg_log2FC))
colors <- rep(colors, each = 2)[1:20]
p1=ggplot(deg_sig, aes(x=difference, y=avg_log2FC)) + 
  geom_point(size=2, color="grey60") + 
  geom_text_repel(data = top20, aes(label=label),max.overlaps = 200000,
                  color="black",fontface="italic")+
  geom_point(data=deg_sig[which(deg_sig$p_val<0.05 & deg_sig$avg_log2FC>0.25),],
             aes(x=difference, y=avg_log2FC),
             size=2, color="#E41A1C",alpha=0.5)+
  geom_point(data=deg_sig[which(deg_sig$p_val<0.05 & deg_sig$avg_log2FC< -0.25),],
             aes(x=difference, y=avg_log2FC),
             size=2, color="#377EB8",alpha=0.5)+
  theme_classic()+
  theme(axis.text.x = element_text(colour = 'black',size = 12),
        axis.text.y = element_text(colour = 'black',size = 12),
        axis.title = element_text(colour = 'black',size = 15),
        axis.line = element_line(color = 'black', size = 1))+
  geom_hline(yintercept = 0,lty=2,lwd = 1)+
  geom_vline(xintercept = 0,lty=2,lwd = 1)+
  ylab("Log-fold Change")+
  xlab("Delta Percent")
############fig4c##############
Idents(scRNA)<-"group"
Positive <- FindMarkers(scRNA, ident.1 = "TRMs", min.pct = 0.1)
sig_dge.all <- subset(Positive, p_val<0.05&abs(avg_log2FC)>0.25)
sig_dge.up <- subset(Positive, p_val<0.05&avg_log2FC>0.25)
sig_dge.up <- sig_dge.up[order(sig_dge.up$avg_log2FC,decreasing = T),]
up_go_BP <- enrichGO(gene          = row.names(sig_dge.up),
                     OrgDb         = 'org.Hs.eg.db',
                     keyType       = 'SYMBOL',
                     ont           = "ALL",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05) 
up_go_BP <- data.frame(up_go_BP)
all2<- up_go_BP[order(up_go_BP$Count,decreasing = T),]
all2$Count<- all2$Count
num=c("")
select=which(all2$ID%in%paste0("GO:",num))
all2=all2[select,]
all2<- all2[order(all2$Count,decreasing = T),]
all2$Description<- factor(all2$Description,levels = rev(all2$Description))
ggplot(data = all2, aes(x = Count, y = Description, fill = pvalue)) +
  scale_fill_gradientn(colors = c("#E83D7E", "#EF84AD","#EFABC5")) +
  geom_bar(stat = "identity", width = 0.9,color = "gray",size = 0) +
  theme_classic2() +
  labs(x = "Number of Gene", y = "", title = "CXCL13-TRM vs Others-TRM")
############fig4d##############
table(scRNA$CXCL13TRM)
gene=c("IFI6","BST2","PLSCR1","OAS1","HERC5",
       "IFIT3","MX1","ISG15","STAT1","IRF7")
A <- singlecell_gene_test(scRNA, 
                          genes.use =gene,
                          group.by = 'CXCL13TRM', 
                          comp = c("CXCL13TRM", "Others-TRM"))
anno_pvalue <- format(A$p_val, scientific = T,digits = 3) 
anno_sig <- A$sig
plots_violins <- VlnPlot(scRNA, 
                         cols = c("#E24840", "#1F78B4"),
                         pt.size = 0,
                         group.by = "CXCL13TRM",
                         features = gene, 
                         ncol = 2, 
                         log = FALSE,
                         combine = FALSE)
for(i in 1:length(plots_violins)) {
  data <- plots_violins[[i]]$data
  colnames(data)[1] <- 'gene'
  plots_violins[[i]] <- plots_violins[[i]] + 
    theme_classic() + 
    theme(axis.text.x = element_text(size = 10,color="black"),
          axis.text.y = element_text(size = 10,color="black"),
          axis.title.y= element_text(size=12,color="black"),
          axis.title.x = element_blank(),
          legend.position='none')+
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+
    scale_x_discrete(labels = c("CXCL13TRM", "Others-TRM"))+
    geom_signif(annotations = anno_sig[i],
                y_position = max(data$gene)+0.5,
                xmin = 1,
                xmax = 2,
                tip_length = 0)
}
CombinePlots(plots_violins)
############fig4e##############
deg<- FindMarkers(scRNA, ident.1 = "CXCL13TRM", min.pct = 0.05)
deg <- subset(deg, p_val<0.05 & abs(avg_log2FC)>0.25)
nrDEG=deg[,c(2,1)]
colnames(nrDEG)=c('log2FoldChange','pvalue')
head(nrDEG) 
gene <- bitr(rownames(nrDEG),    
             fromType = "SYMBOL",     
             toType =  "ENTREZID",   
             OrgDb = org.Hs.eg.db)   
gene$logFC <- nrDEG$log2FoldChange[match(gene$SYMBOL,rownames(nrDEG))] 
geneList=gene$logFC
names(geneList)=gene$ENTREZID 
geneList=sort(geneList,decreasing = T)
head(geneList)
go_gse <- gseGO(geneList     = geneList,
                OrgDb        = "org.Hs.eg.db",
                ont          = "BP",
                nPerm        = 1000,  
                minGSSize    = 10,
                maxGSSize    = 500,
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH",
                verbose      = FALSE)  
go <- DOSE::setReadable(go_gse, OrgDb = 'org.Hs.eg.db', keyType = 'ENTREZID')
sortgo <- go[order(go$enrichmentScore, decreasing = TRUE),]
library(enrichplot)
library(GseaVis)
A1="GO:"
gseaNb(object = go_gse,
       geneSetID = A1,
       newGsea = T,
       addPval = T,
       rmHt = T,
       pvalX = 0.9,
       pvalY = 0.5,
       newCurveCol = c("#377EB8","white","#E41A1C"))
############fig4f##############
deg<- FindMarkers(scRNA, ident.1 = "Tumor", min.pct = 0.1)
deg <- deg[which(deg$p_val<0.05),]
deg <- deg[!grepl("^RP[SL]", rownames(deg), ignore.case = F),]
deg <- deg[!grepl("^MT-", rownames(deg), ignore.case = F),]
deg$difference <- deg$pct.1 - deg$pct.2
deg_sig <- deg[which(deg$p_val<0.05 & abs(deg$avg_log2FC) >0.25),]
deg_sig$label <- rownames(deg_sig)
top20=deg_sig %>% top_n(30,abs(avg_log2FC))
colors <- rep(colors, each = 2)[1:20]
p1=ggplot(deg_sig, aes(x=difference, y=avg_log2FC)) + 
  geom_point(size=2, color="grey60") + 
  geom_text_repel(data = top20, aes(label=label),max.overlaps = 200000,
                  color="black",fontface="italic")+
  geom_point(data=deg_sig[which(deg_sig$p_val<0.05 & deg_sig$avg_log2FC>0.25),],
             aes(x=difference, y=avg_log2FC),
             size=2, color="#E41A1C",alpha=0.5)+
  geom_point(data=deg_sig[which(deg_sig$p_val<0.05 & deg_sig$avg_log2FC< -0.25),],
             aes(x=difference, y=avg_log2FC),
             size=2, color="#377EB8",alpha=0.5)+
  theme_classic()+
  theme(axis.text.x = element_text(colour = 'black',size = 12),
        axis.text.y = element_text(colour = 'black',size = 12),
        axis.title = element_text(colour = 'black',size = 15),
        axis.line = element_line(color = 'black', size = 1))+
  geom_hline(yintercept = 0,lty=2,lwd = 1)+
  geom_vline(xintercept = 0,lty=2,lwd = 1)+
  ylab("Log-fold Change")+
  xlab("Delta Percent")
############fig4h##############
Positive <- FindMarkers(scRNA, ident.1 = "Tumor", min.pct = 0.1)
sig_dge.all <- subset(Positive, p_val<0.05&abs(avg_log2FC)>0.25)
sig_dge.up <- subset(Positive, p_val<0.05&avg_log2FC>0.25)
sig_dge.up <- sig_dge.up[order(sig_dge.up$avg_log2FC,decreasing = T),]
up_go_BP <- enrichGO(gene          = row.names(sig_dge.up),
                     OrgDb         = 'org.Hs.eg.db',
                     keyType       = 'SYMBOL',
                     ont           = "ALL",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05) 
up_go_BP <- data.frame(up_go_BP)
all2<- up_go_BP[order(up_go_BP$Count,decreasing = T),]
all2$Count<- all2$Count
num=c("")
select=which(all2$ID%in%paste0("GO:",num))
all2=all2[select,]
all2<- all2[order(all2$Count,decreasing = T),]
all2$Description<- factor(all2$Description,levels = rev(all2$Description))
ggplot(data = all2, aes(x = Count, y = Description, fill = pvalue)) +
  scale_fill_gradientn(colors = c("#E83D7E", "#EF84AD","#EFABC5")) +
  geom_bar(stat = "identity", width = 0.9,color = "gray",size = 0) +
  theme_classic2() +
  labs(x = "Number of Gene", y = "", title = "GO of CXCL13+CD8+TRMs in tumor")
############fig4i##############
deg<- FindMarkers(scRNA, ident.1 = "Tumor", min.pct = 0.05)
deg <- subset(deg, p_val<0.05 & abs(avg_log2FC)>0.25)
nrDEG=deg[,c(2,1)]
colnames(nrDEG)=c('log2FoldChange','pvalue')
head(nrDEG) 
gene <- bitr(rownames(nrDEG),    
             fromType = "SYMBOL",     
             toType =  "ENTREZID",   
             OrgDb = org.Hs.eg.db)   
gene$logFC <- nrDEG$log2FoldChange[match(gene$SYMBOL,rownames(nrDEG))] 
geneList=gene$logFC
names(geneList)=gene$ENTREZID 
geneList=sort(geneList,decreasing = T)
head(geneList)
go_gse <- gseGO(geneList     = geneList,
                OrgDb        = "org.Hs.eg.db",
                ont          = "BP",
                nPerm        = 1000,  
                minGSSize    = 10,
                maxGSSize    = 500,
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH",
                verbose      = FALSE)  
go <- DOSE::setReadable(go_gse, OrgDb = 'org.Hs.eg.db', keyType = 'ENTREZID')
sortgo <- go[order(go$enrichmentScore, decreasing = TRUE),]
A1="GO:"
gseaNb(object = go_gse,
       geneSetID = A1,
       newGsea = T,
       addPval = T,
       rmHt = T,
       pvalX = 0.9,
       pvalY = 0.5,
       newCurveCol = c("#377EB8","white","#E41A1C"))
############fig4g##############
gene=c("HLA−DRB5","HLA−DRA","CCL4L2","CCL5")
A <- singlecell_gene_test(scRNA, 
                          genes.use =gene,
                          group.by = 'group', 
                          comp = c("Tumor", "Normal"))
anno_pvalue <- format(A$p_val, scientific = T,digits = 3) 
anno_sig <- A$sig
plots_violins <- VlnPlot(scRNA, 
                         cols = c("#E24840", "#1F78B4"),
                         pt.size = 0,
                         group.by = "group",
                         features = gene, 
                         ncol = 2, 
                         log = FALSE,
                         combine = FALSE)
for(i in 1:length(plots_violins)) {
  data <- plots_violins[[i]]$data
  colnames(data)[1] <- 'gene'
  plots_violins[[i]] <- plots_violins[[i]] + 
    theme_classic() + 
    theme(axis.text.x = element_text(size = 10,color="black"),
          axis.text.y = element_text(size = 10,color="black"),
          axis.title.y= element_text(size=12,color="black"),
          axis.title.x = element_blank(),
          legend.position='none')+
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+
    scale_x_discrete(labels = c("Tumor", "Normal"))+
    geom_signif(annotations = anno_sig[i],
                y_position = max(data$gene)+0.5,
                xmin = 1,
                xmax = 2,
                tip_length = 0)
}
CombinePlots(plots_violins)
############fig4j&5J##############
###data from TCGA database#####
library(limma)
library(GSVA)
library(GSEABase) 
library(clusterProfiler)
library(org.Hs.eg.db)
library(data.table)
library(dplyr)
library(GSVA)
library(GSEABase) 
exprSet=read.table("CCsymbol.txt", header=T, sep="\t",row.names=1, check.names=F)
geneSet=getGmt("CXCL13TRM marker.gmt", geneIdType=SymbolIdentifier())
exprSet=as.matrix(exprSet)
ssgseaScore=gsva(exprSet, geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
write.table(ssgseaScore,"CXCL13TRM marker score.txt",sep = "\t")

input1="TCGA_CESC_Sur.txt"
input2="CXCL13TRM marker score.txt"
clinical=read.table(input1, header=T, sep="\t") 
CXCL13TRMscore=read.table(input2, header=T, sep="\t") 
a=merge(clinical,CXCL13TRMscore,by = 'id')
write.table (a, file ="expTime_CC.csv", sep ="\t",row.names = F)

svdata <- read.csv("expTime_CC.csv",header = T,row.names = 1)
library(survival)
library(survminer)
res.cut <- surv_cutpoint(svdata, time = "futime", 
                         event = "fustat", 
                         variables = names(svdata)[3:ncol(svdata)], 
                         minprop = 0.3) 
res.cat <- surv_categorize(res.cut)
my.surv <- Surv(res.cat$futime, res.cat$fustat)
pl<-list()
for (i in colnames(res.cat)[3:ncol(svdata)]) {
  group <- res.cat[,i] 
  survival_dat <- data.frame(group = group)
  fit <- survfit(my.surv ~ group)
  group <- factor(group, levels = c("low", "high"))
  data.survdiff <- survdiff(my.surv ~ group)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  
  HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
  CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
  svsort <- svdata[order(svdata[,i]),]
  pl[[i]]<-ggsurvplot(fit, data = survival_dat ,
                      conf.int = F, 
                      censor = T, 
                      palette = c("#ad1d32","#0c7bac"), 
                      legend.title = i,
                      font.legend = 11,
                      legend.labs=c(paste0(">",round(svsort[fit$n[2],i],2),"(",fit$n[1],")"),
                                    paste0("<",round(svsort[fit$n[2],i],2),"(",fit$n[2],")")),
                      pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                                 paste("p = ",round(p.val,3), sep = "")),
                                   HR, CI, sep = "\n"))
  
}
res <- arrange_ggsurvplots(pl,print = T,ncol = 2, nrow = 1)
############figs4a##############
a="celltype"
pathways <- c("Defense response to virus", 
              "Response to type I interferon")
plot_list <- list()
for (i in 1:length(pathways)) {
  Set <- pathways[i]
  signatures <- paste0(Set,".gmt")
  rds1 <- run_UCell(rds = scRNA, signatures = signatures, seed = 1, assay = "RNA", slot = "data")
  metadata <- rds1@meta.data
  p <- ggviolin(rds1@meta.data, x = a, y = Set, color = a, add = 'mean_sd', fill = a,
                add.params = list(color = "black")) + 
    theme(axis.text.x.bottom = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
    NoLegend() + labs(x = '')
  plot_list[[i]] <- p
}
grid.arrange(grobs = plot_list, ncol = 3)
############figs4b##############
a="group"
pathways <- c("Defense response to virus", 
              "Response to type I interferon")
plot_list <- list()
for (i in 1:length(pathways)) {
  p <- ggviolin(rds1@meta.data, x = a, y = Set, color = a, add = 'mean_sd', fill = a,
                add.params = list(color = "black")) + 
    stat_compare_means(comparisons = my_comparisons, label = "p.signif") + 
    scale_color_manual(values = my9color) + 
    scale_fill_manual(values = my9color) +
    theme(axis.text.x.bottom = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
    NoLegend() + labs(x = '')
  plot_list[[i]] <- p
}
grid.arrange(grobs = plot_list, ncol = 3)

############fig5a##############
Idents(scRNA)="PLAC8TRM"
DimPlot(scRNA, reduction = "tsne", label = T,label.size = 6)
############fig5b##############
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(ggrepel)
library(knitr)
library(RColorBrewer)
deg<- FindMarkers(scRNA, ident.1 = "PLAC8TRM", min.pct = 0.1)
deg <- deg[which(deg$p_val<0.05),]
deg <- deg[!grepl("^RP[SL]", rownames(deg), ignore.case = F),]
deg <- deg[!grepl("^MT-", rownames(deg), ignore.case = F),]
deg$difference <- deg$pct.1 - deg$pct.2
deg_sig <- deg[which(deg$p_val<0.05 & abs(deg$avg_log2FC) >0.25),]
deg_sig$label <- rownames(deg_sig)
top20=deg_sig %>% top_n(30,abs(avg_log2FC))
colors <- rep(colors, each = 2)[1:20]
p1=ggplot(deg_sig, aes(x=difference, y=avg_log2FC)) + 
  geom_point(size=2, color="grey60") + 
  geom_text_repel(data = top20, aes(label=label),max.overlaps = 200000,
                  color="black",fontface="italic")+
  geom_point(data=deg_sig[which(deg_sig$p_val<0.05 & deg_sig$avg_log2FC>0.25),],
             aes(x=difference, y=avg_log2FC),
             size=2, color="#E41A1C",alpha=0.5)+
  geom_point(data=deg_sig[which(deg_sig$p_val<0.05 & deg_sig$avg_log2FC< -0.25),],
             aes(x=difference, y=avg_log2FC),
             size=2, color="#377EB8",alpha=0.5)+
  theme_classic()+
  theme(axis.text.x = element_text(colour = 'black',size = 12),
        axis.text.y = element_text(colour = 'black',size = 12),
        axis.title = element_text(colour = 'black',size = 15),
        axis.line = element_line(color = 'black', size = 1))+
  geom_hline(yintercept = 0,lty=2,lwd = 1)+
  geom_vline(xintercept = 0,lty=2,lwd = 1)+
  ylab("Log-fold Change")+
  xlab("Delta Percent")
############fig5c##############
table(scRNA$PLAC8TRM)
gene=c("NKG7","KLRD1","KLRF1","KLRG1","GNLY",
       "GZMM","GZMH","PRF1","PDCD1","CTLA4","TIGIT","HAVCR2")
A <- singlecell_gene_test(scRNA, 
                          genes.use =gene,
                          group.by = 'PLAC8TRM', 
                          comp = c("PLAC8TRM", "Others-TRM"))
anno_pvalue <- format(A$p_val, scientific = T,digits = 3) 
anno_sig <- A$sig
plots_violins <- VlnPlot(scRNA, 
                         cols = c("#E24840", "#1F78B4"),
                         pt.size = 0,
                         group.by = "PLAC8TRM",
                         features = gene, 
                         ncol = 2, 
                         log = FALSE,
                         combine = FALSE)
for(i in 1:length(plots_violins)) {
  data <- plots_violins[[i]]$data
  colnames(data)[1] <- 'gene'
  plots_violins[[i]] <- plots_violins[[i]] + 
    theme_classic() + 
    theme(axis.text.x = element_text(size = 10,color="black"),
          axis.text.y = element_text(size = 10,color="black"),
          axis.title.y= element_text(size=12,color="black"),
          axis.title.x = element_blank(),
          legend.position='none')+
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+
    scale_x_discrete(labels = c("PLAC8TRM", "Others-TRM"))+
    geom_signif(annotations = anno_sig[i],
                y_position = max(data$gene)+0.5,
                xmin = 1,
                xmax = 2,
                tip_length = 0)
}
CombinePlots(plots_violins)
############fig5d##############
Idents(scRNA)<-"group"
Positive <- FindMarkers(scRNA, ident.1 = "TRMs", min.pct = 0.1)
sig_dge.all <- subset(Positive, p_val<0.05&abs(avg_log2FC)>0.25)
sig_dge.up <- subset(Positive, p_val<0.05&avg_log2FC>0.25)
sig_dge.up <- sig_dge.up[order(sig_dge.up$avg_log2FC,decreasing = T),]
up_go_BP <- enrichGO(gene          = row.names(sig_dge.up),
                     OrgDb         = 'org.Hs.eg.db',
                     keyType       = 'SYMBOL',
                     ont           = "ALL",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05) 
up_go_BP <- data.frame(up_go_BP)
all2<- up_go_BP[order(up_go_BP$Count,decreasing = T),]
all2$Count<- all2$Count
num=c("")
select=which(all2$ID%in%paste0("GO:",num))
all2=all2[select,]
all2<- all2[order(all2$Count,decreasing = T),]
all2$Description<- factor(all2$Description,levels = rev(all2$Description))
ggplot(data = all2, aes(x = Count, y = Description, fill = pvalue)) +
  scale_fill_gradientn(colors = c("#E83D7E", "#EF84AD","#EFABC5")) +
  geom_bar(stat = "identity", width = 0.9,color = "gray",size = 0) +
  theme_classic2() +
  labs(x = "Number of Gene", y = "", title = "PLAC8-TRM vs Others-TRM")
############fig5e##############
a="PLAC8TRM"
pathways <- c("Cytoxicity", 
              "Cell activation involved in immune response", 
              "Ribosome biogenesis", 
              "Cytosolic ribosome")
plot_list <- list()
for (i in 1:length(pathways)) {
  Set <- pathways[i]
  signatures <- paste0(Set,".gmt")
  rds1 <- run_UCell(rds = scRNA, signatures = signatures, seed = 1, assay = "RNA", slot = "data")
  metadata <- rds1@meta.data
  rds1$type <- factor(rds1$type, ordered = TRUE, levels = c("PLAC8-TRM", "Others-TRM"))
  my_comparisons <- list(c("PLAC8-TRM", "Others-TRM"))
  p <- ggviolin(rds1@meta.data, x = a, y = Set, color = a, add = 'mean_sd', fill = a,
                add.params = list(color = "black")) + 
    stat_compare_means(comparisons = my_comparisons, label = "p.signif") + 
    scale_color_manual(values = my9color) + 
    scale_fill_manual(values = my9color) +
    theme(axis.text.x.bottom = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
    NoLegend() + labs(x = '')
  plot_list[[i]] <- p
}
grid.arrange(grobs = plot_list, ncol = 3)
############fig5f##############
deg<- FindMarkers(scRNA, ident.1 = "Tumor", min.pct = 0.1)
deg <- deg[which(deg$p_val<0.05),]
deg <- deg[!grepl("^RP[SL]", rownames(deg), ignore.case = F),]
deg <- deg[!grepl("^MT-", rownames(deg), ignore.case = F),]
deg$difference <- deg$pct.1 - deg$pct.2
deg_sig <- deg[which(deg$p_val<0.05 & abs(deg$avg_log2FC) >0.25),]
deg_sig$label <- rownames(deg_sig)
top20=deg_sig %>% top_n(30,abs(avg_log2FC))
colors <- rep(colors, each = 2)[1:20]
p1=ggplot(deg_sig, aes(x=difference, y=avg_log2FC)) + 
  geom_point(size=2, color="grey60") + 
  geom_text_repel(data = top20, aes(label=label),max.overlaps = 200000,
                  color="black",fontface="italic")+
  geom_point(data=deg_sig[which(deg_sig$p_val<0.05 & deg_sig$avg_log2FC>0.25),],
             aes(x=difference, y=avg_log2FC),
             size=2, color="#E41A1C",alpha=0.5)+
  geom_point(data=deg_sig[which(deg_sig$p_val<0.05 & deg_sig$avg_log2FC< -0.25),],
             aes(x=difference, y=avg_log2FC),
             size=2, color="#377EB8",alpha=0.5)+
  theme_classic()+
  theme(axis.text.x = element_text(colour = 'black',size = 12),
        axis.text.y = element_text(colour = 'black',size = 12),
        axis.title = element_text(colour = 'black',size = 15),
        axis.line = element_line(color = 'black', size = 1))+
  geom_hline(yintercept = 0,lty=2,lwd = 1)+
  geom_vline(xintercept = 0,lty=2,lwd = 1)+
  ylab("Log-fold Change")+
  xlab("Delta Percent")
############fig5h##############
Positive <- FindMarkers(scRNA, ident.1 = "Tumor", min.pct = 0.1)
sig_dge.all <- subset(Positive, p_val<0.05&abs(avg_log2FC)>0.25)
sig_dge.up <- subset(Positive, p_val<0.05&avg_log2FC>0.25)
sig_dge.up <- sig_dge.up[order(sig_dge.up$avg_log2FC,decreasing = T),]
up_go_BP <- enrichGO(gene          = row.names(sig_dge.up),
                     OrgDb         = 'org.Hs.eg.db',
                     keyType       = 'SYMBOL',
                     ont           = "ALL",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05) 
up_go_BP <- data.frame(up_go_BP)
all2<- up_go_BP[order(up_go_BP$Count,decreasing = T),]
all2$Count<- all2$Count
num=c("")
select=which(all2$ID%in%paste0("GO:",num))
all2=all2[select,]
all2<- all2[order(all2$Count,decreasing = T),]
all2$Description<- factor(all2$Description,levels = rev(all2$Description))
ggplot(data = all2, aes(x = Count, y = Description, fill = pvalue)) +
  scale_fill_gradientn(colors = c("#E83D7E", "#EF84AD","#EFABC5")) +
  geom_bar(stat = "identity", width = 0.9,color = "gray",size = 0) +
  theme_classic2() +
  labs(x = "Number of Gene", y = "", title = "GO of PLAC8+CD8+TRMs in tumor")
############fig5i##############
deg<- FindMarkers(scRNA, ident.1 = "Tumor", min.pct = 0.05)
deg <- subset(deg, p_val<0.05 & abs(avg_log2FC)>0.25)
nrDEG=deg[,c(2,1)]
colnames(nrDEG)=c('log2FoldChange','pvalue')
head(nrDEG) 
gene <- bitr(rownames(nrDEG),    
             fromType = "SYMBOL",     
             toType =  "ENTREZID",   
             OrgDb = org.Hs.eg.db)   
gene$logFC <- nrDEG$log2FoldChange[match(gene$SYMBOL,rownames(nrDEG))] 
geneList=gene$logFC
names(geneList)=gene$ENTREZID 
geneList=sort(geneList,decreasing = T)
head(geneList)
go_gse <- gseGO(geneList     = geneList,
                OrgDb        = "org.Hs.eg.db",
                ont          = "BP",
                nPerm        = 1000,  
                minGSSize    = 10,
                maxGSSize    = 500,
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH",
                verbose      = FALSE)  
go <- DOSE::setReadable(go_gse, OrgDb = 'org.Hs.eg.db', keyType = 'ENTREZID')
sortgo <- go[order(go$enrichmentScore, decreasing = TRUE),]
A1="GO:"
gseaNb(object = go_gse,
       geneSetID = A1,
       newGsea = T,
       addPval = T,
       rmHt = T,
       pvalX = 0.9,
       pvalY = 0.5,
       newCurveCol = c("#377EB8","white","#E41A1C"))
############fig5g##############
gene=c("HLA−DRB5","HLA−DRA","CD74","HLA−DQA2")
A <- singlecell_gene_test(scRNA, 
                          genes.use =gene,
                          group.by = 'group', 
                          comp = c("Tumor", "Normal"))
anno_pvalue <- format(A$p_val, scientific = T,digits = 3) 
anno_sig <- A$sig
plots_violins <- VlnPlot(scRNA, 
                         cols = c("#E24840", "#1F78B4"),
                         pt.size = 0,
                         group.by = "group",
                         features = gene, 
                         ncol = 2, 
                         log = FALSE,
                         combine = FALSE)
for(i in 1:length(plots_violins)) {
  data <- plots_violins[[i]]$data
  colnames(data)[1] <- 'gene'
  plots_violins[[i]] <- plots_violins[[i]] + 
    theme_classic() + 
    theme(axis.text.x = element_text(size = 10,color="black"),
          axis.text.y = element_text(size = 10,color="black"),
          axis.title.y= element_text(size=12,color="black"),
          axis.title.x = element_blank(),
          legend.position='none')+
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+
    scale_x_discrete(labels = c("Tumor", "Normal"))+
    geom_signif(annotations = anno_sig[i],
                y_position = max(data$gene)+0.5,
                xmin = 1,
                xmax = 2,
                tip_length = 0)
}
CombinePlots(plots_violins)
############figs5##############
a="group"
pathways <- c("T cell activation", 
              "Cytokine mediated signaling pathway",
              "Lymphocyte mediated immunity",
              "Immune response regulating signaling pathwa")
plot_list <- list()
for (i in 1:length(pathways)) {
  p <- ggviolin(rds1@meta.data, x = a, y = Set, color = a, add = 'mean_sd', fill = a,
                add.params = list(color = "black")) + 
    stat_compare_means(comparisons = my_comparisons, label = "p.signif") + 
    scale_color_manual(values = my9color) + 
    scale_fill_manual(values = my9color) +
    theme(axis.text.x.bottom = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
    NoLegend() + labs(x = '')
  plot_list[[i]] <- p
}
grid.arrange(grobs = plot_list, ncol = 3)
############fig6a##############
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(SingleR)
library(cowplot)
library(tidyverse)
library(CellChat)
library(Seurat)
library(SeuratData)
library(RColorBrewer)
library(dplyr)
library(magrittr)
library(CellChat)
library(patchwork)
chatTumor <- readRDS("chatTumor.rds")
chatNormal <- readRDS("chatNormal.rds")
object.list <- list(Normal = chatNormal, Tumor = chatTumor)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
p1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
p2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
p1 + p2
############fig6b##############
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, 
                          sources.use=c(1),
                          weight.scale = T)
netVisual_diffInteraction(cellchat, 
                          sources.use=c(1),
                          weight.scale = T, measure = "weight")
p3 <- netVisual_heatmap(cellchat)
p4 <- netVisual_heatmap(cellchat, measure = "weight")
p3 + p4
############fig6c##############
netVisual_bubble(cellchat,
                 sources.use = c(1),
                 targets.use = c(2,3),
                 comparison = c(1, 2),
                 angle.x = 45)
############fig6d##############
groupSize <- as.numeric(table(cellchatT@idents))
CellChatDB = CellChatDB.human 
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use = CellChatDB
cellchatT@DB = CellChatDB.use
cellchatT <- subsetData(cellchatT)
cellchatT <- identifyOverExpressedGenes(cellchatT)
cellchatT <- identifyOverExpressedInteractions(cellchatT)
cellchatT <- computeCommunProb(object=cellchatT,raw.use = TRUE)
cellchatT <- computeCommunProbPathway(cellchatT)
cellchatT <- aggregateNet(cellchatT)
pathways.show="MHC-I"
netAnalysis_contribution(cellchatT, signaling = pathways.show)
pairLR.MHC <- extractEnrichedLR(cellchatT, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.MHC[1,] 
vertex.receiver = c(1:3) 
netVisual_individual(cellchatT, signaling = pathways.show,  layout = "circle", pairLR.use = LR.show, vertex.receiver = vertex.receiver)


pathways.show <- c("CXCL") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
pathways.show <- c("CXCL") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot
par(mfrow=c(1,1))
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
