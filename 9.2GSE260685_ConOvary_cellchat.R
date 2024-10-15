library(Seurat)
library(dplyr)
library(Matrix)
library(patchwork)
library(ggplot2)

OVdata<-readRDS(file="../Control_Ovary_analysis_RPCA_byseurat.rds")
OVdata@meta.data$newtype<-Idents(OVdata)
OVdata<-JoinLayers(OVdata)
input<- GetAssayData(OVdata,slot = "data")
meta<-data.frame(labels=Idents(OVdata),row.names=names(Idents(OVdata)))

cellchat<-createCellChat(object=input,meta=meta,group.by="labels")
summary(cellchat)
levels(cellchat@idents)
groupSize<-as.numeric(table(cellchat@idents))
groupSize
CellChatDB<-CellChatDB.human
str(CellChatDB)
colnames(CellChatDB$interaction)
head(CellChatDB$interaction)
head(CellChatDB$cofactor)
head(CellChatDB$geneInfo)
showDatabaseCategory(CellChatDB)
unique(CellChatDB$interaction$annotation)

CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")  
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)

future::plan("multisession", workers =8) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)

cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

cellchat<-updateClusterLabels(cellchat,new.order=c("Stroma",'Endothelia','Pericyte','macrophage','T/NK','mast'))
pcount<-netVisual_heatmap(cellchat,color.heatmap="Reds",measure=c("count")) 
pweight<-netVisual_heatmap(cellchat,color.heatmap="Reds",measure="weight")
pcount+pweight 
pdf("OVdata_circle_plot_each_weight.pdf",width = 10,height = 10)
mat <- cellchat@net$weight   #mat <- cellchat@net$count
par(mfrow =c(3,3),xpd=T)
for (i in 1:nrow(mat)){
    mat2 <- matrix(0,nrow = nrow(mat),ncol = ncol(mat),dimnames = dimnames(mat))
    mat2[i,] <- mat[i,] 
    netVisual_circle(mat2,vertex.weight = groupSize, weight.scale = T,label.edge=T,arrow.size = 0.8, edge.width.max = 2.5,title.name = rownames(mat)[i])
}
dev.off()
pdf(file="OVdata_interaction_strength_heatmap.pdf",width=5,height=5)
pweight
dev.off()
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing",height=15)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming",height=15)
pdf(file="OVdata_signaling_netp_heatmap.pdf",width=18,height=20)
ht1
ht2
dev.off()
saveRDS(cellchat, file = "OVdata_cellchat.rds")




