library(Seurat)
library(dplyr)
library(Matrix)
library(patchwork)
library(ggplot2)

datafiles<-dir("/home/liusong/scRNA/endometriosis/V3data/Downloaded_Ovary_Control/matrix")
Alldata<-list()
for (i in 1:length(datafiles)){
	myfile=paste0("/home/liusong/scRNA/endometriosis/V3data/Downloaded_Ovary_Control/matrix/",datafiles[i])
	pbmc.data <-Read10X(data.dir=myfile)
	rb.genes <-rownames(pbmc.data)[grep("^RP[SL]",rownames(pbmc.data))]
	pbmc.data<-pbmc.data[!(rownames(pbmc.data) %in% rb.genes),]
	Alldata[[i]] <- CreateSeuratObject(count = pbmc.data)
}
combined = merge(Alldata[[1]],y = Alldata[-1],add.cell.ids = datafiles)
combined@meta.data$sampleID<-substr(rownames(combined@meta.data),1,6)
combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^MT-")

combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined,nfeatures=2000)
combined <- ScaleData(combined,vars.to.regress = c("nCount_RNA", "percent.mt"))
combined <- RunPCA(combined)

options(future.globals.maxSize = 8000 * 1024^2)
combined <- IntegrateLayers(object = combined,method = RPCAIntegration,orig.reduction = "pca",new.reduction = "integrated.RPCA",verbose = FALSE)
combined <- FindNeighbors(combined, reduction = "integrated.RPCA", dims = 1:30)
combined <- FindClusters(combined, resolution = 0.5, cluster.name = "rpca_clusters")
combined <- RunUMAP(combined, reduction = "integrated.RPCA", dims = 1:30, reduction.name = "umap.rpca")

saveRDS(combined,file="Control_Ovary_analysis_RPCA_byseurat.rds")
adata<-combined
pc<-DimPlot(adata)
ps<-DimPlot(adata,group.by="sampleID")
FeaturePlot(adata,features=c("LUM","PDGFRA","CD68","RGS5","VWF","TPSB2","CD3D",'GNLY'),pt.size=0.4)
adata<-RenameIdents(adata,"0"="Stroma","1"="Stroma","2"="Stroma","3"="Stroma","4"="Stroma","5"="Stroma","6"="Macrophage","7"="T/NK","8"="Stroma","9"="Endothelia","10"="Pericyte","11"="Stroma","12"="Mast","13"="Stroma","14"="Pericyte")


