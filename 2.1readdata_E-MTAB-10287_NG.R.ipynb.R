library(Seurat)
library(SeuratDisk)

options(Seurat.object.assay.version = "v3")
datafiles<-dir("/home/liusong/scRNA/endometriosis/V3data/Downloaded_NG_EUdata/matrix_data")
Alldata<-list()
for (i in 1:length(datafiles)){
        mymtx=paste0("/home/liusong/scRNA/endometriosis/V3data/Downloaded_NG_EUdata/matrix_data/",datafiles[i],"/matrix.mtx.gz")
        mybarcodes=paste0("/home/liusong/scRNA/endometriosis/V3data/Downloaded_NG_EUdata/matrix_data/",datafiles[i],"/barcodes.tsv.gz")
        myfeatures=paste0("/home/liusong/scRNA/endometriosis/V3data/Downloaded_NG_EUdata/matrix_data/",datafiles[i],"/features.tsv.gz")
        pbmc.data <-ReadMtx(mtx=mymtx,cells=mybarcodes,features=myfeatures,cell.column=1,feature.column=1,skip.cell=1)
        rb.genes <-rownames(pbmc.data)[grep("^RP[SL]",rownames(pbmc.data))]
        pbmc.data<-pbmc.data[!(rownames(pbmc.data) %in% rb.genes),]
        pbmc <- CreateSeuratObject(count = pbmc.data)
	metainfo<-read.csv(file=mybarcodes,sep="\t")
	mymeta<-metainfo[,c("DonorID","Binary.Stage","broad_celltypes","general_celltypes","fine_celltypes")]
	pbmc<- AddMetaData(pbmc,mymeta,col.name=c("DonorID","Phase","CellTypeA","CellTypeB","CellTypeC"))
	Alldata[[i]]<-pbmc
}
combined = merge(Alldata[[1]],y = Alldata[-1])
library(SeuratDisk)
SaveH5Seurat(combined,filename="NG_EU_rawdata.h5seurat")
Convert("NG_EU_rawdata.h5seurat",dest="h5ad")
saveRDS(combined,file="NG_EU_rawdata_noRPL.rds")

