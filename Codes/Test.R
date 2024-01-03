
##Kanton Data"Organoid single-cell genomic atlas uncovers human-specific features of brain development"



SOT<-readRDS("/files/scSeq_PA/analysis/timecourse_human_pseudocells_consensusGenome.rds")


pdf("tesplot.pdf")
DimPlot(object = SO, reduction = "spring", group.by = "celltype_reannot", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
dev.off()

##############################################################################################################################
#Plotting ScRNA.sce separately

NullmM<-SplitObject(seurat_integrated, split.by = "orig.ident")

ZeromM<-NullmM$S05
FivemM<-NullmM$S06

#scRNA.sce <- as.SingleCellExperiment(seurat_integrated)
ZeromM.sce <- as.SingleCellExperiment(ZeromM)

#Reference
ZeromM_SO<-SingleR(test = ZeromM.sce,ref = AR,labels = SOT$lineage)

#pred.hneurons_SO<-SingleR(test = scRNA.sce,ref = AR,labels = SOT$lineage)


cdata_Zero = cbind(colData(ZeromM.sce), 
                 data.frame(HumanNeurons.label2 = ZeromM_SO$pruned.labels))

rescounts = table(ZeromM_SO$pruned.labels, useNA="ifany")
othercells = names(rescounts)[rescounts<100]

cdata_Zero$Neurons100.label = ZeromM_SO$pruned.labels
cdata_Zero$Neurons100.label[cdata_Zero$Neurons100.label %in%othercells] = "other"

Neuronpalette = colorRampPalette(brewer.pal(19, "Paired"))(length(unique(cdata_Zero$Neurons100.label)))
names(Neuronpalette) = unique(cdata_Zero$Neurons100.label)
Neuronpalette["other"] = "gray30"

colData(ZeromM.sce) = cdata_Zero

# Adding reduced dims to the seurat converted object

umap<-Embeddings(ZeromM,"umap")[,]
reducedDim(ZeromM.sce, "UMAP")<-umap

pca<-Embeddings(ZeromM,"pca")[,]
reducedDim(ZeromM.sce, "PCA")<-pca

tsne<-Embeddings(ZeromM,"tsne")[,]
reducedDim(ZeromM.sce, "TSNE")<-tsne


colData(ZeromM.sce)$UMAP1<-data.frame(reducedDim(ZeromM.sce)[,1])
colData(ZeromM.sce)$UMAP2<-data.frame(reducedDim(ZeromM.sce)[,2])

ZeromM.sce$UMAP1<-as.numeric(ZeromM.sce$UMAP1$reducedDim.ZeromM.sce....1.)
ZeromM.sce$UMAP2<-as.numeric(ZeromM.sce$UMAP2$reducedDim.ZeromM.sce....2.)


plotScoreHeatmap(ZeromM_SO)

colsname<-unname(alphabet())[1:16]

pdf("plot_zero.pdf")
cowplot::plot_grid(
  scater::plotReducedDim(ZeromM.sce, 'UMAP', colour_by = 'Neurons100.label') 
  + scale_fill_manual(values=colsname, aesthetics="colour", name="Celltype"),
  scater::plotReducedDim(ZeromM.sce, 'UMAP', colour_by = 'seurat_clusters', text_by = 'seurat_clusters')
  + scale_fill_manual(values=colsname,aesthetics="colour",  name="Clusters"))
dev.off()
################################################################################################################


FivemM.sce <- as.SingleCellExperiment(FivemM)

#Reference
FivemM_SO<-SingleR(test = FivemM.sce,ref = AR,labels = SOT$lineage)

cdata_Five = cbind(colData(FivemM.sce), 
                   data.frame(HumanNeurons.label2 = FivemM_SO$pruned.labels))

rescounts2 = table(FivemM_SO$pruned.labels, useNA="ifany")
othercells2 = names(rescounts2)[rescounts2<100]

cdata_Five$Neurons100.label = FivemM_SO$pruned.labels
cdata_Five$Neurons100.label[cdata_Five$Neurons100.label %in%othercells2] = "other"

Neuronpalette = colorRampPalette(brewer.pal(19, "Paired"))(length(unique(cdata_Five$Neurons100.label)))
names(Neuronpalette) = unique(cdata_Five$Neurons100.label)
Neuronpalette["other"] = "gray30"

colData(FivemM.sce) = cdata_Five

# Adding reduced dims to the seurat converted object

umap<-Embeddings(FivemM,"umap")[,]
reducedDim(FivemM.sce, "UMAP")<-umap

pca<-Embeddings(FivemM,"pca")[,]
reducedDim(FivemM.sce, "PCA")<-pca

tsne<-Embeddings(FivemM,"tsne")[,]
reducedDim(FivemM.sce, "TSNE")<-tsne


colData(FivemM.sce)$UMAP1<-data.frame(reducedDim(FivemM.sce)[,1])
colData(FivemM.sce)$UMAP2<-data.frame(reducedDim(FivemM.sce)[,2])

FivemM.sce$UMAP1<-as.numeric(FivemM.sce$UMAP1$reducedDim.FivemM.sce....1.)
FivemM.sce$UMAP2<-as.numeric(FivemM.sce$UMAP2$reducedDim.FivemM.sce....2.)


plotScoreHeatmap(FivemM_SO)

colsname<-unname(alphabet())[1:16]

pdf("plot_five.pdf")
cowplot::plot_grid(
  scater::plotReducedDim(FivemM.sce, 'UMAP', colour_by = 'Neurons100.label') 
  + scale_fill_manual(values=colsname, aesthetics="colour", name="Celltype"),
  scater::plotReducedDim(FivemM.sce, 'UMAP', colour_by = 'seurat_clusters', text_by = 'seurat_clusters')
  + scale_fill_manual(values=colsname,aesthetics="colour",  name="Clusters"))
dev.off()

####################################################################################################################################
library(hdf5r)
library(rhdf5)
library(Seurat)
library(SeuratDisk)


library(rhdf5)
file = H5Fopen("HumanFetalBrainPool.h5")
file



#h5ls("HumanFetalBrainPool.h5")
CellClass<- h5read(file = "HumanFetalBrainPool.h5",name = "/shoji/CellClass")
Accession<-h5read(file = "HumanFetalBrainPool.h5",name = "/shoji/Accession")
SampleID<-h5read(file = "HumanFetalBrainPool.h5",name = "/shoji/SampleID")
Region<-h5read(file = "HumanFetalBrainPool.h5",name = "/shoji/Region")
Ngenes<-h5read(file = "HumanFetalBrainPool.h5",name = "/shoji/NGenes")
Gene<-h5read(file = "HumanFetalBrainPool.h5",name = "/shoji/Gene")
ClusterID<-h5read(file = "HumanFetalBrainPool.h5",name = "/shoji/ClusterID")
Clusters<-h5read(file = "HumanFetalBrainPool.h5",name = "/shoji/Clusters")
CellID<-h5read(file = "HumanFetalBrainPool.h5",name = "/shoji/CellID")
Tissue<-h5read(file = "HumanFetalBrainPool.h5",name = "/shoji/Tissue")
Subregion<-h5read(file = "HumanFetalBrainPool.h5",name = "/shoji/Subregion")
MeanExpression<-h5read(file = "HumanFetalBrainPool.h5",name = "/shoji/MeanExpression")
Cell_Genes<-h5read(file = "HumanFetalBrainPool.h5",name = "/shoji/Expression",index = list(NULL, 1:10000))


colnames(Cell_Genes)<-CellID[1:100]
rownames(Cell_Genes)<-Accession

test<-CreateSeuratObject(Cell_Genes)

did <- H5Dopen(file,"A")
h5writeAttribute(did, attr="shoji",name="colnames")
Attributes<-h5readAttributes("HumanFetalBrainPool.h5","shoji")



h5info=h5ls("HumanFetalBrainPool.h5")
names=h5info$name



h5writeAttribute(did, attr=names(df_A),name="colnames")

write.table(names[[3]],"Accession.txt",sep="\t")



test<-Convert(h5f$A)



A = matrix(,nrow=5,ncol=2)
h5write(A, "HumanFetalBrainPool.h5","A/H5I_DATASET")



D = h5read("HumanFetalBrainPool.h5","A/H5I_DATASET")


obj<-CreateSeuratObject(h5f)


A = matrix(1:10,nr=5,nc=2)
h5write(A, "HumanFetalBrainPool.h5","A")

fname <- tempfile(fileext = ".h5")
file <- H5File$new(fname, mode = "a")

h5_files <- list.files(pattern = "*.h5")
CB_h5<-ReadCB_h5(h5_files)


h5_read <- Read10X_h5(h5_files,use.names = F,unique.features = T)

h5_seurat <- lapply(h5_read, CreateSeuratObject)


#####################################################
