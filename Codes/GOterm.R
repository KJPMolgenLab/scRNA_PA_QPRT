library(org.Hs.eg.db)
library(Seurat)
library(topGO)
library(gprofiler2)
library("biomaRt")

seurat_integrated<-readRDS("seurat_integrated.rds") # Read the scRNA data
scRNA.sce <- as.SingleCellExperiment(seurat_integrated)

Idents(seurat_integrated) <- seurat_integrated$seurat_clusters

#Test
cluster15<- subset(seurat_integrated,idents="15")

expr <- cluster15[["RNA"]]$counts
# Select genes that are expressed > 0 in at least 75% of cells (somewhat arbitrary definition)
n.gt.0 <- apply(expr, 1, function(x)length(which(x > 0)))
expressed.genes <- rownames(expr)[which(n.gt.0/ncol(expr) >= 0.5)]
all.genes <- rownames(expr)

# define geneList as 1 if gene is in expressed.genes, 0 otherwise
geneList <- ifelse(all.genes %in% expressed.genes, 1, 0)
names(geneList) <- all.genes


markers<-read.table("../../Cluster_markers (2).txt",sep="\t")



genelist = list()

for(i in 1:15){
  genelist[[i]] = markers$gene[markers$cluster == i ]
}


# Gene universe file
bg_genes=rownames(seurat_integrated)
length(bg_genes)

# Read in genes of interest
candidate_list=genelist[[15]]
length(candidate_list)

# create GO db for genes to be used using biomaRt - please note that this takes a while
db= useMart('ENSEMBL_MART_ENSEMBL',dataset='hsapiens_gene_ensembl', host="www.ensembl.org")
go_ids= getBM(attributes=c('go_id', 'external_gene_name', 'namespace_1003'), filters='external_gene_name', values=bg_genes, mart=db)

# build the gene 2 GO annotation list (needed to create topGO object)
gene_2_GO=unstack(go_ids[,c(1,2)])

# remove any candidate genes without GO annotation
keep = candidate_list %in% go_ids[,2]
keep =which(keep==TRUE)
candidate_list=candidate_list[keep]
# make named factor showing which genes are of interest
geneList=factor(as.integer(bg_genes %in% candidate_list))
names(geneList)= bg_genes

GOdata=new('topGOdata', ontology='BP', allGenes = geneList, annot = annFUN.gene2GO, gene2GO = gene_2_GO)
resultFisher <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
GenTable(GOdata, Fisher = resultFisher, topNodes = 20, numChar = 60)

write.table(GenTable(GOdata, Fisher = resultFisher, topNodes = 20, numChar = 60),"GO_terms_clust16_BP.txt",sep="\t") # since starting from 0



###################################################################################


genebackground = rownames(seurat_integrated)
gostres = gost(query = genelist, evcodes = T, custom_bg = genebackground)

names(gostres)

tablegostres = gostres$result[grepl("GO:", gostres$result$source), c("query", "source", "term_name","p_value", "intersection")]










              