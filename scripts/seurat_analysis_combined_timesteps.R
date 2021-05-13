#load some libraries that we will need
library(Seurat)
library(Matrix)
library(stringr)
library(entropy)
library(cluster)
library(RColorBrewer)
library(ggplot2)
library('monocle')

#load the data, /Users/Dan/projects/Chi_10x/age_E825/data/ is where the e8.25 data lives
pbmc.data <- Read10X("/Users/Dan/projects/Chi_10x/age_E825/data/")
pbmc_E825 <- CreateSeuratObject(pbmc.data)
pbmc.data <- Read10X("/Users/Dan/projects/Chi_10x/age_E775/data/")
pbmc_E775 <- CreateSeuratObject(pbmc.data)
pbmc.data <- Read10X("/Users/Dan/projects/Chi_10x/age_E750/data/")
pbmc_E750 <- CreateSeuratObject(pbmc.data)
pbmc.data <- Read10X("/Users/Dan/projects/Chi_10x/age_E720/data/")
pbmc_E720 <- CreateSeuratObject(pbmc.data)

pbmc.combined <- MergeSeurat(object1 = pbmc_E825, object2 = pbmc_E775, add.cell.id1 = "E825", add.cell.id2 = "E775", project = "Mesp1")
pbmc.combined <- MergeSeurat(object1 = pbmc_E750, object2 = pbmc.combined, add.cell.id1 = "E750", project = "Mesp1")
pbmc.combined <- MergeSeurat(object1 = pbmc_E720, object2 = pbmc.combined, add.cell.id1 = "E720", project = "Mesp1")
pbmc<-pbmc.combined

#rm(pbmc.combined,pbmc.data,pbmc_E825,pbmc_E775,pbmc_E750,pbmc_E720)

#find the mito genes
mito.genes <- grep("^mt-", rownames(pbmc@data), value = T)
percent.mito <- Matrix::colSums(expm1(pbmc@data[mito.genes, ])) / Matrix::colSums(expm1(pbmc@data))

#AddMetaData adds columns to object@data.info, and is a great place to stash QC stats
pbmc <- AddMetaData(pbmc, percent.mito, "percent.mito")

#if you want a violin plot of the stats, uncomment this
#VlnPlot(pbmc, c("nGene", "nUMI", "percent.mito"), nCol = 3)

#get the batch info
cell_names<-pbmc@cell.names
time_batch<-factor(str_extract(cell_names,"E[:digit:]+"))
pbmc@meta.data$batch<-time_batch

#normalize, find variable genes
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

#Read in cell cycle genes
s.genes<-readLines(con='/Users/Dan/Data/mouse_cell_cycle/s_phase.txt')
g2m.genes<-readLines(con='/Users/Dan/Data/mouse_cell_cycle/G2M_phase.txt')

#We filter out cells that have  > 5% mitochondrial percentage, nUMI < 25000
pbmc <- SubsetData(pbmc, subset.name = "percent.mito", accept.high = 0.05)
pbmc <- SubsetData(pbmc, subset.name='nUMI', accept.low = 25000)

#cell cycle scoring
pbmc<-CellCycleScoring(pbmc,s.genes=s.genes,g2m.genes=g2m.genes)

#here we split the batch correction versus not
pbmc_scaled_batch <- ScaleData(pbmc,vars.to.regress = c("percent.mito", "nUMI","batch"))
pbmc_scaled_noBatch <- ScaleData(pbmc,vars.to.regress = c("percent.mito", "nUMI"))

pbmc_scaled_batch <- RunPCA(pbmc_scaled_batch, pc.genes = pbmc@var.genes)
pbmc_scaled_noBatch<- RunPCA(object = pbmc_scaled_noBatch, pc.genes = pbmc@var.genes)

PCAPlot(object = pbmc_scaled_batch, group.by='batch')
PCAPlot(object = pbmc_scaled_noBatch, group.by='batch')

pbmc_scaled_batch <- RunTSNE(object = pbmc_scaled_batch, dims.use = 1:10, do.fast = TRUE)
pbmc_scaled_noBatch <- RunTSNE(object = pbmc_scaled_noBatch, dims.use = 1:10, do.fast = TRUE)

#attempt to justify batch correction using k means

pbmc_scaled_noBatch<-KClustDimension(pbmc_scaled_noBatch, dims.use = 1:10, reduction.use = "pca", k.use = 10, set.ident = TRUE, seed.use = 1)
pbmc_scaled_batch<-KClustDimension(pbmc_scaled_batch, dims.use = 1:10, reduction.use = "pca", k.use = 10, set.ident = TRUE, seed.use = 1)

noBatch_forMI<-mi.plugin(table(c(pbmc_scaled_noBatch@meta.data$kdimension.ident),c(pbmc_scaled_noBatch@meta.data$batch)))
batch_forMI<-mi.plugin(table(c(pbmc_scaled_batch@meta.data$kdimension.ident),c(pbmc_scaled_batch@meta.data$batch)))

#> batch_forMI
#[1] 0.3745987
#> noBatch_forMI
#[1] 0.4041218

#how many k? use silhouette

x = GetCellEmbeddings(object = pbmc_scaled_noBatch, reduction.type = "pca", dims.use = 1:10)
i=1
for (k in 8:25){

pbmc_scaled_noBatch<-KClustDimension(pbmc_scaled_noBatch, dims.use = 1:10, reduction.use = "pca", k.use = k, set.ident = TRUE, seed.use = 1)
s<-silhouette(pbmc_scaled_noBatch@meta.data$kdimension.ident,dist(x,method = 'euclidean'))
s_mean<-mean(s[,'sil_width'])
if (i==1){
s_means<-c(k,s_mean)
} else{
s_means<-rbind(s_means,c(k,s_mean))
}
i=i+1
}

plot(s_means, type='l')

#final k-means here, k=15
pbmc_scaled_noBatch<-KClustDimension(pbmc_scaled_noBatch, dims.use = 1:10, reduction.use = "pca", k.use = 15, set.ident = TRUE, seed.use = 1)
#pbmc_scaled_noBatch<-KClustDimension(pbmc_scaled_batch, dims.use = 1:10, reduction.use = "pca", k.use = 13, set.ident = TRUE, seed.use = 1)


# cluster tsne plot
#set the colors


#darkcols <- c(brewer.pal(9, "Set1"),brewer.pal(6,"Set2")) 
darkcols <- c('#E41A1C','#377EB8','#4DAF4A','#984EA3','#FF7F00',
	 '#FFFF33','#A65628','#F781BF','#999999','#66C2A5',
	 '#FC8D62','#8DA0CB','#000000','#A6D854','#FFD92F')

TSNEPlot(pbmc_scaled_noBatch,do.label=TRUE, colors.use=darkcols)
#3D

library(scatterplot3d)

tsne_1 <- pbmc_scaled_noBatch@dr$tsne@cell.embeddings[,1]

tsne_2 <- pbmc_scaled_noBatch@dr$tsne@cell.embeddings[,2]

tsne_3 <- pbmc_scaled_noBatch@dr$tsne@cell.embeddings[,3]
scatterplot3d(x = tsne_1, y = tsne_2, z = tsne_3, col=pbmc_scaled_noBatch@ident)



#get the markers
pbmc.markers <- FindAllMarkers(object = pbmc_scaled_batch, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.table(pbmc.markers,file='combined_cluster_markers.txt', quote=FALSE, sep='\t',col.names=NA)

#hierarchical clustering, x is the first ten PCA reduction
x = GetCellEmbeddings(object = pbmc_scaled_batch, reduction.type = "pca", dims.use = 1:10)
clusters<-hclust(dist(x,method = 'euclidean'))
cluster_colors<-sapply(pbmc_scaled_batch@meta.data$kdimension.ident,function(x)darkcols[x])
#save.image('combined.Rdata')


#get the top 100 most variable genes for the heatmap
top_var<-sort(apply(as.matrix(pbmc_scaled_batch@data[pbmc@var.genes,]),1,var),decreasing=TRUE)
most_variable_genes<-names(top_var)[1:100]
pdf('E825_heatmap.pdf', height=9, pointsize=9)
heatmap.2(as.matrix(pbmc_scaled_batch@data[most_variable_genes,]),Colv=as.dendrogram(clusters), ColSideColors=cluster_colors, trace='none',labCol = FALSE)
dev.off()

#You can look at any set of genes on the tsne plot
markers_from_josh_clustering<-c('Tbx4','Meox1','Lefty2','Tbx18','Trim10','Hba-x','Six2','Tbx1')
FeaturePlot(object = pbmc_scaled_batch, features.plot = markers_from_josh_clustering, cols.use = c("grey", "blue"), reduction.use = "tsne")


#here is the code for marking the gene set PC scoring
#this returns a per-gene score on the first PC
pcGenesetSignal<-function(obj, genelist)
{
genes<-readLines(genelist)
overlap<-genes[genes %in% rownames(obj@data)]
overlap_non_zero<-overlap[rowSums(as.matrix(obj@data[overlap,]))!=0]
col_non_zero<-colSums(as.matrix(obj@data[overlap,]))!=0

pc<-prcomp(obj@data[overlap_non_zero,col_non_zero], scale.=TRUE)
correct_direction<-sum(pc$x[,'PC1']>0)>(length(pc$x[,'PC1'])/2)

outscore<-rep(0,length(colnames(obj@data)))
names(outscore)<-colnames(obj@data)

if (!correct_direction){
        outscore[col_non_zero]<--1*pc$rotation[,'PC1']
        outscore[!col_non_zero]<--1*max(pc$rotation[,'PC1'])
}else{
        outscore[col_non_zero]<-pc$rotation[,'PC1']
        outscore[!col_non_zero]<-min(pc$rotation[,'PC1'])
}
return(outscore)
}

#so then you can use this to read in gene lists and score them, adding a vector in pbmc_scaled_batch@meta.data
#that you can visualize on the tsne plot
prefix<-'GenesetsFromJosh'
gene_lists<-c('YSendoderm.txt','YSmesoderm.txt','branchial_arch.txt','cardiac.txt','somites.txt')

for (gl in gene_lists){

path_file<-paste(prefix,gl, sep='/')
pbmc_scaled_batch@meta.data[[gl]]<-pcGenesetSignal(pbmc_scaled_batch,path_file)

}




library('monocle')

#beta tools
#install.packages("devtools")
#devtools::install_github("cole-trapnell-lab/monocle-release@develop")

cds<-importCDS(pbmc_scaled_noBatch)

#possible re normalize?

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
disp_table <- dispersionTable(cds)
ordering_genes <- subset(disp_table, mean_expression >= 0.1)
cds <- setOrderingFilter(cds, ordering_genes)
cds <- reduceDimension(cds)
cds <- orderCells(cds)