library(Seurat)

load('/Users/Dan/projects/Chi_10x/cleaned_code_for_Josh/complete_urd_minimal.Rdata')
setwd('/Users/Dan/projects/Chi_10x/cardiac_tip_urd/cardiac_progenitors_and_cardiomyocytes_urd_branch/urd_v3')

pbmc_scaled_noBatch<-UpdateSeuratObject(pbmc_scaled_noBatch)


pbmc_scaled_cardiac<-subset(pbmc_scaled_noBatch, cells=axial.tree@tree$`cells.in.segment`$`1`)
pbmc_scaled_cardiac<-StashIdent(pbmc_scaled_cardiac, 'original_cluster')

pbmc_scaled_cardiac <- FindVariableFeatures(pbmc_scaled_cardiac)
pbmc_scaled_cardiac <- RunPCA(pbmc_scaled_cardiac)
pbmc_scaled_cardiac <- RunTSNE(pbmc_scaled_cardiac)
pbmc_scaled_cardiac<-FindNeighbors(pbmc_scaled_cardiac)
pbmc_scaled_cardiac<-FindClusters(pbmc_scaled_cardiac)

#pull the cardiac sub-clusters
cardiac_subclusters<-pbmc_scaled_cardiac@meta.data$seurat_clusters
cardiac_subclusters<-sapply(cardiac_subclusters, paste, "_cb", sep="")
names(cardiac_subclusters)<-rownames(pbmc_scaled_cardiac@meta.data)

#now to URD

counts<-pbmc_scaled_noBatch@assays$RNA@counts
meta<-pbmc_scaled_noBatch@meta.data

for_tips<-pbmc_scaled_noBatch@meta.data$kdimension.ident
names(for_tips)<-rownames(pbmc_scaled_noBatch@meta.data)

for_tips[names(cardiac_subclusters)]<-cardiac_subclusters

unloadNamespace("Suerat")

#previously use 1.0.2
#update to
library(URD)

axial <- createURD(count.data = counts, meta = meta, min.cells=3, min.counts=3)
axial<-findVariableGenes(axial)
axial <- calcDM(axial, knn = 100, sigma=NULL)

root.cells <-cellsInCluster(axial, "init", "E720")

axial.floods <- floodPseudotime(axial, root.cells = root.cells, n=50, minimum.cells.flooded = 2, verbose=F)

axial <- floodPseudotimeProcess(axial, axial.floods, floods.name="pseudotime")

#URD_cb_v2: Tips 3, 1, 2 - THF (16); 0, 5 - SHF(18); 6 - FHF(20)

for_tips_num<-sub("0_cb","18",for_tips)
for_tips_num<-sub("1_cb","16",for_tips_num) #
for_tips_num<-sub("2_cb","16",for_tips_num) #
for_tips_num<-sub("3_cb","16",for_tips_num) #THF
for_tips_num<-sub("4_cb","19",for_tips_num) #
for_tips_num<-sub("5_cb","18",for_tips_num) #SHF
for_tips_num<-sub("6_cb","20",for_tips_num) #FHF

for_tips_num<-as.numeric(for_tips_num)
names(for_tips_num)<-names(for_tips)

f1<-function(x){x %in% c(16,20,18,3,4,5,6,7,8,9,14,15)}
tips_cells<-names(for_tips_num)[sapply(for_tips_num,f1)]

axial.tips <- urdSubset(axial, cells.keep=tips_cells)

axial.tips@group.ids$`tip_clusters`<-for_tips_num[tips_cells]
axial@group.ids[rownames(axial.tips@group.ids), "tip.clusters"]<-axial.tips@group.ids$`tip_clusters`

#biased random walks

axial.ptlogistic <- pseudotimeDetermineLogistic(axial, "pseudotime", optimal.cells.forward=20, max.cells.back=40)

# Bias the transition matrix acording to pseudotime
axial.biased.tm <- as.matrix(pseudotimeWeightTransitionMatrix(axial, "pseudotime", logistic.params=axial.ptlogistic))

# Simulate the biased random f1<-function(x){x %in% c(16,20,18,3,4,5,6,7,8,9,14,15)}
tips_cells<-names(for_tips_num)[sapply(for_tips_num,f1)]

axial.tips <- urdSubset(axial, cells.keep=tips_cells)

axial.tips@group.ids$`tip_clusters`<-for_tips_num[tips_cells]
axial@group.ids[rownames(axial.tips@group.ids), "tip.clusters"]<-axial.tips@group.ids$`tip_clusters`

#biased random walks

axial.ptlogistic <- pseudotimeDetermineLogistic(axial, "pseudotime", optimal.cells.forward=20, max.cells.back=40)

# Bias the transition matrix acording to pseudotime
axial.biased.tm <- as.matrix(pseudotimeWeightTransitionMatrix(axial, "pseudotime", logistic.params=axial.ptlogistic))

# Simulate the biased random walks from each tip
axial.walks <- simulateRandomWalksFromTips(axial, tip.group.id="tip.clusters", root.cells=root.cells, transition.matrix = axial.biased.tm, n.per.tip = 25000, root.visits = 1, max.steps = 5000)

#axial.walks <- simulateRandomWalksFromTips(axial, tip.group.id="tip.clusters", root.cells=root.cells, transition.matrix = axial.biased.tm)

# Process the biased random walks into visitation frequencies
axial <- processRandomWalksFromTips(axial, axial.walks, verbose = F)

# Load the cells used for each tip into the URD object
axial.tree <- loadTipCells(axial, "tip.clusters")

tips_use<-c(16,20,18,3,4,5,6,7,8,9,14,15)

axial.tree <- buildTree(axial.tree, tips.use=tips_use, pseudotime = "pseudotime", save.all.breakpoint.info = T, divergence.method='preference' )

tip_names<-c('Cardiac THF', 'Cardiac FHF', 'Cardiac SHF','Cranial-pharyngial','Somites',
                          'Mature endothelium','Blood', 'Lateral plate mesoderm/somites','Epithelial',
                          'Allantois','Lateral plate mesoderm','Extraembryonic')

axial.tree@meta$`kdimension.ident_d`<-factor(axial.tree@meta$`kdimension.ident`)

darkcols <- c('#E41A1C','#377EB8','#4DAF4A','#984EA3','#FF7F00',
         '#FFFF33','#A65628','#F781BF','#999999','#66C2A5',
         '#FC8D62','#8DA0CB','#000000','#A6D854','#FFD92F')
axial.tree <- nameSegments(axial.tree, segments=tips_use, segment.names = tip_names, short.names=tip_names)

png(file='binary_urd.png')
print(plotTree(axial.tree, label='kdimension.ident_d', discrete.colors=darkcols, cell.alpha=0.8,cell.size=2))
dev.off()

stage_cols<-c(brewer.pal(4,"Dark2"))
png(file='binary_urd_stage.png')
print(plotTree(axial.tree, label='batch', discrete.colors=stage_cols, cell.alpha=0.8,cell.size=2))
dev.off()

#save.image('urd_cb_v3.Rdata')

axial.tree <- nameSegments(axial.tree, segments=tips_use, segment.names = tip_names, short.names=tip_names)
axial.tree <- treeForceDirectedLayout(axial.tree, start.temp=105,cut.unconnected.segments=7)

#plotTreeForce(axial.tree, label='kdimension.ident_d', discrete.colors=darkcols, label.tips=FALSE, symmetric.color.scale =FALSE)
#plotTreeForce(axial.tree, label='batch', discrete.colors=darkcols, label.tips=FALSE, symmetric.color.scale =FALSE)


