

load('combined.Rdata')


#pull in the timewise?
time_clusters<-read.table('all_timepoint_clusters.txt', sep='\t', row.names=1, header=1)

counts<-pbmc_scaled_noBatch@raw.data
meta<-pbmc_scaled_noBatch@meta.data

meta_m<-merge(meta, time_clusters, by='row.names')
row.names(meta_m)<-meta_m[,'Row.names']

genes<-row.names(pbmc_scaled_noBatch@scale.data)

unloadNamespace("Suerat")

library(URD)
library(entropy)

exp_entropy<-function(M){

logM<-log10(M+1)
d_f<-function(x){discretize(x, numBins=20, r=c(0,max(x)))}
logM_d=apply(logM,2,d_f)

entOut<-apply(logM_d,2,entropy)

}

toPlot<-merge(entOut,meta_m,by='row.names')

p<-ggplot(toPlot, aes(x=x,color=complete_clusters))+ geom_density()

p

axial <- createURD(count.data = counts, meta = meta_m, min.cells=3, min.counts=3)

axial<-findVariableGenes(axial)

#diffusion maps
axial <- calcDM(axial, knn = 100, sigma=NULL)

#pseudotime
root.cells <- cellsInCluster(axial, "init", "E720")

axial.floods <- floodPseudotime(axial, root.cells = root.cells, n=50, minimum.cells.flooded = 2, verbose=F)

axial <- floodPseudotimeProcess(axial, axial.floods, floods.name="pseudotime")

plotDimArray(axial, reduction.use = "dm", dims.to.plot = 1:8, outer.title = "Difffusion Map", label="init", plot.title="", legend=F)

plotDists(axial, "pseudotime", "init", plot.title="Pseudotime by stage")

#find the tips

#set manually? maybe 1,3,4,5,6,7,8,9,15

f1<-function(x){x %in% c(1,3,4,5,6,7,8,9,14,15)}
tips_cells<-row.names(meta_m)[sapply(meta_m['kdimension.ident'],f1)]
#selected from 

all_cluster<-unique(meta_m['complete_clusters'])

import seaborn as sns


toPlot<-merge(entOut,meta_m, by='row.names')
row.names(toPlot)<-toPlot[,'Row.names']
toPlot<-merge(axial@pseudotime,toPlot, by='row.names')
toPlot['kdimension.ident']<-as.character(toPlot['kdimension.ident'])
row.names(toPlot)<-toPlot[,'Row.names']

mean_ent<-c()
mean_pseudotime<-c()
for (cl in as.vector(t(all_cluster['complete_clusters'])))
{
	cells<-meta_m['Row.names'][meta_m['complete_clusters']==cl]

	mean_ent<-c(mean_ent,mean(entOut[cells]))
	mean_pseudotime<-c(mean_pseudotime, mean(toPlot[cells,'pseudotime']))

}

ggplot()

names(mean_ent)<-as.vector(t(all_cluster[]))
names(mean_pseudotime)<-as.vector(t(all_cluster['complete_clusters']))

#tips
axial.tips <- urdSubset(axial, cells.keep=tips_cells)

axial.tips@group.ids$`tip_clusters`<-meta_m[tips_cells,'kdimension.ident']

axial@group.ids[rownames(axial.tips@group.ids), "tip.clusters"] <- axial.tips@group.ids$`tip_clusters`

#biased random walks

axial.ptlogistic <- pseudotimeDetermineLogistic(axial, "pseudotime", optimal.cells.forward=20, max.cells.back=40, do.plot = T)

# Bias the transition matrix acording to pseudotime
axial.biased.tm <- as.matrix(pseudotimeWeightTransitionMatrix(axial, "pseudotime", logistic.params=axial.ptlogistic))

# Simulate the biased random walks from each tip
axial.walks <- simulateRandomWalksFromTips(axial, tip.group.id="tip.clusters", root.cells=root.cells, transition.matrix = axial.biased.tm, n.per.tip = 25000, root.visits = 1, max.steps = 5000)

# Process the biased random walks into visitation frequencies
axial <- processRandomWalksFromTips(axial, axial.walks, verbose = F)

# Load the cells used for each tip into the URD object
axial.tree <- loadTipCells(axial, "tip.clusters")

#tips use?
E825_tips=c('E825_1','E825_2','E825_3','E825_4','E825_5','E825_6','E825_7','E825_8','E825_9','E825_10','E825_11','E825_12','E825_13')
tips.use=c(1,3,4,5,6,7,8,9,14,15)

# Build the tree
axial.tree <- buildTree(axial.tree, tips.use=tips.use, pseudotime = "pseudotime", save.all.breakpoint.info = T, divergence.method='preference' )


complete_clusters_ordered<-c('E720_1','E720_2','E720_3','E720_4','E720_5','E720_6','E720_7','E720_8','E720_9','E750_1','E750_10','E750_2','E750_3','E750_4','E750_5','E750_6','E750_7','E750_8','E750_9','E775_1','E775_2','E775_3','E775_4','E775_5','E775_6','E825_1','E825_10','E825_11','E825_12','E825_13','E825_2','E825_3','E825_4','E825_5','E825_6','E825_7','E825_8','E825_9')


complete_clusters_colors<-c('#8A5275','#FC8D62','#FDC4AD','#6A58B3','#3c5283',
	'#8DA0CB','#66C2A5','#998DCB','#5874B3','#F781BF','#FDC4AD','#66C2A5',
	'#3F9D7F','#FC8D62','#E78AC3','#5874B3','#8A5275','#8DA0CB','#998DCB',
	'#E78AC3','#A65628','#FFD92F','#FFFF33','#FF7F00','#4DAF4A','#FFFF33','#A65628','#4DAF4A',
	'#FF7F00','#999999','#377EB8','#FFD92F','#C36910','#C1C1C1','#636363','#FF7F00','#984EA3','#E41A1C')

library(RColorBrewer)
darkcols_E825 <- c(rep("#D9D9D9",25), brewer.pal(9, "Set1"),brewer.pal(4,"Set2"))

axial.tree@meta$`kdimension.ident_d`<-factor(axial.tree@meta$`kdimension.ident`)

darkcols <- c('#E41A1C','#377EB8','#4DAF4A','#984EA3','#FF7F00',
         '#FFFF33','#A65628','#F781BF','#999999','#66C2A5',
         '#FC8D62','#8DA0CB','#000000','#A6D854','#FFD92F')

lineagecols <- c('#377EB8',
'#8DA0CB',
'#A6D854',
'#FFD92F',
'#984EA3',
'#E41A1C',
'#FC8D62',
'#FF7F00',
'#000000',
'#66C2A5',
'#FFFF33',
'#F781BF',
'#A65628',
'#999999',
'#4DAF4A')

plotTree(axial.tree, label='kdimension.ident_d', discrete.colors=darkcols, cell.alpha=0.8,cell.size=2)

#plotTree(axial.tree, label='kdimension.ident_d', discrete.colors=lineagecols, cell.alpha=0.8,cell.size=2)



tip_names<-c('Cardiac','Cranial-pharyngial','Somites',
			  'Mature endothelium','Blood', 'Lateral plate mesoderm/somites','Epithelial',
			  'Allantois','Lateral plate mesoderm','Extraembryonic')

axial.tree <- nameSegments(axial.tree, segments=tips.use, segment.names = tip_names, short.names=tip_names)

plotTree(axial.tree, label='batch', cell.alpha=0.8,cell.size=2)
plotTree(axial.tree, label='complete_clusters', discrete.colors=complete_clusters_colors, cell.alpha=0.8,cell.size=2)


save.image('complete_urd.Rdata')

axial.tree <- nameSegments(axial.tree, segments=tips.use, segment.names = tip_names, short.names=tip_names)
axial.tree <- treeForceDirectedLayout(axial.tree, start.temp=105)

plotTreeForce(axial.tree, label='kdimension.ident_d', discrete.colors=darkcols)
