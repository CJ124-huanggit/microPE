
#package
library(WGCNA)
library(tidyverse)

##data
ps = base::readRDS("./data/dataNEW/ps_ITS.rds")

##filter
ps0<-ps %>% scale_micro() %>%
  filter_taxa(function(x) sum(x ) > 0.01 , TRUE)
ps0
#
source("C:\\SANREN6\\micro\\total_amplicon.R")
#environment

res = theme_my()
mytheme1 = res[[1]];mytheme2 = res[[2]]; colset1 = res[[3]];colset2 = res[[4]];colset3 = res[[5]];colset4 = res[[6]]
result<- dir.amp(ps0 = ps0)
id = result[[2]];id

env = read.csv("./data/dataNEW/ITS_env.csv")
head(env)
envRDA = env
head(env)
row.names(envRDA) = env$ID
envRDA$ID = NULL
head(envRDA)

traitData = envRDA
allTraits = envRDA
datTraits = envRDA


if (is.na(match("Fungi",id))) {
  res1path <- "./result_and_plot/Micro_and_other_index_16s_230410"
  
} else if(is.na(match("Bacteria",id))) {
  res1path <- "./result_and_plot/Micro_and_other_index_ITS"
}
dir.create(res1path)Top = 

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
psnet  = ps0 %>%
  ggClusterNet::scale_micro(method = "TMM") %>%
  ggClusterNet::filter_OTU_ps(Top)
datExpr0 = ps0 %>%
  ggClusterNet::scale_micro(method = "TMM") %>%
  ggClusterNet::filter_OTU_ps(Top) %>%
  ggClusterNet::vegan_otu() %>% 
  as.data.frame()
datExpr0[1:5,1:5]

sampleTree = hclust(dist(datExpr0), method = "average")
pdf(file = paste(WGCNApath,"/","1_sampleClustering.pdf",sep = ""), width = 12, height = 9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
### Plot a line to show the cut
abline(h = 10000, col = "red")
dev.off()


sampleTree2 = hclust(dist(datExpr0), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(allTraits, signed = FALSE)
# Plot the sample dendrogram and the colors underneath.
#sizeGrWindow(12,12)


pdf(file=paste(WGCNApath,"/","2_Sample dendrogram and trait heatmap.pdf",sep = ""),width=12,height=12)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(allTraits),
                    main = "Sample dendrogram and trait heatmap")
dev.off()




#############################network constr########################################
# See note above.
enableWGCNAThreads()
# Choose a set of soft-thresholding powers
powers = c(1:30)
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)



pdf(file=paste(WGCNApath,"/","3_Scale independence.pdf",sep = ""),width=9,height=5)
par(mfrow = c(1,2))
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()



######chose the softPower
# softPower =sft$powerEstimate
softPower = 18
#
adjacency = adjacency(datExpr0, power = softPower)

##### Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency); 
#TOM
dissTOM = 1-TOM
# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)

#sizeGrWindow(12,9)
pdf(file=paste(WGCNApath,"/","4_Micro clustering on TOM-based dissimilarity.pdf",sep = ""),width=12,height=9)
plot(geneTree, xlab="", sub="", main = "Micro clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()

#
minModuleSize = 30
# Module identification using dynamic tree cut
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
#
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
#sizeGrWindow(8,6)

pdf(file=paste(WGCNApath,"/","5_Dynamic Tree Cut.pdf",sep = ""),width=8,height=6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

#
# Calculate eigengenes
datExpr0[1:5,1:5]

MEList = moduleEigengenes(datExpr0, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")
# Plot the result
#sizeGrWindow(7, 6)
pdf(file=paste(WGCNApath,"/","6_Clustering of module eigengenes.pdf",sep = ""),width=7,height=6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.05######
# Plot the cut line into the dendrogram，
abline(h=MEDissThres, col = "red")
dev.off()


# Call an automatic merging function
merge = mergeCloseModules(datExpr0, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs

write.csv(mergedMEs,file = "its_WGCNA_DXAL.csv")

#sizeGrWindow(12, 9)
#
pdf(file=paste(WGCNApath,"/","7_merged dynamic.pdf",sep = ""), width = 9, height = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# Rename to moduleColors
moduleColors = mergedColors
table(moduleColors)
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
table(moduleLabels)
MEs = mergedMEs

#
# Save module colors and labels for use in subsequent parts
save(MEs, TOM, dissTOM,  moduleLabels, moduleColors, geneTree, sft, file = paste(WGCNApath,"/","networkConstruction-stepByStep.RData",sep = ""))

load(paste(WGCNApath,"/","networkConstruction-stepByStep.RData",sep = ""))


#
nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

#sizeGrWindow(10,6)
pdf(file=paste(WGCNApath,"/","8_Module-trait relationships.pdf",sep = ""),width=15,height=6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")

dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()


#--12
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^7
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA


if (F) {
  pdf(file= paste(WGCNApath,"/","12_Network heatmap plot_all gene.pdf",sep = ""),width=9, height=9)
  TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")
  dev.off()
}


#
nSelect = 400
# For reproducibility, we set the random seed
set.seed(10)
select = sample(nGenes, size = nSelect)
selectTOM = dissTOM[select, select]
# There's no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select]

# Open a graphical window
#sizeGrWindow(9,9)
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing
# the color palette; setting the diagonal to NA also improves the clarity of the plot
plotDiss = selectTOM^7
diag(plotDiss) = NA


pdf(file=paste(WGCNApath,"/","13_Network heatmap plot_selected genes.pdf",sep = ""),width=9, height=9)
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")
dev.off()


#


gru = data.frame(ID =  names(datExpr0), group = moduleColors )
head(gru)

cor = TOM
cor[1:10,1:10]
row.names(cor) = names(datExpr0)
colnames(cor) = names(datExpr0)


library(ggClusterNet)
node = PolygonMaptreeG(cor = cor,seed = 12)
head(node)

otu_table = as.data.frame(t(ggClusterNet::vegan_otu(psnet)))
tax_table = as.data.frame(ggClusterNet::vegan_tax(psnet))

nodes = nodeadd(plotcord =node,otu_table = otu_table,tax_table = tax_table)
head(nodes)
#-----edge
edge = edgeBuild(cor = cor,node = node)
colnames(edge)[8] = "cor"
head(edge)



nodeG = nodes %>% 
  dplyr::left_join(gru,by = c( "elements" = "ID"))
head(nodeG)

write.csv(nodeG,file = "ITS_nodeG_DXAL.csv")

p1 <- ggplot() +
  # geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = cor),
  #              data = edge, size = 0.5,alpha = 0.6) +
  # geom_curve(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = cor),
  #            data = edge, size = 0.5,alpha = 0.3,curvature = -0.2) +
  geom_point(aes(X1, X2,size = mean),pch = 21, data = nodeG, fill = nodeG$group,color = "gray40") +
  # geom_text_repel(aes(X1, X2,label = elements),pch = 21, data = nodes) +
  # geom_text(aes(X1, X2,label = elements),pch = 21, data = nodes) +
  scale_colour_manual(values = c("#377EB8","#E41A1C")) +
  scale_size(range = c(5, 15)) +
  scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
  theme(panel.background = element_blank()) +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank()
        
        ) +
  theme(legend.background = element_rect(colour = NA)) +
  theme(panel.background = element_rect(fill = "white",  colour = NA)) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())

filename = paste(WGCNApath,"/","14_Network_all_1.pdf",sep = "")
ggsave(filename,p1,width = 42,height = 42,limitsize = FALSE)


