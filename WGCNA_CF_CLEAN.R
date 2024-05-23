#######################################
# Functions taken from PloGO2 package
#######################################

# This functions are for plotting

plotErrorBarsLines <- function (v, barSizes, lines, labels = NULL, col = "blue", 
                                ylim = c(min(lines), max(lines)), ...) 
{
  barSizes[is.na(barSizes)] <- 0
  topBars <- v + 0.5 * barSizes
  bottomBars <- v - 0.5 * barSizes
  N <- length(v)
  if (is.null(labels)) 
    labels <- 1:N
  ylims <- c(min(bottomBars, ylim[1], min(lines)), max(topBars, 
                                                       ylim[2], max(lines)))
  par(pch = 19, xaxt = "n")
  plot(as.numeric(labels), v, ylim = ylims, col = col, type = "b", 
       lwd = 3, ...)
  par(xaxt = "s")
  
  for (i in 1:N) {
    lines(c(i, i), c(topBars[i], bottomBars[i]))
  }
  for (i in 1:ncol(lines)) {
    lines(as.numeric(labels), lines[, i], lwd = 0.5, lty = "dotted", 
          col = "gray")
  }
}

plotClusterProfileWGCNA <- function(cluster.data, moduleColors, group, MEs=NULL, 
                                    ylab="Abundance", 
                                    file="ClusterPatterns.png", ...) {
  
  gp = group
  noClusters <- nlevels(as.factor(moduleColors))
  
  r.temp <- aggregate(t(cluster.data), by=list(gp=gp), FUN=mean)
  ag.sample <- r.temp[,-1]
  rownames(ag.sample) <- r.temp[,1]
  ag.genes <- aggregate(t(ag.sample), by=list(Cluster=moduleColors), FUN=mean)
  ag.sd <- aggregate(t(ag.sample), by=list(Cluster=moduleColors), FUN=sd)
  ag.matrix <- as.matrix(ag.genes[,-1])
  
  if(!is.null(MEs) ) {		
    r.temp <- aggregate(MEs, by=list(gp=gp), FUN=mean)
    ag.matrix <- t(r.temp[,-1])
    colnames(ag.matrix) <- r.temp[,1]
  }
  ag.counts <- summary(as.factor(moduleColors))
  ag.bars <- as.matrix(ag.sd[,-1])
  
  fScale = max(8,noClusters)/8
  
  
  png(file, 2000, 3000*fScale, res=300)
  par(bg=gray(.95), fg=gray(0.3), mar= c(8, 6, 2, 1) + 0.1, col.main="black", col.sub="black", col.lab="black", col.axis="black")
  layout(matrix(1:(ceiling(noClusters/2)*2), ncol=2, byrow=TRUE))
  NSig <- noClusters
  cols = levels(as.factor(moduleColors) )
  for(i in 1:NSig) {
    gname <-  paste(levels(as.factor(moduleColors))[i], "(", ag.counts[i], "proteins )")
    lines <- ag.sample[, moduleColors==levels(as.factor(moduleColors))[i], drop=FALSE]
    plotErrorBarsLines(ag.matrix[i,], 2*ag.bars[i,], lines, 
                       labels=1:ncol(ag.matrix), 
                       col=cols[i],  main=gname, # bgcol="gray", split=split,
                       ylab=ylab, xlab="",
                       ylim=c(min(ag.matrix), max(ag.matrix)), ...)
    axis(1,at=1:ncol(ag.matrix), las=2, labels=colnames(ag.matrix), col="black", ...)
    abline(h=0, lty="dotted")
  }
  
  dev.off()
  
}

#-----------------------------------------------------------------#
########### Weight Gene Co-expression Network Analysis ############
#-----------------------------------------------------------------#

# Load the libraries
library(WGCNA)
library(readxl)
library(heatmap3)
enableWGCNAThreads()

# Load in the count dataframe
protein_cf <- as.data.frame(read_excel("210-CF Izq_Final Report.xlsx"))
rownames(protein_cf) <- protein_cf$UniProt
protein_cf <- protein_cf[, -c(1,2,3,4)]

# Load in experimental data
design <- read.csv("../Experimental design proteomics.csv")
design$Mouse <- colnames(protein_cf) 
datTraits <- design
rownames(datTraits) <- datTraits$Mouse
datTraits <- datTraits[,-1]

# Check for genes and samples with too many missing values
gsg <- goodSamplesGenes(protein_cf, verbose = 3)
gsg$allOK

# Cluster the samples and check for outliers
sampleTree <- hclust(dist(t(protein_cf)), method = "average")
sizeGrWindow(12,9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

# Convert traits to numeric
datTraits$Genotype <- as.factor(datTraits$Genotype)
datTraits$Diet <- as.factor(datTraits$Diet)
datTraits$Genotype <- as.numeric(datTraits$Genotype)
datTraits$Diet <- as.numeric(datTraits$Diet)
# Convert KO to 2 and WT to 1
datTraits$Genotype <- c(2,2,2,1,1,1,2,2,2,1,1,1)

# Re-cluster samples and convert traits to a color representation
sampleTree2 = hclust(dist(t(protein_cf)), method = "average")
traitColors = numbers2colors(datTraits, signed = FALSE)

# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")

# Network creation

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(t(protein_cf), powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(12, 9)
par(mfrow = c(1,2));
cex1 = 0.9;
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

# Adj
softPower = 8
adjacency = adjacency(t(protein_cf), power = softPower, type = "signed")

# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")


# Calculate module eigengenes
MEList = moduleEigengenes(t(protein_cf), colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

MEDissThres = 0.2 # CORRESPONDS TO CORRELATION >0.8 to merge modules
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(t(protein_cf), dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs

sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


# Plot module profiles with eigenproteins overlaid
WGCNAClusterID = mergedColors
Group = c("HKO", "HKO", "HKO", "HWT", "HWT", "HWT", "NKO", "NKO", "NKO", "NWT",
          "NWT", "NWT")

#Group = c("KO", "KO", "KO", "WT", "WT", "WT", "KO", "KO", "KO", "WT", "WT", "WT")

plotClusterProfileWGCNA(protein_cf, WGCNAClusterID, Group,  MEs= MEs,
                        ylab="Average log ratio", file="WGCNAClusterPattenME.png",
                        cex.main=1.8, cex.lab=1.7, cex.axis=1.5)


plotClusterProfileWGCNA(protein_cf, WGCNAClusterID, Group,  
                        ylab="Average log ratio", file="WGCNAClusterPattenAve.png",
                        cex.main=1.8, cex.lab=1.7, cex.axis=1.5)

##########################################################################

# dendrogram and heatmap for eigenproteins
png("Dendrogram eigenproteins.png", 2000,2000,res=300)					
plotEigengeneNetworks(MEs, "Eigenprotein Network", marHeatmap = c(3,4,2,2), marDendro = c(3,4,2,5),
                      plotDendrograms = TRUE, xLabelsAngle = 90,heatmapColors=blueWhiteRed(50))	
dev.off()

##########################################################################

png("Heatmap eigenproteins.png", 550,500,res=100)
heatmap3(t(MEs), #distfun = function(x) dist(x, method="euclidean"), 
         ColSideColors=rainbow(nlevels(as.factor(Group)))[as.factor(Group)],
         method = "average", 
         main="Module eigenproteins")
legend("topleft", fill=rainbow(nlevels(as.factor(Group)))[1:nlevels(as.factor(Group))],
       legend=levels(as.factor(Group)), cex=.6, xpd=TRUE, inset=-.1 )
dev.off()

##########################################################################

# Boxplot for eigenproteins

sample <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12")
ag.temp = aggregate(MEs, by=list(Group=Group), FUN=mean)
ag.eigengenes = t(ag.temp[,-1])
colnames(ag.eigengenes) = ag.temp[,1]

fScale = max(8,nlevels(as.factor(mergedColors)))/8

png("Boxplot eigenproteins.png", 2000, 3000*fScale, res=300)
par(mar= c(7, 4, 2, 1) + 0.1)
layout(matrix(1:(ceiling(nlevels(as.factor(mergedColors))/2)*2), ncol=2, byrow=TRUE))
cols = levels(as.factor(mergedColors))
for(ii in 1:ncol(MEs))	
  boxplot(MEs[,ii] ~ Group, las=2, col=cols[ii], ylab = "log ratio",
          main=paste(colnames(MEs)[ii], table(mergedColors)[ii] ), cex.main=1.7, cex.lab=1.7, cex.axis=1.5 )

dev.off()


##########################################################################

# Relate modules to clinical traits
# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs

# Define numbers of genes and samples
nGenes = ncol(t(protein_cf));
nSamples = nrow(t(protein_cf));
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(t(protein_cf), moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
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

##########################################################################

## STUDY OF EACH CLINICAL TRAIT RELATIONSHIP WITH CpGs

# Calculate module membership for each protein
modNames = substring(names(MEs), 3) # names (colors) of the modules
geneModuleMembership = as.data.frame(cor(t(protein_cf), MEs, use = "p")) # MM
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

# Calculate Intramodular Connectivity for each protein
kIn <- intramodularConnectivity(adjacency, mergedColors, scaleByMax = FALSE)

# Specify the variable containing the weight column of datTrait
gen = as.data.frame(datTraits$Genotype)
names(gen) = "Genotype"

# Calculate CpG-Trait significance for each clinical trait
geneTraitSignificance = as.data.frame(cor(t(protein_cf), gen, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(gen), sep="")
names(GSPvalue) = paste("p.GS.", names(gen), sep="")

# Get full dataframe of probes annotated
annot <- rownames(protein_cf)

# Create the starting data frame
geneInfo0 = data.frame(Uniprot = annot ,
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue,
                       kIn,
                       geneModuleMembership)

# Order modules by their significance for weight
modOrder = order(-abs(cor(MEs, gen, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.Genotype));
geneInfo = geneInfo0[geneOrder, ]

write.csv(geneInfo, "geneInfo_cortex_Diet.csv")

# Plot module membership vs. gene significance for a specific module
module = "blue"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Diet",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)


## Extract info manually by module
# Extract genes with > 0.4 significance with the trait and > 0.8 MM
hubs <- as.data.frame(rownames(geneTraitSignificance))
hubs <- cbind(hubs, geneModuleMembership$MMblue)
hubs <- cbind(hubs, geneTraitSignificance$GS.Genotype)
hubs <- cbind(hubs, kIn)
hubs <- cbind(hubs, GSPvalue)
hubs <- hubs[,c(1,2,3,8,4,5,6,7)]
colnames(hubs) <- c("Uniprot", "MM", "GS", "p-value", "kTotal", "kWithin", "kOut", "kDiff")
hubs <- hubs[abs(hubs$GS) > 0.4 & hubs$MM > 0.8,]
write.csv(hubs, "lightcyan_hubs.csv", row.names = FALSE)

## Extract hubs for each module with WGCNA's built in function
hubs <- chooseTopHubInEachModule(t(protein_cf), 
                                 colorh = dynamicColors, 
                                 omitColors = "grey", 
                                 power = 8, 
                                 type = "signed")
