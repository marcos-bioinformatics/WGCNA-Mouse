#-----------------------------------------------------------------#
# Weighted protein co-expression network analysis (WGCNA) of the
# frontal cortex proteome of control and GLUT1-deficient mice
#
# Author: Marcos Elizalde-Horcada
# License: GPL-2 (see LICENSE)
#
# Code for fig. S3 of: Ardanaz CG, et al. Astrocytic GLUT1 reduction
# paradoxically improves central and peripheral glucose homeostasis.
# Sci Adv. 2024;10(42):eadp1115. https://doi.org/10.1126/sciadv.adp1115
#
# The network construction follows the WGCNA tutorials by Peter Langfelder
# and Steve Horvath (Langfelder & Horvath, BMC Bioinformatics 2008;9:559).
# plotErrorBarsLines(), plotClusterProfileWGCNA() and the eigenprotein
# plots come from the WGCNA_proteomics.R script of the PloGO2 package
# (Wu et al., J Proteome Res 2020;19:2898-2906; license GPL-2).
#
# Revision notes (October 2026)
# Parameters and results are unchanged except where noted.
#  - Eigenprotein plots (cluster profiles, dendrogram, heatmap, boxplots)
#    now use the eigenproteins of the merged modules, ordered by colour,
#    as in the PloGO2 script. They used the pre-merge eigenproteins with
#    the merged colours, which mislabels the plots if any module merges.
#  - pickSoftThreshold() now evaluates a signed network, the type built
#    with adjacency(). The power used (8) is unchanged, so the modules are
#    the same; only the scale-free fit plot changes.
#  - goodSamplesGenes() now receives samples in rows.
#  - Output names and labels refer to genotype, the trait analysed
#    (geneInfo_cortex_Genotype.csv; it was named "Diet").
#  - Hub proteins are taken only from the proteins assigned to the module
#    and saved as <module>_hubs.csv (it was lightcyan_hubs.csv while the
#    code used the blue module). This can change the hub list.
#  - chooseTopHubInEachModule() uses the merged module colours.
#  - Sample groups and genotype are read from the column names of the
#    report instead of hard-coded vectors, and the script stops if the
#    design file does not list the mice in the same order.
#  - Axis labels say "log abundance" (the input values are log-transformed
#    abundances, not ratios) and "Eigenprotein value" for the eigenproteins.
#  - Removed duplicated module membership columns from geneInfo, unused
#    variables and comments carried over from earlier projects; plot titles
#    refer to proteins instead of genes.
#  - License changed from GPL-3 to GPL-2, the license of PloGO2.
#-----------------------------------------------------------------#

#######################################
# Functions taken from PloGO2 package
# (inst/script/WGCNA_proteomics.R, J. Wu, 2020)
#######################################

# Plotting functions

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
######### Weighted Protein Co-expression Network Analysis #########
#-----------------------------------------------------------------#

# Load the libraries
library(WGCNA)
library(readxl)
library(heatmap3)
enableWGCNAThreads()

# Load the protein abundance table (SWATH-MS report, first sheet): proteins
# in rows, mice in columns after the four annotation columns (UniProt,
# Protein name, Gene, Peptides). Abundances are already log-transformed and
# centred on the median of each sample.
protein_cf <- as.data.frame(read_excel("210-CF Izq_Final Report.xlsx"))
rownames(protein_cf) <- protein_cf$UniProt
protein_cf <- protein_cf[, -c(1,2,3,4)]

# Load the experimental design
# (one row per mouse, in the same order as the columns of protein_cf)
design <- read.csv("../Experimental design proteomics.csv")
design$Mouse <- colnames(protein_cf)
datTraits <- design
rownames(datTraits) <- datTraits$Mouse
datTraits <- datTraits[,-1]

# Sample groups from the column names, e.g. "HFD_KO_CF Izq-118": diet
# (HFD = high-fat diet, NCD = normal chow diet) and genotype (KO, WT)
diet = sub("_.*", "", colnames(protein_cf))
genotype = sub("^[^_]*_([^_]*)_.*", "\\1", colnames(protein_cf))
Group = paste0(substr(diet, 1, 1), genotype) # HKO, HWT, NKO, NWT
table(Group)
stopifnot(all(Group %in% c("HKO", "HWT", "NKO", "NWT")))

# Alternative grouping by genotype only, for the plots:
#Group = genotype

# Check that the design file lists the mice in the same order as the columns
stopifnot(all(rowSums(table(design$Diet, diet) > 0) == 1),
          all(rowSums(table(design$Genotype, genotype) > 0) == 1))

# Check for proteins and samples with too many missing values
# (WGCNA functions expect samples in rows and proteins in columns)
gsg <- goodSamplesGenes(t(protein_cf), verbose = 3)
gsg$allOK

# Cluster the samples and check for outliers
sampleTree <- hclust(dist(t(protein_cf)), method = "average")
sizeGrWindow(12,9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

# Convert traits to numeric
# Genotype: WT is 1, KO is 2
datTraits$Genotype <- ifelse(genotype == "KO", 2, 1)
# Diet: levels are numbered in alphabetical order, check the coding
datTraits$Diet <- as.numeric(as.factor(datTraits$Diet))
table(design$Diet, datTraits$Diet)

# Convert traits to a color representation
traitColors = numbers2colors(datTraits, signed = FALSE)

# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")

# Network creation

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
# (networkType must match the signed adjacency built below)
sft = pickSoftThreshold(t(protein_cf), powerVector = powers,
                        networkType = "signed", verbose = 5)
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
# R^2 cut-off of 0.90
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# Build the signed adjacency matrix with the chosen power
softPower = 8
adjacency = adjacency(t(protein_cf), power = softPower, type = "signed")

# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Protein clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)

# Minimum module size: 30 proteins
minModuleSize = 30;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)

# Convert numeric labels into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Protein dendrogram and module colors")


# Calculate the eigenproteins of the modules before merging
MEList = moduleEigengenes(t(protein_cf), colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigenproteins
MEDiss = 1-cor(MEs);
# Cluster module eigenproteins
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigenproteins",
     xlab = "", sub = "")

MEDissThres = 0.2 # Merge modules whose eigenproteins correlate > 0.8
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(t(protein_cf), dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigenproteins of the new merged modules:
mergedMEs = merge$newMEs

sizeGrWindow(12, 9)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# Use the eigenproteins of the merged modules, ordered by module colour so
# that they match levels(as.factor(mergedColors)) in the plots below
MEs = mergedMEs
MEs = MEs[, order(colnames(MEs))]

# Plot module profiles with eigenproteins overlaid
WGCNAClusterID = mergedColors

plotClusterProfileWGCNA(protein_cf, WGCNAClusterID, Group,  MEs= MEs,
                        ylab="Average log abundance", file="WGCNAClusterPattenME.png",
                        cex.main=1.8, cex.lab=1.7, cex.axis=1.5)


plotClusterProfileWGCNA(protein_cf, WGCNAClusterID, Group,
                        ylab="Average log abundance", file="WGCNAClusterPattenAve.png",
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

ag.temp = aggregate(MEs, by=list(Group=Group), FUN=mean)
ag.eigengenes = t(ag.temp[,-1])
colnames(ag.eigengenes) = ag.temp[,1]

fScale = max(8,nlevels(as.factor(mergedColors)))/8

png("Boxplot eigenproteins.png", 2000, 3000*fScale, res=300)
par(mar= c(7, 4, 2, 1) + 0.1)
layout(matrix(1:(ceiling(nlevels(as.factor(mergedColors))/2)*2), ncol=2, byrow=TRUE))
cols = levels(as.factor(mergedColors))
for(ii in 1:ncol(MEs))
  boxplot(MEs[,ii] ~ Group, las=2, col=cols[ii], ylab = "Eigenprotein value",
          main=paste(colnames(MEs)[ii], table(mergedColors)[ii] ), cex.main=1.7, cex.lab=1.7, cex.axis=1.5 )

dev.off()


##########################################################################

# Relate modules to traits (genotype and diet)
# Rename to moduleColors
moduleColors = mergedColors

# Number of samples
nSamples = nrow(t(protein_cf));
# Recalculate the eigenproteins of the merged modules, ordered by similarity
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

## RELATIONSHIP OF EACH PROTEIN WITH GENOTYPE

# Calculate module membership for each protein
modNames = substring(names(MEs), 3) # names (colors) of the modules
geneModuleMembership = as.data.frame(cor(t(protein_cf), MEs, use = "p")) # MM
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

# Calculate Intramodular Connectivity for each protein
kIn <- intramodularConnectivity(adjacency, moduleColors, scaleByMax = FALSE)

# Trait of interest: genotype
gen = as.data.frame(datTraits$Genotype)
names(gen) = "Genotype"

# Calculate protein significance (GS) for genotype
geneTraitSignificance = as.data.frame(cor(t(protein_cf), gen, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(gen), sep="")
names(GSPvalue) = paste("p.GS.", names(gen), sep="")

# Protein IDs (UniProt accessions)
annot <- rownames(protein_cf)

# Create the starting data frame
geneInfo0 = data.frame(Uniprot = annot ,
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue,
                       kIn)

# Order modules by their significance for genotype
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
# Order the proteins first by module color, then by |GS|
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.Genotype));
geneInfo = geneInfo0[geneOrder, ]

write.csv(geneInfo, "geneInfo_cortex_Genotype.csv")

# Plot module membership vs. protein significance for the module of interest
module = "blue"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Protein significance for Genotype",
                   main = paste("Module membership vs. protein significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)


## Extract hub proteins of the module
# Proteins of the module with |GS| > 0.4 for genotype and MM > 0.8
hubs <- data.frame(Uniprot = rownames(geneTraitSignificance),
                   MM = geneModuleMembership[, column],
                   GS = geneTraitSignificance$GS.Genotype,
                   "p-value" = GSPvalue$p.GS.Genotype,
                   kIn, check.names = FALSE)
hubs <- hubs[moduleGenes & abs(hubs$GS) > 0.4 & hubs$MM > 0.8,]
write.csv(hubs, paste0(module, "_hubs.csv"), row.names = FALSE)

## Top hub protein of each module with WGCNA's built in function
topHubs <- chooseTopHubInEachModule(t(protein_cf),
                                    colorh = moduleColors,
                                    omitColors = "grey",
                                    power = softPower,
                                    type = "signed")
topHubs
