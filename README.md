# WGCNA-Mouse

Weighted protein co-expression network analysis (WGCNA) of the frontal cortex proteome of mice with astrocyte-specific GLUT1 reduction (GLUT1<sup>ΔGFAP</sup>) and their controls.

This repository contains the R code for the network analysis in fig. S3 of:

> Ardanaz CG, de la Cruz A, Minhas PS, Hernández-Martín N, Pozo MÁ, Valdecantos MP, Valverde ÁM, Villa-Valverde P, **Elizalde-Horcada M**, Puerta E, Ramírez MJ, Ortega JE, Urbiola A, Ederra C, Ariz M, Ortiz-de-Solórzano C, Fernández-Irigoyen J, Santamaría E, Karsenty G, Brüning JC, Solas M. **Astrocytic GLUT1 reduction paradoxically improves central and peripheral glucose homeostasis.** *Science Advances* 10(42): eadp1115 (2024). https://doi.org/10.1126/sciadv.adp1115

I wrote and ran this analysis (in the paper: data curation, formal analysis, methodology, software and validation).

## Study design

| | |
|---|---|
| Tissue | Frontal cortex |
| Proteomics | SWATH-MS (Sciex TripleTOF 5600+), 2,329 proteins quantified |
| Samples | 12 mice: 6 controls and 6 GLUT1<sup>ΔGFAP</sup>; in each genotype, 3 on normal chow diet and 3 on high-fat diet |
| Traits | Genotype (WT = 1, GLUT1<sup>ΔGFAP</sup> = 2) and diet |

## Analysis

`WGCNA_CF.R` runs the whole workflow:

| Step | WGCNA functions | Settings |
|---|---|---|
| Sample QC | `goodSamplesGenes`, `hclust` | Average-linkage clustering of samples to detect outliers |
| Soft-thresholding power | `pickSoftThreshold` | Powers 1–10 and 12–20, signed network |
| Network | `adjacency`, `TOMsimilarity` | Signed network, β = 8, topological overlap (TOM) |
| Module detection | `hclust`, `cutreeDynamic` | Average linkage on 1 − TOM, `deepSplit = 2`, `minClusterSize = 30` |
| Module merging | `mergeCloseModules` | `cutHeight = 0.2` (eigenprotein correlation > 0.8) |
| Module–trait relationships | `moduleEigengenes`, `cor`, `corPvalueStudent` | Pearson correlation of each module eigenprotein with genotype and diet |
| Protein-level measures | `intramodularConnectivity` | Module membership (MM), protein significance for genotype (GS), intramodular connectivity |
| Hub proteins | | Proteins of the module with \|GS\| > 0.4 and MM > 0.8 |

## Main result

The network had 27 modules. Only the **blue module** (179 proteins) correlated with genotype at p < 0.01 (r = −0.78, p = 0.003), with a lower eigenprotein in GLUT1<sup>ΔGFAP</sup> mice. Gene Ontology enrichment of the blue module (ShinyGO 0.76, run outside this script, with the 2,329 detected proteins as background) pointed to energy metabolism: tricarboxylic acid cycle, aerobic respiration and oxidative phosphorylation.

## Data

The input data are not included in this repository:

- Raw mass spectrometry data: ProteomeXchange [PXD035376](https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD035376) / jPOST [JPST001793](https://repository.jpostdb.org/entry/JPST001793).
- Source data for the paper's figures: Zenodo, https://doi.org/10.5281/zenodo.12740848.

The script expects two files:

- `210-CF Izq_Final Report.xlsx`: SWATH quantification report. The script reads the first sheet, with one row per protein, four annotation columns (`UniProt`, `Protein name`, `Gene`, `Peptides`) and then one column per mouse with log-transformed abundances centred on the median of each sample. The mouse columns are named `<diet>_<genotype>_<sample>`, for example `HFD_KO_CF Izq-118`; the script takes the groups from these names.
- `../Experimental design proteomics.csv`: one row per mouse, in the same order as the columns of the report, with the mouse ID in the first column and the `Genotype` and `Diet` columns. The script stops if this order does not match the column names.

## Outputs

| File | Content |
|---|---|
| `WGCNAClusterPattenME.png`, `WGCNAClusterPattenAve.png` | Abundance profile of each module per group, with the eigenprotein or the mean overlaid |
| `Dendrogram eigenproteins.png`, `Heatmap eigenproteins.png`, `Boxplot eigenproteins.png` | Relationships between module eigenproteins and their values per group |
| `geneInfo_cortex_Genotype.csv` | Module, GS, connectivity and MM of every protein |
| `blue_hubs.csv` | Hub proteins of the blue module |

The script also draws the sample dendrogram, scale-free fit, module dendrogram, module–trait heatmap and the MM vs GS plot of the blue module in the graphics device.

## Requirements

R (≥ 4.0) with the packages WGCNA, readxl and heatmap3:

```r
install.packages("BiocManager")
BiocManager::install(c("WGCNA", "readxl", "heatmap3"))
```

## Usage

Put the input files where the script expects them, set the working directory to the folder of the report and run:

```r
source("WGCNA_CF.R")
```

## Credits

- The network construction and module analysis follow the [WGCNA tutorials](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/) by Peter Langfelder and Steve Horvath.
- `plotErrorBarsLines()`, `plotClusterProfileWGCNA()` and the eigenprotein plots come from the `WGCNA_proteomics.R` script of the PloGO2 package (Jemma Wu and Dana Pascovici; GPL-2), described in Wu et al. (2020).

## References

1. Langfelder P, Horvath S. WGCNA: an R package for weighted correlation network analysis. *BMC Bioinformatics* 9: 559 (2008). https://doi.org/10.1186/1471-2105-9-559
2. Zhang B, Horvath S. A general framework for weighted gene co-expression network analysis. *Statistical Applications in Genetics and Molecular Biology* 4: Article17 (2005). https://doi.org/10.2202/1544-6115.1128
3. Langfelder P, Zhang B, Horvath S. Defining clusters from a hierarchical cluster tree: the Dynamic Tree Cut package for R. *Bioinformatics* 24(5): 719–720 (2008). https://doi.org/10.1093/bioinformatics/btm563
4. Wu JX, Pascovici D, Wu Y, Walker AK, Mirzaei M. Workflow for rapidly extracting biological insights from complex, multicondition proteomics experiments with WGCNA and PloGO2. *Journal of Proteome Research* 19(7): 2898–2906 (2020). https://doi.org/10.1021/acs.jproteome.0c00198

## Revision (October 2026)

The script was reviewed after publication. The main fixes:

- The eigenprotein plots now use the eigenproteins of the merged modules.
- The scale-free fit is now evaluated for a signed network, the type that is built. β = 8 is kept, so the modules are unchanged.
- Hub proteins are now restricted to the proteins assigned to the module.

The full list of changes is in the header of `WGCNA_CF.R`.

## License

GNU General Public License v2.0 (see [LICENSE](LICENSE)). The repository uses GPL-2 because it includes code from PloGO2, which is distributed under that license.
