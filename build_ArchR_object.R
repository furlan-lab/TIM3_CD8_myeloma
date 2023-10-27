setwd("/fh/fast/furlan_s/user/owalt/exhaustion/mm4")
RES_DIR  <- file.path("res")     # SPECIFY HERE
RMD_DIR  <- file.path("rmd")     # SPECIFY HERE
FIG_DIR <- file.path("figs")

suppressPackageStartupMessages({
  library(monocle3)
  library(m3addon)
  library(reticulate)
  library(viridis)
  library(openxlsx)  
  library(dplyr)
  library(Matrix)
  library(ggplot2)
  library(xfun)
  library(pals)
  library(RColorBrewer)
  library(Signac)
  library(Seurat)
  library(ggplot2)
  library(future)
  library(GenomicRanges)
  library(EnsDb.Mmusculus.v79)
  library(ComplexHeatmap)
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(JASPAR2020)
  library(TFBSTools)
  library(patchwork)
  library(stringr)
  library(ggrepel)
  library(ggsignif)
  library(ggpubr)
  library(ArchR)
  library(harmony)
  library(viridis)
  library(scCustomize)
  library(SeuratDisk)
  library(parallel)
})

set.seed(1234)
MO <- readRDS(file.path("MO_111521.rds"))
dyn.load('/app/software/ArrayFire/3.8.1/lib64/libaf.so.3')
library(RcppArrayFire)

Sys.setenv(MODULEPATH="/app/modules/all", MODULEPATH_ROOT="/app/modules/all", MODULESHOME="/app/lmod/lmod")
source("/app/lmod/lmod/init/R", echo=FALSE, max=Inf)
module("load", "MACS2/2.2.6-foss-2019b-Python-3.7.4")

samps<-c("MM3_2", "MM4_1", "MM4_2")

inputFiles <- sapply(1:length(samps), function(i){
  
  path<- file.path(samps[i],  "atac_fragments-Reformat.tsv.gz")
  path
})

addArchRGenome("mm10")

names(inputFiles)<- samps
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  minTSS = 3, 
  minFrags = 1000,
  subThreading = F,
  force = T
)

ArrowFiles<-paste0( "ArrowFiles/",paste0(samps, ".arrow"))
file.exists(ArrowFiles)

doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
  LSIMethod = 1
)

mm4 <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "/fh/fast/furlan_s/user/owalt/exhaustion/mm4",
  copyArrows = T #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)

saveArchRProject(mm4)
files<-paste(samps, "/filtered_feature_bc_matrix.h5", sep="")
Sys.setenv(HDF5_USE_FILE_LOCKING=FALSE,RHDF5_USE_FILE_LOCKING=FALSE)
seRNA <- import10xFeatureMatrix(input = files, names = samps)

mm4<-addGeneExpressionMatrix(input=mm4, seRNA=seRNA, force = T, excludeChr = "ChrY")
table(is.na(mm4$Gex_nUMI))
mm4 <- mm4[!is.na(mm4$Gex_nUMI)]
names(mm4@cellColData)

mm4 <- addIterativeLSI(
  ArchRProj = mm4, 
  clusterParams = list(
    resolution = 0.1, 
    sampleCells = 10000,
    n.start = 10
  ),
  saveIterations = FALSE,
  useMatrix = "GeneExpressionMatrix", 
  depthCol = "Gex_nUMI",
  varFeatures = 1500,
  firstSelection = "variable",
  binarize = F,
  name = "LSI_RNA",
  dimsToUse = 1:15,
  force = T
  
)

mm4 <- addIterativeLSI(
  ArchRProj = mm4,
  useMatrix = "TileMatrix", 
  name = "ATAC_LSI", 
  iterations = 4, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.1, 0.2, 0.4), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 10000, 
  dimsToUse = 2:15,
  force =T
)

mm4 <- addHarmony(
  ArchRProj = mm4,
  reducedDims = "ATAC_LSI",
  name = "Harmony",
  groupBy = "Sample",
  scaleDims = T,
  force = T,
  corCutOff = 0.6
)

mm4 <- addCombinedDims(mm4, reducedDims = c("LSI_RNA", "Harmony"), name =  "mo")
mm4<- addImputeWeights(mm4, reducedDims = "mo")
mm4<- addUMAP(mm4, reducedDims = "mo", force = T, minDist = 0.4, corCutOff = 0.5, scaleDims = T, nNeighbors = 15)
plotEmbedding(mm4, reduction = "UMAP")

MO<-readRDS("MO_111521.rds")

DimPlot(MO, group.by = "sub.cluster", reduction = "wnn.umap.harmony")
mm4@embeddings$mo_UMAP<-mm4@embeddings$UMAP

umap<- data.frame(MO@reductions$wnn.umap.harmony@cell.embeddings)

rownames(umap)<- gsub("MM3_2_", "MM3_2#",rownames(umap))
rownames(umap)<- gsub("MM4_2_", "MM4_2#", rownames(umap))
rownames(umap) <- gsub("MM4_1_", "MM4_1#", rownames(umap))

umap

colnames(mm4)

umap<- umap[which(rownames(umap) %in% rownames(mm4@cellColData)),]

umap<- umap[order(match(rownames(umap), rownames(mm4@cellColData))),]

names(umap)<- names(mm4@embeddings@listData$mo_UMAP@listData$df)

mm4@embeddings@listData$mo_UMAP@listData$df <- umap

df<-MO@meta.data
rownames(df) <- gsub("MM3_2_", "MM3_2#", rownames(df))
rownames(df) <- gsub("MM4_2_", "MM4_2#", rownames(df))
rownames(df) <- gsub("MM4_1_", "MM4_1#", rownames(df))

mm4<-mm4[which(rownames(mm4@cellColData) %in% rownames(df)),]

meta<- data.frame(mm4@cellColData)

df<- df[which(rownames(df) %in% rownames(mm4@cellColData)),]

df<- df[order(match(rownames(df), rownames(mm4@cellColData))),]

mm4$cluster<- df$sub.cluster

mm4$state <- mm4$Sample
mm4$state[mm4$Sample == "MM3_2"]<-"3_Relapsed"
mm4$state[mm4$Sample == "MM4_1"]<-"2_Controlled"
mm4$state[mm4$Sample == "MM4_2"]<-"1_Free"

mm4$cell_type<-mm4$cluster

mm4$cell_type[mm4$cluster == "2_0"]<-"TEX_1"
mm4$cell_type[mm4$cluster == "1_0"]<-"TEX_2"
mm4$cell_type[mm4$cluster == "1_1"]<-"TEX_3"
mm4$cell_type[mm4$cluster == "0"]<-"Naive"
mm4$cell_type[mm4$cluster == "4"]<-"TEX_Cycling"
mm4$cell_type[mm4$cluster == "2_3"]<-"TEM_Cytotoxic"
mm4$cell_type[mm4$cluster == "2_1"]<-"TEM_Il18r+"
mm4$cell_type[mm4$cluster == "2_2"]<-"TEM"
mm4$cell_type[mm4$cluster == "3"]<-"TCM"

mm4$group <- mm4$cell_type
mm4$group[mm4$group == "TEX_Cycling"]<-"Cycling"
mm4$group[grep("^TEX", mm4$group)]<-"TEX"

mm4$group[grep("^TEM", mm4$group)]<-"TEM"
mm4<- addImputeWeights(mm4, reducedDims = "mo")
plotEmbedding(mm4, name = "group", embedding = "mo_UMAP")

saveArchRProject(mm4)

pathToMacs2 <- findMacs2()

mm4<- addGroupCoverages(mm4, groupBy = "cell_type")

mm4 <- addReproduciblePeakSet(
  ArchRProj = mm4, 
  groupBy = "cell_type", 
  pathToMacs2 = pathToMacs2
)

mm4 <- addPeakMatrix(mm4)
getAvailableMatrices(mm4)
saveArchRProject(mm4)

mm4 <- addBgdPeaks(mm4)
mm4 <- addMotifAnnotations(ArchRProj = mm4, motifSet = "cisbp", name = "Motif")

markersPeaks <- getMarkerFeatures(
  ArchRProj = mm4, 
  useMatrix = "PeakMatrix", 
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

enrichMotifs <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = mm4,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)

heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 7, transpose = TRUE)
ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")

mm4<-addDeviationsMatrix(mm4, force = T)
saveArchRProject(mm4)
