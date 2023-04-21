##### script for installing and loading packages ####

# install Catalyst

install.packages("BiocManager")
BiocManager::install("CATALYST")


# Load packages

library(Cairo)
library (Matrix)
library(CATALYST)
library(carData)
library(cowplot)
library(flowCore)
library(diffcyt)
library(scater)
library(SingleCellExperiment)
library(umap)
library(uwot)
library(destiny)
library(readxl)
#### set working directory
getwd()


# upload anything with FCS in the wd, read FCS files flowset and label as silvia_fs
fcs_files <- list.files(pattern = ".fcs")
MvHV_fs <- read.flowSet(fcs_files, transformation = FALSE, truncate_max_range = FALSE)


## Read in panel as excel file and metadata
##
panel <- read_excel("Cytof_panel.xlsx")

## Read in metadata excel ##
md <- read_excel("metadata.xlsx")


# spot check that all panel columns are in the flowSet object - checks that columns in FCS files vs panel match up
all(panel$fcs_colname %in% colnames(MvHV_fs))

# specify levels for conditions & sample IDs to assure desired ordering - storing  more info about FCS files - grouping things by tox vs no tox
md$condition <- factor(md$condition, levels = c("Healthy Volunteer", "Melanoma"))        
md$sample_id <- factor(md$sample_id,                                    
                       levels = md$sample_id[order(md$condition)])  


# construct SingleCellExperiment   - command for bringing together data matrix                                      
sce <- prepData(MvHV_fs, panel, md, features = panel$fcs_colname)


# cell population identification / clustering by phenotype 
set.seed(1220)
sce <- cluster(sce, features = type_markers(sce),
               xdim = 10, ydim = 10, maxK = 25, seed = 4321)


# visualisation of clustering - plot heatmap 
plotClusterHeatmap(sce,
                   hm2 = NULL, k = "meta25", m = NULL, cluster_anno = TRUE, draw_freqs = TRUE)


# diffusion map code 
set.seed(1234)
sce <- runDR(sce, "DiffusionMap", cells = 30000, features = "type")

#PLOT DR MAPS
plotDR(sce, "DiffusionMap", color_by = "meta25")

#Split by Condition
plotDR(sce, "DiffusionMap", color_by = "meta25", facet_by = "condition")

