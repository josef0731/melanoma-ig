library(CATALYST)
library(cowplot)
library(flowCore)
library(diffcyt)
library(scater)
library(SingleCellExperiment)
library(readxl)
library(uwot)
library(reshape2) 

input_dir <- "C:\\Path_to_directory\\Phenograph_Output\\"
output_dir <- paste0(input_dir)

data_files = choose.files(default=paste0(input_dir, "*.fcs"),filters = Filters[c(".fcs")])
metadata_file <- choose.files(default=paste0(input_dir, "*.xlsx"),filters = Filters[c(".xlxs")])
panel_data_file <- choose.files(default=paste0(input_dir, "*.xlsx"),filters = Filters[c(".xlsx")])

md <- read_excel(metadata_file)
panel <- read_excel(panel_data_file)

fs <- read.flowSet(data_files, transformation = TRUE, truncate_max_range = FALSE)
sce <- prepData(fs, panel, md, features = panel$fcs_colname)

#When merging
merging_table_file <- choose.files(default=paste0(input_dir, "*.xlsx"),filters = Filters[c(".xlsx")])
mt <- read_excel(merging_table_file)

foo <- sce
foo$cluster_id <- sce@int_colData$Cluster_ID
metadata(foo)$cluster_codes <- data.frame(
  custom = factor(c(1:max( sce@int_colData$Cluster_ID))))

#When merging
foo <- mergeClusters(foo, k = "custom", table = mt, id = "merged") 

names(cluster_codes(foo))
table(cluster_ids(foo, "merged")) 
delta_area(foo)

hm1 <- plotExprHeatmap(foo, features = "type", 
                       by = "cluster_id", k = "custom", 
                       scale = "never", q = 0.01, perc = TRUE, 
                       col_dend = FALSE, col_clust = FALSE, col_anno = TRUE,
                       row_clust = FALSE, fun = "median",
                       row_anno="condition", bars = TRUE, bin_anno = TRUE)

hm3 <- plotExprHeatmap(foo, features = "type",
                       by = "cluster_id", k = "custom",
                       scale = "last", q = 0.01, perc = TRUE, 
                       col_dend = FALSE, col_clust = FALSE, col_anno = TRUE,
                       row_clust = FALSE, fun = "median", row_anno="condition", bars = TRUE, bin_anno = TRUE)

write.csv(hm1@matrix,paste0(output_dir,"Heatmap_Median_Table_All.csv"))
#For merging
hm2 <- plotExprHeatmap(foo, features = "type",
                       by = "cluster_id", k = "merged",
                       scale = "last", q = 0.01, perc = TRUE, 
                       col_dend = FALSE, col_clust = FALSE, col_anno = TRUE,
                       row_clust = FALSE, fun = "median",row_anno="condition", bars = TRUE, bin_anno = TRUE)

write.csv(hm2@matrix,paste0(output_dir,"Heatmap_Table_Phenograph_Median_Merged_Norm.csv"))