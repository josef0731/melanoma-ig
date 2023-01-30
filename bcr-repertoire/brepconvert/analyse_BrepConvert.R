#_________________________________
# BrepConvert to look for receptor revision

#_________________________________
# Heavy chain 
brepconvert_results <- list(
  list.files(pattern = "IGHV[1-6].*results.*.rds",
             path = "BrepConvert/BrepConvert_healthy/",
             full.names = TRUE),
  list.files(pattern = "IGHV[1-6].*results.*.rds",
             path = "BrepConvert/BrepConvert/BrepConvert_tumour/",
             full.names = TRUE)
)
brepconvert_results <- list(
  do.call("rbind", lapply(brepconvert_results[[1]], readRDS)),
  do.call("rbind", lapply(brepconvert_results[[2]], readRDS))
)

# remove events without identifiable donor gene
brepconvert_results[[1]] <- brepconvert_results[[1]][which(brepconvert_results[[1]]$gene != "NA"), ]
brepconvert_results[[2]] <- brepconvert_results[[2]][which(brepconvert_results[[2]]$gene != "NA"), ]

# remove events entirely within cdr regions (since we can't be
# certain if they are just mutations)
cdr_pos <- list( # IMGT definitions
  "CDR1" = c(26 * 3 + 1, 38 * 3),
  "CDR2" = c(55 * 3 + 1, 66 * 3),
  "CDR3" = c(104 * 3 + 1, 117 * 3)
)
brepconvert_results[[1]]$cdr <- apply(brepconvert_results[[1]][, c("start", "end")], MARGIN = 1, function(x){
  sum(sapply(cdr_pos, function(y){
    x[1] >= y[1] & x[2] <= y[2]
  })) > 0
})
brepconvert_results[[2]]$cdr <- apply(brepconvert_results[[2]][, c("start", "end")], MARGIN = 1, function(x){
  sum(sapply(cdr_pos, function(y){
    x[1] >= y[1] & x[2] <= y[2]
  })) > 0
})
brepconvert_results[[1]] <- brepconvert_results[[1]][which(!brepconvert_results[[1]]$cdr), ]
brepconvert_results[[2]] <- brepconvert_results[[2]][which(!brepconvert_results[[2]]$cdr), ]

# size of events
brepconvert_results[[1]]$size <- apply(brepconvert_results[[1]][, c("start", "end")], MARGIN = 1, function(x) x[2] - x[1] + 1)
brepconvert_results[[2]]$size <- apply(brepconvert_results[[2]][, c("start", "end")], MARGIN = 1, function(x) x[2] - x[1] + 1)

all_seqs <- readRDS("HV_EB_COV_MEL_comb_Oct22.RDS")
all_seqs2 <- list(
  all_seqs[which(all_seqs$SampleType == "Healthy"), ],
  all_seqs[which(all_seqs$SampleType == "melanoma"), ]
)
all_seqs2[[1]]$converted <- (all_seqs2[[1]]$Seq_ID %in% brepconvert_results[[1]][brepconvert_results[[1]]$size > 30, "SeqID"])
all_seqs2[[2]]$converted <- (all_seqs2[[2]]$Seq_ID %in% brepconvert_results[[2]][brepconvert_results[[2]]$size > 30, "SeqID"])

library(plyr)
converted_counts <- list(
  ddply(all_seqs2[[1]][which(grepl("^productive", all_seqs2[[1]]$V.DOMAIN.Functionality) & !is.na(all_seqs2[[1]]$Subclass)), ], 
        c("PatientID", "Subclass", "converted"),
        nrow, .drop = FALSE),
  ddply(all_seqs2[[2]][which(grepl("^productive", all_seqs2[[2]]$V.DOMAIN.Functionality) & !is.na(all_seqs2[[2]]$Subclass)), ], 
        c("PatientID", "Subclass", "converted"),
        nrow, .drop = FALSE)
)
converted_counts[[1]] <- converted_counts[[1]][which(converted_counts[[1]]$Subclass %in% c("IgM", "IgG3", "IgG1", "IgA1",
                                                                                           "IgG2", "IgG4", "IgE", "IgA2")), ]
converted_counts[[2]] <- converted_counts[[2]][which(converted_counts[[2]]$Subclass %in% c("IgM", "IgG3", "IgG1", "IgA1",
                                                                                           "IgG2", "IgG4", "IgE", "IgA2")), ]
converted_counts[[1]]$SampleType <- "Healthy"
converted_counts[[2]]$SampleType <- "melanoma"
colnames(converted_counts[[1]]) <- c("PatientID", "Subclass", "converted", "n", "SampleType")
colnames(converted_counts[[2]]) <- c("PatientID", "Subclass", "converted", "n", "SampleType")
converted_counts[[1]] <- converted_counts[[1]][, c("PatientID", "Subclass", "converted", "SampleType", "n")]
converted_counts[[2]] <- converted_counts[[2]][, c("PatientID", "Subclass", "converted", "SampleType", "n")]
converted_counts <- rbind(converted_counts[[1]], converted_counts[[2]])
converted_counts <- merge(
  converted_counts, 
  ddply(converted_counts, c("PatientID", "SampleType", "Subclass"),
        summarise, total = sum(n))
)
converted_counts$proportion <- converted_counts$n / converted_counts$total

library(ggplot2)
hv <- ggplot(converted_counts[which(converted_counts$converted & 
                                      ! converted_counts$Subclass %in% c("IgE", "IgG4")), ], 
             aes(x = Subclass, y = proportion, fill = SampleType)) + 
  geom_boxplot(outlier.shape = NA, 
               position = position_dodge(width = 1), color = "black") + 
  geom_point(position = position_dodge(width = 1), size = 1.3) + 
  cowplot::theme_cowplot() + 
  scale_y_continuous(labels = scales::percent, 
                     name = "% sequences with\nreplacements > 30 nt") +
  ggpubr::stat_compare_means(label = "p.signif") +
  scale_fill_manual(values = c("Healthy" = "#d6d6d6", "melanoma" = "#ff9300"), 
                    labels = c("HV", "MP"), name = "")
print(hv)
write.table(converted_counts[which(converted_counts$converted ), ], 
            "proportion_BrepConvert_heavy.csv",
            sep = ",", row.names = FALSE, col.names = TRUE)

#_________________________________
# BrepConvert light chains
brepconvert_results <- list(
  list.files(pattern = "IG[KL]V.*results.*.rds",
             path = "BrepConvert/BMLight/",
             full.names = TRUE),
  list.files(pattern = "IG[KL]V.*results.*.rds",
             path = "BrepConvert/BrepConvert_tumour/",
             full.names = TRUE)
)
brepconvert_results <- list(
  do.call("rbind", lapply(brepconvert_results[[1]], readRDS)),
  do.call("rbind", lapply(brepconvert_results[[2]], readRDS))
)

# remove events without identifiable donor gene
brepconvert_results[[1]] <- brepconvert_results[[1]][which(brepconvert_results[[1]]$gene != "NA"), ]
brepconvert_results[[2]] <- brepconvert_results[[2]][which(brepconvert_results[[2]]$gene != "NA"), ]

# remove events entirely within cdr regions (since we can't be
# certain if they are just mutations)
cdr_pos <- list( # IMGT definitions
  "CDR1" = c(26 * 3 + 1, 38 * 3),
  "CDR2" = c(55 * 3 + 1, 66 * 3),
  "CDR3" = c(104 * 3 + 1, 117 * 3)
)
brepconvert_results[[1]]$cdr <- apply(brepconvert_results[[1]][, c("start", "end")], MARGIN = 1, function(x){
  sum(sapply(cdr_pos, function(y){
    x[1] >= y[1] & x[2] <= y[2]
  })) > 0
})
brepconvert_results[[2]]$cdr <- apply(brepconvert_results[[2]][, c("start", "end")], MARGIN = 1, function(x){
  sum(sapply(cdr_pos, function(y){
    x[1] >= y[1] & x[2] <= y[2]
  })) > 0
})
brepconvert_results[[1]] <- brepconvert_results[[1]][which(!brepconvert_results[[1]]$cdr), ]
brepconvert_results[[2]] <- brepconvert_results[[2]][which(!brepconvert_results[[2]]$cdr), ]

# size of events
brepconvert_results[[1]]$size <- apply(brepconvert_results[[1]][, c("start", "end")], MARGIN = 1, function(x) x[2] - x[1] + 1)
brepconvert_results[[2]]$size <- apply(brepconvert_results[[2]][, c("start", "end")], MARGIN = 1, function(x) x[2] - x[1] + 1)

all_seqs2 <- list(
  read.csv("BM_Master_Light_Useasreftrue_140325_withPepStats_8.csv",
           stringsAsFactors = FALSE),
  all_seqs[which(all_seqs$SampleType == "melanoma"), ]
)
all_seqs2[[1]]$converted <- (all_seqs2[[1]]$Seq_ID %in% brepconvert_results[[1]][brepconvert_results[[1]]$size > 30, "SeqID"])
all_seqs2[[2]]$converted <- (all_seqs2[[2]]$Seq_ID %in% brepconvert_results[[2]][brepconvert_results[[2]]$size > 30, "SeqID"])
all_seqs2[[1]]$Subclass <- substr(all_seqs2[[1]]$Vfamily, 1, 3)

library(plyr)
converted_counts <- list(
  ddply(all_seqs2[[1]][which(grepl("^productive", all_seqs2[[1]]$Functionality) & !is.na(all_seqs2[[1]]$Subclass)), ], 
        c("Patient_ID", "Subclass", "converted"),
        nrow),
  ddply(all_seqs2[[2]][which(grepl("^productive", all_seqs2[[2]]$V.DOMAIN.Functionality) & !is.na(all_seqs2[[2]]$Subclass)), ], 
        c("PatientID", "Subclass", "converted"),
        nrow)
)
converted_counts[[2]] <- converted_counts[[2]][which(converted_counts[[2]]$Subclass %in% c("IgK", "IgL")), ]
converted_counts[[2]]$Subclass <- toupper(converted_counts[[2]]$Subclass)
converted_counts[[1]]$SampleType <- "Healthy"
converted_counts[[2]]$SampleType <- "melanoma"
colnames(converted_counts[[1]]) <- c("PatientID", "Subclass", "converted", "n", "SampleType")
colnames(converted_counts[[2]]) <- c("PatientID", "Subclass", "converted", "n", "SampleType")
converted_counts[[1]] <- converted_counts[[1]][, c("PatientID", "Subclass", "converted", "SampleType", "n")]
converted_counts[[2]] <- converted_counts[[2]][, c("PatientID", "Subclass", "converted", "SampleType", "n")]
converted_counts <- rbind(converted_counts[[1]], converted_counts[[2]])
converted_counts <- merge(
  converted_counts, 
  ddply(converted_counts, c("PatientID", "SampleType", "Subclass"),
        summarise, total = sum(n))
)
converted_counts$proportion <- converted_counts$n / converted_counts$total
lv <- ggplot(converted_counts[which(converted_counts$converted), ], 
             aes(x = Subclass, y = proportion, fill = SampleType)) + 
  geom_boxplot(outlier.shape = NA, 
               position = position_dodge(width = 1), color = "black") + 
  geom_point(position = position_dodge(width = 1), size = 1.3) + 
  cowplot::theme_cowplot() + xlab("") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1),
                     name = "% sequences with\nrearrangement > 30 nt") +
  ggpubr::stat_compare_means(label = "p.signif", 
                             symnum.args = list(cutpoints = c(0, 0.01, 1),
                                                symbols = c("*", "ns"))) +
  scale_fill_manual(values = c("Healthy" = "grey50", "melanoma" = "#ff9300"), 
                    labels = c("HV", "MP"), name = "")
print(lv)

cowplot::plot_grid(hv + ggtitle("Heavy chain"), 
                   lv + ggtitle("Light chain"), 
                   nrow = 1, rel_widths = c(3.5, 2.5), 
                   align = "h", axis = "tb")

write.table(converted_counts[which(converted_counts$converted ), ], 
            "proportion_BrepConvert_light.csv",
            sep = ",", row.names = FALSE, col.names = TRUE)
