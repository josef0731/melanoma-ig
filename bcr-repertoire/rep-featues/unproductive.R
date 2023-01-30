#_________________________________
# find % of unproductive sequences
tonsil_repseq <- list.files(path = ".", 
                            full.names = TRUE, 
                            pattern = "REPseq")
tonsil_repseq <- lapply(tonsil_repseq, read.table, sep = "\t", 
                        header= TRUE, stringsAsFactors = FALSE)
tonsil_repseq <- do.call("rbind", tonsil_repseq)

all_seqs <- readRDS("HV_EB_COV_MEL_comb_Oct22.RDS")

library(plyr)
counts <- list(
  ddply(all_seqs, c("PatientID", "SampleType", 
                    "V.DOMAIN.Functionality", "Subclass"),
        nrow),
  ddply(tonsil_repseq[which(!duplicated(tonsil_repseq$SEQUENCE_VDJ)), ], 
        c("DONOR", "FUNCTIONAL", "ISOTYPE"),
        nrow)
)
counts[[1]] <- counts[[1]][which(grepl("productive", counts[[1]][, "V.DOMAIN.Functionality"])), ]
counts[[1]] <- counts[[1]][which(counts[[1]]$Subclass %in% c("IgM", "IgG3", "IgG1", "IgA1",
                                                             "IgG2", "IgG4", "IgE", "IgA2")), ]
counts[[2]] <- counts[[2]][which(counts[[2]]$ISOTYPE %in% c("IgM", "IgG3", "IgG1", "IgA1",
                                                            "IgG2", "IgG4", "IgE", "IgA2")), ]
counts[[1]]$productive <- (grepl("^productive", counts[[1]][, "V.DOMAIN.Functionality"]))
counts[[2]]$SampleType <- "Tonsils_GC"
colnames(counts[[1]]) <- c("PatientID", "SampleType", "productive", "Subclass", "n", "FUNCTIONAL")
colnames(counts[[2]]) <- c("PatientID", "FUNCTIONAL", "Subclass", "n", "SampleType")
counts[[1]] <- counts[[1]][, c("PatientID", "SampleType", "Subclass", "FUNCTIONAL", "n")]
counts[[2]] <- counts[[2]][, c("PatientID", "SampleType", "Subclass", "FUNCTIONAL", "n")]
counts <- rbind(counts[[1]], counts[[2]])
counts <- merge(
  counts, ddply(counts, c("PatientID", "SampleType", "Subclass"),
                summarise, total = sum(n))
)
counts$proportion <- counts$n / counts$total
counts <- ddply(counts, c("PatientID", "SampleType", "Subclass", "FUNCTIONAL"),
                summarise, n = sum(n), total = unique(total), proportion = sum(proportion))
counts$SampleType <- factor(counts$SampleType,
                            levels = c("Healthy", "melanoma",
                                       "EB", "COVID", "Tonsils_GC"),
                            labels = c("HV", "MP", "EB", "Covid",
                                       "Tonsils"))
write.table(counts[which(counts$FUNCTIONAL == "FALSE" ), ], 
            "proportion_unproductive.csv",
            sep = ",", row.names = FALSE, col.names = TRUE)

stat_test <- ggpubr::compare_means(
  proportion ~ SampleType, 
  data = counts[which(counts$FUNCTIONAL == "FALSE" & #counts$SampleType %in% c("HV", "MP") &
                        ! counts$Subclass %in% c("IgE", "IgG4")), ], 
  group.by = "Subclass", ref.group = "HV", p.adjust.method = "bonferroni"
)
stat_test$p.signif <- sapply(stat_test$p.adj, function(x){
  if(x < 0.0001) return("****")
  if(x < 0.001) return("***")
  if(x < 0.01) return("**")
  if(x < 0.05) return("*")
  return("ns")
})
library(ggplot2)
svg("unproductive_percentage.svg",
    width = 6.2, height = 3.2)
ggplot(counts[which(counts$FUNCTIONAL == "FALSE" & 
                      ! counts$Subclass %in% c("IgE", "IgG4")), ], 
       aes(x = Subclass, y = proportion, fill = SampleType)) + 
  geom_boxplot(outlier.shape = NA, 
               position = position_dodge(width = 0.8), color = "black") + 
  geom_point(position = position_dodge(width = 0.8), size = 0.6) + 
  scale_y_continuous(labels= scales::percent, 
                     name = "% sequences unproductive") +
  scale_fill_manual(values = c("HV" = "#d6d6d6", "MP" = "#ff9300", 
                               "EB" = "#fffb00", "Covid" = "#00f900",
                               "Tonsils" = "#E76BF3"), name = "") +
  cowplot::theme_cowplot()
dev.off()

