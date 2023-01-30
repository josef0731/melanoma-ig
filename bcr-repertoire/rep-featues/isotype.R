#_________________________________
# isotype distribution
all_seqs <- readRDS("HV_EB_COV_MEL_comb_Oct22.RDS")

isotype_distribution <- data.frame(
  table(all_seqs[which(all_seqs$SampleType == "melanoma"), 
                 "PatientID"], 
        all_seqs[which(all_seqs$SampleType == "melanoma"), 
                 "Subclass"])
)
isotype_distribution <- isotype_distribution[
  which(isotype_distribution$Var2 %in% 
          paste0("Ig", c("M", "D", "G3", "G1", "A1", 
                         "G2", "G4", "E", "A2"))), 
]
library(plyr)
isotype_distribution <- merge(
  isotype_distribution, 
  ddply(isotype_distribution, "Var1", summarise, 
        Freq_total = sum(Freq)), by = "Var1", 
  suffixes = c("", "_total")
)
isotype_distribution$prop <- isotype_distribution$Freq / isotype_distribution$Freq_total
library(ggplot2)
ggplot(isotype_distribution, 
       aes(y = Var1, x = Var2, fill = prop, 
           label = scales::percent(prop, accuracy = 0.01))) + 
  geom_tile() + cowplot::theme_cowplot() + 
  geom_text(aes(color = (prop < 0.5))) + 
  scale_fill_gradient2(name = "% repertoire", 
                       labels = scales::percent) + 
  scale_color_manual(values = c("TRUE" = "black", 
                                "FALSE" = "white")) + 
  xlab("") + ylab("Tumour") + guides(color = "none")

