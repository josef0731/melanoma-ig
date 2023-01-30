library(plyr)
library(BrepConvert)

library(argparse, quietly = TRUE, verbose = FALSE, warn.conflicts = FALSE)

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument("-f", "--fasta", type="character",
                    help="filepath of input fasta file. Relative to working directory")
parser$add_argument("-v", "--vgenes", type="character",
                    help="filepath of input fasta containing V gene sequences. Relative to working directory.")
parser$add_argument("-g", "--gene", type="character",
                    help="gene of interest")
args <- parser$parse_args()


ref_seq <- Biostrings::readDNAStringSet(args$vgenes)
gene <- args$gene
#test_data <- readRDS("../HV_EB_COV_MEL_comb_Oct22.RDS")
vdj <- Biostrings::readDNAStringSet(args$fasta)
rep_n <- unlist(strsplit(args$fasta, split = "_"))
rep_n <- gsub(".fasta", "", rep_n[length(rep_n)])

cat(paste0(gene, " ...\n"))
# select sequences of the gene as 'functional germline'
functional <- ref_seq[ which(grepl(paste0("|", gene), names(ref_seq), fixed = TRUE))]
pseudogenes <- ref_seq[ which(!grepl(paste0("|", gene), names(ref_seq), fixed = TRUE))]
pseudogenes <- pseudogenes[ which(!grepl(paste0("|", gene), names(pseudogenes), fixed = TRUE))]
if( !file.exists( paste0("functional_", gene, ".fasta") ) ){
  Biostrings::writeXStringSet(functional, paste0("functional_", gene, ".fasta"))
  Biostrings::writeXStringSet(pseudogenes, paste0("pseudogenes_", gene, ".fasta"))
}

repertoire <- vdj
# remove those with >= 50 blanks at the front (those are most likely incomplete sequences)
repertoire <- repertoire[ which( ! grepl("^[\\.]{50}", repertoire ) ) ]
repertoire <- repertoire[ which( repertoire != "" ) ]
# optional blacklist to manually remove those which can't go through BrepConvert
blacklist <- c() # add Seq_ID here
repertoire <- repertoire[ which( ! names(repertoire) %in% blacklist ) ]
repertoire <- as.character(repertoire)

#select repertoire sequences with this Vgene
#repertoire <- test_data[which(test_data$Vfamily == gene & test_data$Seq_ID %in% names(vdj)), "Seq_ID"]
#repertoire <- as.character(vdj[repertoire])
cat(paste0("Total n = ", length(repertoire), " sequences ...\n"))
results <- BrepConvert::batchConvertAnalysis(functional = paste0("functional_", gene, ".fasta"),
                                             pseudogene = paste0("pseudogenes_", gene, ".fasta"),
                                             repertoire = repertoire,
                                             blat_exec = system.file("exe/blat", package = "BrepConvert"),
                                             convertAA = FALSE)
saveRDS(results, paste0(gene, "_BrepConvert_results", rep_n, ".rds"))
