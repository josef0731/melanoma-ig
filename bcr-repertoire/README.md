# Analysis of BCR repertoire features

This folder contains code for the high-throughput BCR repertoire analysis.

## Data

Data generated in this work (repertoire of intratumoural B cells) and those from healthy volunteers, Ebola convalescent patients and hospitalised COVID-19 patients can be found in <https://doi.org/10.5281/zenodo.7584227>.

The analysis also used bulk BCR repertoires from tonsilar B cells from [King et al. Science Immunology 2021](https://doi.org/10.1126/sciimmunol.abe6291). The data is downloaded from [ArrayExpress entry E-MTAB-8999](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-8999).

## Contents

`rep-features`: Analysis of repertoire features (CDR3, gene usage etc.). Files: 

* `cdr3_analysis_Oct22.Rmd`: Analysis of CDR3 characteristics (Kidera factors, amino acid usage etc.). *Note: R markdown format*.

* `gene_usage_Oct22.Rmd`: Analysis of VDJ germline gene usage. *Note: R markdown format*. 

* `isotype.R`: Counting isotype distributions.

* `unproductive.R`: Counting unproductive sequences.

`brepconvert`: Analysis of gene conversion events using [BrepConvert](https://github.com/Fraternalilab/BrepConvert/). See [Mallaby et al.](https://doi.org/10.1093/discim/kyad002) for the associated manuscript of the BrepConvert software. Files:

* `run_BrepConvert.R`: Code to read in sequence chunks and run the BrepConvert tool. *Note: for running a large number of BCR sequences over BrepConvert, we recommend separating the data into chunks and submit jobs on high-performance computing (HPC) resources in parallel for optimal runtime, if you have access to HPC resources.*

* `analyse_BrepConvert.R`: Analyse the BrepConvert results.
