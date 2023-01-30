# Analysis of gene conversion events 

Using [BrepConvert](https://github.com/Fraternalilab/BrepConvert/). 

See [Mallaby et al.](https://doi.org/10.1093/discim/kyad002) for the associated manuscript of the BrepConvert software. 

Files:

* `run_BrepConvert.R`: Code to read in sequence chunks and run the BrepConvert tool. *Note: for running a large number of BCR sequences over BrepConvert, we recommend separating the data into chunks and submit jobs on high-performance computing (HPC) resources in parallel for optimal runtime, if you have access to HPC resources.*

* `analyse_BrepConvert.R`: Analyse the BrepConvert results.
