Compute per-barcode binding to vaccinated mouse sera samples
================
Tyler Starr
5/12/2022

This notebook reads in per-barcode counts from `count_variants.ipynb`
for sera-binding titration experiments, computes functional scores for
RBD binding values via delta-AUC metrics, and does some basic QC on
variant binding functional scores.

``` r
#list of packages to install/load
packages = c("yaml","data.table","tidyverse","gridExtra")
#install any packages not already installed
installed_packages <- packages %in% rownames(installed.packages())
if(any(installed_packages == F)){
  install.packages(packages[!installed_packages],
                   lib=c(paste("/uufs/chpc.utah.edu/common/home/",Sys.getenv("USER"),"/RLibs/",Sys.getenv("R_VERSION"),sep="")),
                   repos=c("http://cran.us.r-project.org"))
}
#load packages
invisible(lapply(packages, library, character.only=T))

knitr::opts_chunk$set(echo = T)
knitr::opts_chunk$set(dev.args = list(png = list(type = "cairo")))

#read in config file
config <- read_yaml("config.yaml")

#make output directory
if(!file.exists(config$mAb_EC50_dir)){
  dir.create(file.path(config$mAb_EC50_dir))
}
```

Session info for reproducing environment:

``` r
sessionInfo()
```

    ## R version 4.1.3 (2022-03-10)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Rocky Linux 8.8 (Green Obsidian)
    ## 
    ## Matrix products: default
    ## BLAS/LAPACK: /uufs/chpc.utah.edu/sys/spack/linux-rocky8-nehalem/gcc-8.5.0/intel-oneapi-mkl-2021.4.0-h43nkmwzvaltaa6ii5l7n6e7ruvjbmnv/mkl/2021.4.0/lib/intel64/libmkl_rt.so.1
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] gridExtra_2.3     forcats_0.5.1     stringr_1.4.0     dplyr_1.0.8      
    ##  [5] purrr_0.3.4       readr_2.1.2       tidyr_1.2.0       tibble_3.1.6     
    ##  [9] ggplot2_3.4.1     tidyverse_1.3.1   data.table_1.14.2 yaml_2.3.5       
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_1.1.2 xfun_0.30        haven_2.4.3      colorspace_2.0-3
    ##  [5] vctrs_0.5.2      generics_0.1.2   htmltools_0.5.2  utf8_1.2.2      
    ##  [9] rlang_1.0.6      pillar_1.7.0     glue_1.6.2       withr_2.5.0     
    ## [13] DBI_1.1.2        dbplyr_2.1.1     modelr_0.1.8     readxl_1.3.1    
    ## [17] lifecycle_1.0.3  munsell_0.5.0    gtable_0.3.0     cellranger_1.1.0
    ## [21] rvest_1.0.2      evaluate_0.15    knitr_1.37       tzdb_0.2.0      
    ## [25] fastmap_1.1.0    fansi_1.0.2      broom_0.7.12     Rcpp_1.0.11     
    ## [29] backports_1.4.1  scales_1.2.1     jsonlite_1.8.7   fs_1.5.2        
    ## [33] hms_1.1.1        digest_0.6.29    stringi_1.7.6    grid_4.1.3      
    ## [37] cli_3.6.0        tools_4.1.3      magrittr_2.0.2   crayon_1.5.0    
    ## [41] pkgconfig_2.0.3  ellipsis_0.3.2   xml2_1.3.3       reprex_2.0.1    
    ## [45] lubridate_1.8.0  rstudioapi_0.13  assertthat_0.2.1 rmarkdown_2.13  
    ## [49] httr_1.4.7       R6_2.5.1         compiler_4.1.3

## Setup

First, we will read in metadata on our sort samples, the table giving
number of reads of each barcode in each of the sort bins, and the
barcode-variant lookup tables, and merge these tables together.

``` r
#read dataframe with list of barcode runs
barcode_runs <- read.csv(file=config$barcode_runs,stringsAsFactors=F); barcode_runs <- subset(barcode_runs, select=-c(R1))

#read file giving count of each barcode in each sort partition
counts <- data.table(read.csv(file=config$variant_counts_file,stringsAsFactors=F))

#read in barcode-variant lookup tables
dt <- data.table(read.csv(file=config$codon_variant_table_file_lib47,stringsAsFactors=F))
setkey(dt,barcode,library)

dt <- merge(counts, dt, by=c("library","barcode")); rm(counts)

samples_268C61 <- data.frame(sample=sort(unique(paste(rep("268C61",6),formatC(barcode_runs[barcode_runs$sample_type=="268C61","concentration"], width=2,flag="0"),sep="_"))),conc=c(10000, 400, 16, 0.64, 0.0256, 0))

samples_268C183 <- data.frame(sample=sort(unique(paste(rep("268C183",6),formatC(barcode_runs[barcode_runs$sample_type=="268C183","concentration"], width=2,flag="0"),sep="_"))),conc=c(10000, 400, 16, 0.64, 0.0256, 0))

samples_268C185 <- data.frame(sample=sort(unique(paste(rep("268C185",6),formatC(barcode_runs[barcode_runs$sample_type=="268C185","concentration"], width=2,flag="0"),sep="_"))),conc=c(10000, 400, 16, 0.64, 0.0256, 0))

samples_268C3 <- data.frame(sample=sort(unique(paste(rep("268C3",6),formatC(barcode_runs[barcode_runs$sample_type=="268C61","concentration"], width=2,flag="0"),sep="_"))),conc=c(10000, 400, 16, 0.64, 0.0256, 0))
```

Convert from Illumina read counts to estimates of the number of cells
that were sorted into a bin, and add some other useful information to
our data frame.

``` r
#for each bin, normalize the read counts to the observed ratio of cell recovery among bins
for(i in 1:nrow(barcode_runs)){
  lib <- as.character(barcode_runs$library[i])
  bin <- as.character(barcode_runs$sample[i])
  ratio <- sum(dt[library==lib & sample==bin,"count"])/barcode_runs$number_cells[i]
  if(ratio<1){ #if there are fewer reads from a FACS bin than cells sorted
    dt[library==lib & sample==bin, count.norm := as.numeric(count)] #don't normalize cell counts, make count.norm the same as count
    print(paste("reads < cells for",lib,bin,", un-normalized (ratio",ratio,")")) #print to console to inform of undersampled bins
  }else{
    dt[library==lib & sample==bin, count.norm := as.numeric(count/ratio)] #normalize read counts by the average read:cell ratio, report in new "count.norm" column
    print(paste("read:cell ratio for",lib,bin,"is",ratio))
  }
}
```

    ## [1] "read:cell ratio for lib47_SARSr-wts 268C3_01_bin1 is 1.90162924713764"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C3_01_bin2 is 6.56682863590299"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C3_01_bin3 is 6.50915001203949"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C3_01_bin4 is 15.2117138742199"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C3_02_bin1 is 2.49749149432318"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C3_02_bin2 is 8.63144492726824"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C3_02_bin3 is 5.09264033264033"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C3_02_bin4 is 2.60033010136919"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C3_03_bin1 is 1.9846242175918"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C3_03_bin2 is 2.46669694085181"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C3_03_bin3 is 2.26068095402696"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C3_03_bin4 is 1.95356632210581"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C3_04_bin1 is 2.21509202481258"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C3_04_bin2 is 2.73634652613465"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C3_04_bin3 is 2.40807150198209"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C3_04_bin4 is 6.11764705882353"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C3_05_bin1 is 1.62598766841745"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C3_05_bin2 is 3.07014382803856"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C3_05_bin3 is 4.88888888888889"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C3_05_bin4 is 81.5"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C3_06_bin1 is 2.32435937111577"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C3_06_bin2 is 3.25406984496406"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C3_06_bin3 is 2.02941176470588"
    ## [1] "reads < cells for lib47_SARSr-wts 268C3_06_bin4 , un-normalized (ratio 0.1 )"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C183_01_bin1 is 2.73760890680042"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C183_01_bin2 is 2.06527210918463"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C183_01_bin3 is 4.90597416199232"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C183_01_bin4 is 3.61528688108266"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C183_02_bin1 is 2.1960082363765"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C183_02_bin2 is 1.00363422560902"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C183_02_bin3 is 7.23480702906881"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C183_02_bin4 is 6.25334037840733"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C183_03_bin1 is 6.15917054131235"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C183_03_bin2 is 1.53032789664643"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C183_03_bin3 is 5.50470383153946"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C183_03_bin4 is 5.96376811594203"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C183_04_bin1 is 1.28069259549291"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C183_04_bin2 is 1.98679301992539"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C183_04_bin3 is 16.780121444551"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C183_04_bin4 is 23.5625"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C183_05_bin1 is 2.34410096533073"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C183_05_bin2 is 2.37487753162243"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C183_05_bin3 is 2.43333333333333"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C183_05_bin4 is 3.53846153846154"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C183_06_bin1 is 2.32435937111577"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C183_06_bin2 is 3.25406984496406"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C183_06_bin3 is 2.02941176470588"
    ## [1] "reads < cells for lib47_SARSr-wts 268C183_06_bin4 , un-normalized (ratio 0.1 )"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C185_01_bin1 is 1.47608754390705"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C185_01_bin2 is 3.6651059955781"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C185_01_bin3 is 2.41701328685778"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C185_01_bin4 is 2.81545992739089"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C185_02_bin1 is 3.43276220145379"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C185_02_bin2 is 1.52693581253184"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C185_02_bin3 is 2.45255146253331"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C185_02_bin4 is 2.72592926408402"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C185_03_bin1 is 9.09196045513897"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C185_03_bin2 is 3.66981597687372"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C185_03_bin3 is 1.33537577002053"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C185_03_bin4 is 1.49095241554668"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C185_04_bin1 is 3.38683886838868"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C185_04_bin2 is 2.71560236566019"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C185_04_bin3 is 4.74693340431301"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C185_04_bin4 is 11.0441176470588"
    ## [1] "reads < cells for lib47_SARSr-wts 268C185_05_bin1 , un-normalized (ratio 0.663139228542216 )"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C185_05_bin2 is 3.11581033433785"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C185_05_bin3 is 1.56756756756757"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C185_05_bin4 is 2.06"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C185_06_bin1 is 2.32435937111577"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C185_06_bin2 is 3.25406984496406"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C185_06_bin3 is 2.02941176470588"
    ## [1] "reads < cells for lib47_SARSr-wts 268C185_06_bin4 , un-normalized (ratio 0.1 )"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C61_01_bin1 is 11.0870574754237"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C61_01_bin2 is 8.95213773061512"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C61_01_bin3 is 2.41343243513552"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C61_01_bin4 is 2.89317357913373"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C61_02_bin1 is 3.52419757192749"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C61_02_bin2 is 7.68793647756625"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C61_02_bin3 is 2.70205421350941"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C61_02_bin4 is 3.54797775799138"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C61_03_bin1 is 6.66068645663205"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C61_03_bin2 is 6.17342570018434"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C61_03_bin3 is 3.40263000662515"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C61_03_bin4 is 2.78388929816576"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C61_04_bin1 is 2.21340764246287"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C61_04_bin2 is 3.07686814610325"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C61_04_bin3 is 3.77789706668725"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C61_04_bin4 is 1.91304347826087"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C61_05_bin1 is 3.29502268755782"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C61_05_bin2 is 2.86783755659926"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C61_05_bin3 is 2.21818181818182"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C61_05_bin4 is 3.15384615384615"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C61_06_bin1 is 2.32435937111577"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C61_06_bin2 is 3.25406984496406"
    ## [1] "read:cell ratio for lib47_SARSr-wts 268C61_06_bin3 is 2.02941176470588"
    ## [1] "reads < cells for lib47_SARSr-wts 268C61_06_bin4 , un-normalized (ratio 0.1 )"

``` r
#annotate each barcode as to whether it's a homolog variant, wildtype, synonymous muts only, stop, nonsynonymous, >1 nonsynonymous mutations
dt[,variant_class:=as.character(NA)]
dt[n_codon_substitutions==0, variant_class := "wildtype"]
dt[n_codon_substitutions > 0 & n_aa_substitutions==0, variant_class := "synonymous"]
dt[n_aa_substitutions>0 & grepl("*",aa_substitutions,fixed=T), variant_class := "stop"]
dt[n_aa_substitutions == 1 & !grepl("*",aa_substitutions,fixed=T), variant_class := "1 nonsynonymous"]
dt[n_aa_substitutions > 1 & !grepl("*",aa_substitutions,fixed=T), variant_class := ">1 nonsynonymous"]

#cast the data frame into wide format
dt <- dcast(dt, library + barcode + target + variant_class + aa_substitutions + n_aa_substitutions ~ sample, value.var="count.norm")
```

## Calculating mean bin for each barcode at each sample concentration

Next, for each barcode at each of the sera concentrations, calculate the
“mean bin” response variable. This is calculated as a simple mean, where
the value of each bin is the integer value of the bin (bin1=unbound,
bin4=highly bound) – because of how bins are defined, the mean
fluorescence of cells in each bin are equally spaced on a log-normal
scale, so mean bin correlates with simple mean fluorescence.

We do not use the fluorescence boundaries of the FACS bins in our
calculations here, but we provide them for posterity’s sake below.

`(-288, 53), (53, 1099), (1099, 22649), (22649, 262143)`

``` r
#function that returns mean bin and sum of counts for four bins cell counts. Includes cutoffs for bimodal sample splits to filter out
calc.meanbin <- function(vec, split13filter=0.4, split24filter=0.4, split14filter=0.2){
  total <- sum(vec)
  if(is.na(total) | (vec[1] > split13filter*total & vec[3] > split13filter*total) | (vec[2] > split24filter*total & vec[4] > split24filter*total) | (vec[1] > split14filter*total & vec[4] > split14filter*total)){
    return(list(NA,NA))
  }else{
    return( list((vec[1]*1+vec[2]*2+vec[3]*3+vec[4]*4)/(vec[1]+vec[2]+vec[3]+vec[4]), total) )
  }
}
  

#iterate through titration samples, compute mean_bin and total_count for each barcode variant
#268C61
for(i in 1:nrow(samples_268C61)){ #iterate through titeseq sample (concentration)
  meanbin_out <- paste(samples_268C61[i,"sample"],"_meanbin",sep="") #define the header name for the meanbin output for the given concentration sample
  totalcount_out <- paste(samples_268C61[i,"sample"],"_totalcount",sep="") #define the header name for the total cell count output for the given concentration sample
  bin1_in <- paste(samples_268C61[i,"sample"],"_bin1",sep="") #define the header names for the input cell counts for bins1-4 of the given concnetration sample
  bin2_in <- paste(samples_268C61[i,"sample"],"_bin2",sep="")
  bin3_in <- paste(samples_268C61[i,"sample"],"_bin3",sep="")
  bin4_in <- paste(samples_268C61[i,"sample"],"_bin4",sep="")
  dt[,c(meanbin_out,totalcount_out) := calc.meanbin(c(get(bin1_in),get(bin2_in),get(bin3_in),get(bin4_in))),by=c("barcode","library")]
}

#268C183
for(i in 1:nrow(samples_268C183)){ #iterate through titeseq sample (concentration)
  meanbin_out <- paste(samples_268C183[i,"sample"],"_meanbin",sep="") #define the header name for the meanbin output for the given concentration sample
  totalcount_out <- paste(samples_268C183[i,"sample"],"_totalcount",sep="") #define the header name for the total cell count output for the given concentration sample
  bin1_in <- paste(samples_268C183[i,"sample"],"_bin1",sep="") #define the header names for the input cell counts for bins1-4 of the given concnetration sample
  bin2_in <- paste(samples_268C183[i,"sample"],"_bin2",sep="")
  bin3_in <- paste(samples_268C183[i,"sample"],"_bin3",sep="")
  bin4_in <- paste(samples_268C183[i,"sample"],"_bin4",sep="")
  dt[,c(meanbin_out,totalcount_out) := calc.meanbin(c(get(bin1_in),get(bin2_in),get(bin3_in),get(bin4_in))),by=c("barcode","library")]
}

#268C185
for(i in 1:nrow(samples_268C185)){ #iterate through titeseq sample (concentration)
  meanbin_out <- paste(samples_268C185[i,"sample"],"_meanbin",sep="") #define the header name for the meanbin output for the given concentration sample
  totalcount_out <- paste(samples_268C185[i,"sample"],"_totalcount",sep="") #define the header name for the total cell count output for the given concentration sample
  bin1_in <- paste(samples_268C185[i,"sample"],"_bin1",sep="") #define the header names for the input cell counts for bins1-4 of the given concnetration sample
  bin2_in <- paste(samples_268C185[i,"sample"],"_bin2",sep="")
  bin3_in <- paste(samples_268C185[i,"sample"],"_bin3",sep="")
  bin4_in <- paste(samples_268C185[i,"sample"],"_bin4",sep="")
  dt[,c(meanbin_out,totalcount_out) := calc.meanbin(c(get(bin1_in),get(bin2_in),get(bin3_in),get(bin4_in))),by=c("barcode","library")]
}

#268C3
for(i in 1:nrow(samples_268C3)){ #iterate through titeseq sample (concentration)
  meanbin_out <- paste(samples_268C3[i,"sample"],"_meanbin",sep="") #define the header name for the meanbin output for the given concentration sample
  totalcount_out <- paste(samples_268C3[i,"sample"],"_totalcount",sep="") #define the header name for the total cell count output for the given concentration sample
  bin1_in <- paste(samples_268C3[i,"sample"],"_bin1",sep="") #define the header names for the input cell counts for bins1-4 of the given concnetration sample
  bin2_in <- paste(samples_268C3[i,"sample"],"_bin2",sep="")
  bin3_in <- paste(samples_268C3[i,"sample"],"_bin3",sep="")
  bin4_in <- paste(samples_268C3[i,"sample"],"_bin4",sep="")
  dt[,c(meanbin_out,totalcount_out) := calc.meanbin(c(get(bin1_in),get(bin2_in),get(bin3_in),get(bin4_in))),by=c("barcode","library")]
}
```

## Fit binding curves

We will calculate a simple AUC metric across each barcode’s titration
series. We will also include a minimum cell count that is required for a
meanbin estimate to be used in the titration fit, and a minimum number
of concentrations with determined meanbin that is required for a
titration to be reported.

``` r
#For QC and filtering, output columns giving the average number of cells that were sampled for a barcode across the 9 sample concentrations, and a value for the number of meanbin estimates that were removed for being below the # of cells cutoff
cutoff <- 2

dt[,`268C61_avgcount` := mean(c(`268C61_01_totalcount`,`268C61_02_totalcount`,`268C61_03_totalcount`,
                             `268C61_04_totalcount`,`268C61_05_totalcount`,`268C61_06_totalcount`),na.rm=T),by=c("library","barcode")]

#number of concentrations at which meanbin is calculated from < cutoff cells or is missing b/c filtered for bimodality
dt[,`268C61_min_cell_filtered` := sum(c(c(`268C61_01_totalcount`,`268C61_02_totalcount`,`268C61_03_totalcount`,
                                       `268C61_04_totalcount`,`268C61_05_totalcount`,`268C61_06_totalcount`)<cutoff,
                                     is.na(c(`268C61_01_totalcount`,`268C61_02_totalcount`,`268C61_03_totalcount`,
                                             `268C61_04_totalcount`,`268C61_05_totalcount`,`268C61_06_totalcount`))),na.rm=T),by=c("library","barcode")]

dt[,`268C183_avgcount` := mean(c(`268C183_01_totalcount`,`268C183_02_totalcount`,`268C183_03_totalcount`,
                             `268C183_04_totalcount`,`268C183_05_totalcount`,`268C183_06_totalcount`),na.rm=T),by=c("library","barcode")]

#number of concentrations at which meanbin is calculated from < cutoff cells or is missing b/c filtered for bimodality
dt[,`268C183_min_cell_filtered` := sum(c(c(`268C183_01_totalcount`,`268C183_02_totalcount`,`268C183_03_totalcount`,
                                       `268C183_04_totalcount`,`268C183_05_totalcount`,`268C183_06_totalcount`)<cutoff,
                                     is.na(c(`268C183_01_totalcount`,`268C183_02_totalcount`,`268C183_03_totalcount`,
                                             `268C183_04_totalcount`,`268C183_05_totalcount`,`268C183_06_totalcount`))),na.rm=T),by=c("library","barcode")]

dt[,`268C185_avgcount` := mean(c(`268C185_01_totalcount`,`268C185_02_totalcount`,`268C185_03_totalcount`,
                             `268C185_04_totalcount`,`268C185_05_totalcount`,`268C185_06_totalcount`),na.rm=T),by=c("library","barcode")]

#number of concentrations at which meanbin is calculated from < cutoff cells or is missing b/c filtered for bimodality
dt[,`268C185_min_cell_filtered` := sum(c(c(`268C185_01_totalcount`,`268C185_02_totalcount`,`268C185_03_totalcount`,
                                       `268C185_04_totalcount`,`268C185_05_totalcount`,`268C185_06_totalcount`)<cutoff,
                                     is.na(c(`268C185_01_totalcount`,`268C185_02_totalcount`,`268C185_03_totalcount`,
                                             `268C185_04_totalcount`,`268C185_05_totalcount`,`268C185_06_totalcount`))),na.rm=T),by=c("library","barcode")]

dt[,`268C3_avgcount` := mean(c(`268C3_01_totalcount`,`268C3_02_totalcount`,`268C3_03_totalcount`,
                             `268C3_04_totalcount`,`268C3_05_totalcount`,`268C3_06_totalcount`),na.rm=T),by=c("library","barcode")]

#number of concentrations at which meanbin is calculated from < cutoff cells or is missing b/c filtered for bimodality
dt[,`268C3_min_cell_filtered` := sum(c(c(`268C3_01_totalcount`,`268C3_02_totalcount`,`268C3_03_totalcount`,
                                       `268C3_04_totalcount`,`268C3_05_totalcount`,`268C3_06_totalcount`)<cutoff,
                                     is.na(c(`268C3_01_totalcount`,`268C3_02_totalcount`,`268C3_03_totalcount`,
                                             `268C3_04_totalcount`,`268C3_05_totalcount`,`268C3_06_totalcount`))),na.rm=T),by=c("library","barcode")]

#function that fits a nls regression to the titration series, including an option to filter below certain thresholds for average cells across all samples, and number of samples below a cutoff of cells
fit.titration <- function(y.vals,x.vals,count.vals,min.cfu=cutoff,
                          min.means=1,min.average=7,EC50.start=10,
                          a.start=3,a.lower=2,a.upper=3,
                          b.start=1,b.lower=1,b.upper=1.5,
                          n.start=1,n.lower=1/2,n.upper=2){
  indices <- count.vals>min.cfu & !is.na(y.vals)
  y <- y.vals[indices]
  x <- x.vals[indices]
  if((length(y) < min.means*length(y.vals)) | (mean(count.vals,na.rm=T) < min.average)){ #return NAs if < min.means fraction of concentrations have above min.cfu counts or if the average count across all concentrations is below min.average
    return(list(as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA)))
  }else{
    fit <- nls(y ~ (a/(1+(EC50/x)^n))+b,
               start=list(a=a.start,b=b.start,n=n.start,EC50=EC50.start),
               lower=list(a=a.lower,b=b.lower,n=n.lower,EC50=min(x.vals[x.vals>0])/100), #constrain EC50 to be no lower than 1/100x the lowest concentration value
               upper=list(a=a.upper,b=b.upper,n=n.upper,EC50=max(x.vals[x.vals>0])*10), #constrain EC50 to be no higher than the 10x highest concentration value
               algorithm="port")
    y.pred <- predict(fit,newdata=list(x=x))
    resid <- y - y.pred
    resid.norm <- resid/as.numeric(summary(fit)$coefficients["a","Estimate"])
    nMSR <- mean((resid.norm)^2,na.rm=T)
    return(list(as.numeric(summary(fit)$coefficients["EC50","Estimate"]),
                as.numeric(summary(fit)$coefficients["EC50","Std. Error"]),
                as.numeric(summary(fit)$coefficients["a","Estimate"]),
                as.numeric(summary(fit)$coefficients["b","Estimate"]),
                as.numeric(nMSR)))
  }
}

#fit titration to 268C61 data for each barcode
dt[,c("EC50_268C61","EC50_SE_268C61","response_268C61","baseline_268C61","nMSR_268C61") :=
     tryCatch(fit.titration(y.vals=c(`268C61_01_meanbin`,`268C61_02_meanbin`,`268C61_03_meanbin`,
                                     `268C61_04_meanbin`,`268C61_05_meanbin`,`268C61_06_meanbin`),
                            x.vals=samples_268C61$conc,
                            count.vals=c(`268C61_01_totalcount`,`268C61_02_totalcount`,`268C61_03_totalcount`,
                                         `268C61_04_totalcount`,`268C61_05_totalcount`,`268C61_06_totalcount`)),
              error=function(e){list(as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA))}),by=c("library","barcode")]


#fit titration to 268C183 data for each barcode
dt[,c("EC50_268C183","EC50_SE_268C183","response_268C183","baseline_268C183","nMSR_268C183") :=
     tryCatch(fit.titration(y.vals=c(`268C183_01_meanbin`,`268C183_02_meanbin`,`268C183_03_meanbin`,
                                     `268C183_04_meanbin`,`268C183_05_meanbin`,`268C183_06_meanbin`),
                            x.vals=samples_268C183$conc,
                            count.vals=c(`268C183_01_totalcount`,`268C183_02_totalcount`,`268C183_03_totalcount`,
                                         `268C183_04_totalcount`,`268C183_05_totalcount`,`268C183_06_totalcount`)),
              error=function(e){list(as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA))}),by=c("library","barcode")]


#fit titration to 268C185 data for each barcode
dt[,c("EC50_268C185","EC50_SE_268C185","response_268C185","baseline_268C185","nMSR_268C185") :=
     tryCatch(fit.titration(y.vals=c(`268C185_01_meanbin`,`268C185_02_meanbin`,`268C185_03_meanbin`,
                                     `268C185_04_meanbin`,`268C185_05_meanbin`,`268C185_06_meanbin`),
                            x.vals=samples_268C185$conc,
                            count.vals=c(`268C185_01_totalcount`,`268C185_02_totalcount`,`268C185_03_totalcount`,
                                         `268C185_04_totalcount`,`268C185_05_totalcount`,`268C185_06_totalcount`)),
              error=function(e){list(as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA))}),by=c("library","barcode")]


#fit titration to 268C3 data for each barcode
dt[,c("EC50_268C3","EC50_SE_268C3","response_268C3","baseline_268C3","nMSR_268C3") :=
     tryCatch(fit.titration(y.vals=c(`268C3_01_meanbin`,`268C3_02_meanbin`,`268C3_03_meanbin`,
                                     `268C3_04_meanbin`,`268C3_05_meanbin`,`268C3_06_meanbin`),
                            x.vals=samples_268C3$conc,
                            count.vals=c(`268C3_01_totalcount`,`268C3_02_totalcount`,`268C3_03_totalcount`,
                                         `268C3_04_totalcount`,`268C3_05_totalcount`,`268C3_06_totalcount`)),
              error=function(e){list(as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA))}),by=c("library","barcode")]
```

## QC and sanity checks

We will do some QC to make sure we got good titration curves for most of
our library barcodes. We will also spot check titration curves from
across our measurement range, and spot check curves whose fit parameters
hit the different boundary conditions of the fit variables.

We successfully generated EC50 estimates for 0.5135152 of our 268C61 mAb
measurements, 0.7484381 of the 268C183 mAb, 0.5605636 of 268C185, and
0.7181563 of 268C3.

``` r
#make functions that allow me to plot a titration for any given row from the counts data frames, for spot checking curves
plot.titration <- function(row,mAb,output.text=F){
  y.vals <- c();for(sample in get(paste("samples_",mAb,sep=""))$sample){y.vals <- c(y.vals,paste(sample,"_meanbin",sep=""))};y.vals <- unlist(dt[row,y.vals,with=F])
  x.vals <- get(paste("samples_",mAb,sep=""))$conc
  count.vals <- c();for(sample in get(paste("samples_",mAb,sep=""))$sample){count.vals <- c(count.vals,paste(sample,"_totalcount",sep=""))};count.vals <- unlist(dt[row,count.vals,with=F])
  title <- paste(dt[row,target]," ",mAb,sep="")
  indices <- count.vals>cutoff & !is.na(count.vals)
  y.vals <- y.vals[indices]
  x.vals <- x.vals[indices]
  plot(x.vals,y.vals,xlab=paste("[",mAb,"] (ng/mL)",sep=""),
       ylab="mean bin",log="x",ylim=c(1,4),xlim=c(0.01,11000),pch=19,main=title)
  EC50_var <- paste("EC50_",mAb,sep="")
  fit <- nls(y.vals ~ (a/(1+(EC50/x.vals)^n))+b,
             start=list(a=3,b=1,EC50=dt[row,get(EC50_var)],n=1),
             lower=list(a=2,b=1,EC50=0.001,n=1/2),
             upper=list(a=3,b=1.25,EC50=100000,n=2),
             algorithm="port")
  if(!is.na(dt[row,get(EC50_var)])){
    lines(10^c(seq(-3,5,0.25)),predict(fit,newdata=list(x.vals=10^c(seq(-3,5,0.25)))))
    legend("topleft",bty="n",cex=1,legend=paste("EC50",format(dt[row,get(EC50_var)],digits=3),"ng/mL"))
  }
  if(output.text==T){ #for troubleshooting and interactive work, output some info from the counts table for the given row
    vars <- c("barcode","variant_class","wildtype","position","mutant",as.character(paste(mAb,"_avgcount",sep="")),as.character(paste(mAb,"_min_cell_filtered",sep="")),as.character(paste("EC50_",mAb,sep="")),as.character(paste("EC50_SE_",mAb,sep="")),as.character(paste("baseline_",mAb,sep="")),as.character(paste("response_",mAb,sep="")),as.character(paste("nMSR_",mAb,sep="")))
    return(dt[row,..vars])
  }
}
```

``` r
par(mfrow=c(3,3))
plot.titration(which(dt$EC50_268C61>=10000)[1],"268C61")
plot.titration(which(dt$EC50_268C61>=10000)[2],"268C61")
plot.titration(which(dt$EC50_268C61>=10000)[3],"268C61")
plot.titration(which(dt$EC50_268C61>1 & dt$EC50_268C61<10)[1],"268C61")
plot.titration(which(dt$EC50_268C61>1 & dt$EC50_268C61<10)[2],"268C61")
plot.titration(which(dt$EC50_268C61>1 & dt$EC50_268C61<10)[3],"268C61")
plot.titration(which(dt$EC50_268C61<0.5)[1],"268C61")
plot.titration(which(dt$EC50_268C61<0.5)[2],"268C61")
plot.titration(which(dt$EC50_268C61<0.5)[3],"268C61")
```

<img src="compute_EC50_files/figure-gfm/EC50_268C61-1.png" style="display: block; margin: auto;" />

``` r
par(mfrow=c(3,3))
plot.titration(which(dt$EC50_268C183>=10000)[1],"268C183")
plot.titration(which(dt$EC50_268C183>=10000)[2],"268C183")
plot.titration(which(dt$EC50_268C183>=10000)[3],"268C183")
plot.titration(which(dt$EC50_268C183>1 & dt$EC50_268C183<10)[1],"268C183")
plot.titration(which(dt$EC50_268C183>1 & dt$EC50_268C183<10)[2],"268C183")
plot.titration(which(dt$EC50_268C183>1 & dt$EC50_268C183<10)[3],"268C183")
plot.titration(which(dt$EC50_268C183<1)[1],"268C183")
plot.titration(which(dt$EC50_268C183<1)[2],"268C183")
plot.titration(which(dt$EC50_268C183<1)[3],"268C183")
```

<img src="compute_EC50_files/figure-gfm/EC50_268C183-1.png" style="display: block; margin: auto;" />

``` r
par(mfrow=c(3,3))
plot.titration(which(dt$EC50_268C185>=10000)[1],"268C185")
plot.titration(which(dt$EC50_268C185>=10000)[2],"268C185")
plot.titration(which(dt$EC50_268C185>=10000)[3],"268C185")
plot.titration(which(dt$EC50_268C185>1 & dt$EC50_268C185<10)[1],"268C185")
plot.titration(which(dt$EC50_268C185>1 & dt$EC50_268C185<10)[2],"268C185")
plot.titration(which(dt$EC50_268C185>1 & dt$EC50_268C185<10)[3],"268C185")
plot.titration(which(dt$EC50_268C185<0.5)[1],"268C185")
plot.titration(which(dt$EC50_268C185<0.5)[2],"268C185")
plot.titration(which(dt$EC50_268C185<0.5)[3],"268C185")
```

<img src="compute_EC50_files/figure-gfm/EC50_268C185-1.png" style="display: block; margin: auto;" />

``` r
par(mfrow=c(3,3))
plot.titration(which(dt$EC50_268C3>=10000)[1],"268C3")
plot.titration(which(dt$EC50_268C3>=10000)[2],"268C3")
plot.titration(which(dt$EC50_268C3>=10000)[3],"268C3")
plot.titration(which(dt$EC50_268C3>1 & dt$EC50_268C3<10)[1],"268C3")
plot.titration(which(dt$EC50_268C3>1 & dt$EC50_268C3<10)[2],"268C3")
plot.titration(which(dt$EC50_268C3>1 & dt$EC50_268C3<10)[3],"268C3")
plot.titration(which(dt$EC50_268C3<0.5)[1],"268C3")
plot.titration(which(dt$EC50_268C3<0.5)[2],"268C3")
plot.titration(which(dt$EC50_268C3<0.5)[3],"268C3")
```

<img src="compute_EC50_files/figure-gfm/EC50_268C3-1.png" style="display: block; margin: auto;" />

## Data filtering by fit quality

Next, let’s filter out poor fits using the value we previously computed,
the *normalized* mean square residual (nMSR). This metric computes the
residual between the observed response variable and that predicted from
the titration fit, normalizes this residual by the response range of the
titration fit (which is allowed to vary between 1.5 and 3 per the
titration fits above), and computes the mean-square of these normalized
residuals.

Distribution of the nMSR metric in each set of fits

``` r
par(mfrow=c(2,2))
hist(dt$nMSR_268C61,main="268C61",xlab="Response-normalized mean squared residual",col="gray50",breaks=40,xlim=c(0,0.6))
hist(dt$nMSR_268C183,main="268C183",xlab="Response-normalized mean squared residual",col="gray50",breaks=40,xlim=c(0,0.6))
hist(dt$nMSR_268C185,main="268C185",xlab="Response-normalized mean squared residual",col="gray50",breaks=40,xlim=c(0,0.6))
hist(dt$nMSR_268C3,main="268C3",xlab="Response-normalized mean squared residual",col="gray50",breaks=40,xlim=c(0,0.6))
```

<img src="compute_EC50_files/figure-gfm/nMSR_distribution-1.png" style="display: block; margin: auto;" />

As we would expect, the MSR stat decreases with cell count, indicating
that higher cell counts leads to better curve fits. Also show the cutoff
I’m proposing for nMSR (10x median across all fits), legend gives
percent of curve fits eliminated

``` r
median.nMSR <- median(c(dt$nMSR_268C61,dt$nMSR_268C183,dt$nMSR_268C185,dt$nMSR_268C3),na.rm=T)

par(mfrow=c(2,2))
plot(log10(dt$`268C61_avgcount`),dt$nMSR_268C61,main="268C61",pch=19,col="#00000010",xlab="average cell count (log10)",ylab="nMSR",xlim=c(1,3),ylim=c(0,0.6))
abline(h=10*median.nMSR,col="red",lty=2)
legend("topleft",bty="n",cex=1,legend=paste(format(100*nrow(dt[nMSR_268C61 > 10*median.nMSR & !is.na(nMSR_268C61),])/nrow(dt[!is.na(nMSR_268C61),]),digits=3),"%"))

plot(log10(dt$`268C183_avgcount`),dt$nMSR_268C183,main="268C183",pch=19,col="#00000010",xlab="average cell count (log10)",ylab="nMSR",xlim=c(1,3),ylim=c(0,0.6))
abline(h=10*median.nMSR,col="red",lty=2)
legend("topleft",bty="n",cex=1,legend=paste(format(100*nrow(dt[nMSR_268C183 > 10*median.nMSR & !is.na(nMSR_268C183),])/nrow(dt[!is.na(nMSR_268C183),]),digits=3),"%"))

plot(log10(dt$`268C185_avgcount`),dt$nMSR_268C185,main="268C185",pch=19,col="#00000010",xlab="average cell count (log10)",ylab="nMSR",xlim=c(1,3),ylim=c(0,0.6))
abline(h=10*median.nMSR,col="red",lty=2)
legend("topleft",bty="n",cex=1,legend=paste(format(100*nrow(dt[nMSR_268C185 > 10*median.nMSR & !is.na(nMSR_268C185),])/nrow(dt[!is.na(nMSR_268C185),]),digits=3),"%"))

plot(log10(dt$`268C3_avgcount`),dt$nMSR_268C3,main="268C3",pch=19,col="#00000010",xlab="average cell count (log10)",ylab="nMSR",xlim=c(1,3),ylim=c(0,0.6))
abline(h=10*median.nMSR,col="red",lty=2)
legend("topleft",bty="n",cex=1,legend=paste(format(100*nrow(dt[nMSR_268C3 > 10*median.nMSR & !is.na(nMSR_268C3),])/nrow(dt[!is.na(nMSR_268C3),]),digits=3),"%"))
```

<img src="compute_EC50_files/figure-gfm/nMSR_v_cell_count-1.png" style="display: block; margin: auto;" />

Next, we will apply this filtering step on normalized MSR, removing
curves with nMSR \>10x the median across all experiments

``` r
dt[nMSR_268C61 > 10*median.nMSR,c("EC50_268C61","EC50_SE_268C61","response_268C61","baseline_268C61") := list(as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA))]

dt[nMSR_268C183 > 10*median.nMSR,c("EC50_268C183","EC50_SE_268C183","response_268C183","baseline_268C183") := list(as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA))]

dt[nMSR_268C185 > 10*median.nMSR,c("EC50_268C185","EC50_SE_268C185","response_268C185","baseline_268C185") := list(as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA))]

dt[nMSR_268C3 > 10*median.nMSR,c("EC50_268C3","EC50_SE_268C3","response_268C3","baseline_268C3") := list(as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA))]
```

## Final scores

Let’s visualize the EC50 binding measurements as violin plots for the
different wildtype targets, for each serum metric.

``` r
p1 <- ggplot(dt[!is.na(EC50_268C61),],aes(x=target,y=EC50_268C61))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("268C61 EC50")+xlab("target")+theme(axis.text.x=element_text(angle=-90,hjust=0))+
  facet_wrap(~library,nrow=1)+
  scale_y_log10()

grid.arrange(p1,ncol=1)
```

<img src="compute_EC50_files/figure-gfm/binding_distribution_vioplot_268C61-1.png" style="display: block; margin: auto;" />

``` r
#save pdf
invisible(dev.print(pdf, paste(config$mAb_EC50_dir,"/violin-plot_EC50-by-target_268C61.pdf",sep="")))
```

``` r
p1 <- ggplot(dt[!is.na(EC50_268C183),],aes(x=target,y=EC50_268C183))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("268C183 EC50")+xlab("target")+theme(axis.text.x=element_text(angle=-90,hjust=0))+
  facet_wrap(~library,nrow=1)+
  scale_y_log10()

grid.arrange(p1,ncol=1)
```

<img src="compute_EC50_files/figure-gfm/binding_distribution_vioplot_268C183-1.png" style="display: block; margin: auto;" />

``` r
#save pdf
invisible(dev.print(pdf, paste(config$mAb_EC50_dir,"/violin-plot_EC50-by-target_268C183.pdf",sep="")))
```

``` r
p1 <- ggplot(dt[!is.na(EC50_268C185),],aes(x=target,y=EC50_268C185))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("268C185 EC50")+xlab("target")+theme(axis.text.x=element_text(angle=-90,hjust=0))+
  facet_wrap(~library,nrow=1)+
  scale_y_log10()

grid.arrange(p1,ncol=1)
```

<img src="compute_EC50_files/figure-gfm/binding_distribution_vioplot_268C185-1.png" style="display: block; margin: auto;" />

``` r
#save pdf
invisible(dev.print(pdf, paste(config$mAb_EC50_dir,"/violin-plot_EC50-by-target_268C185.pdf",sep="")))
```

``` r
p1 <- ggplot(dt[!is.na(EC50_268C3),],aes(x=target,y=EC50_268C3))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("268C3 EC50")+xlab("target")+theme(axis.text.x=element_text(angle=-90,hjust=0))+
  facet_wrap(~library,nrow=1)+
  scale_y_log10()

grid.arrange(p1,ncol=1)
```

<img src="compute_EC50_files/figure-gfm/binding_distribution_vioplot_268C3-1.png" style="display: block; margin: auto;" />

``` r
#save pdf
invisible(dev.print(pdf, paste(config$mAb_EC50_dir,"/violin-plot_EC50-by-target_268C3.pdf",sep="")))
```

## Save barcode-level metrics

In the next script, we will collapse bcs down to final
mutant/variant-level phenotypes, integrate things like expression
effects of variants, and visualize final phenotypes.

``` r
dt[,.(library,barcode,target,variant_class,
     `268C61_avgcount`,EC50_268C61,
     `268C183_avgcount`,EC50_268C183,
     `268C185_avgcount`,EC50_268C185,
     `268C3_avgcount`,EC50_268C3)] %>%
  mutate_if(is.numeric, round, digits=6) %>%
  write.csv(file=config$mAb_EC50_file, row.names=F)
```
