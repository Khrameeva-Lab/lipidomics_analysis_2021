The Hitchhiker’s Guide to untargeted lipidomics analysis: Practical
guidelines
================
D. Smirnov, P. Mazin, M. Osetrova, E. Stekolshchikova, E. Khrameeva
8/18/2021

## Introduction

## Package import

``` r
library(xcms)
library(ggplot2)
library(DT)
```

## Data import

``` r
#url <- "http://arcuda.skoltech.ru/~d.smirnov/sampledata.tar.gz"
#download.file(url, destfile = 'sampledata.tar.gz')
#untar('sampledata.tar.gz')
```

``` r
mzfiles <- list.files('sampledata/', recursive = TRUE, full.names = TRUE, pattern = '.mzML')
group <- unlist(lapply(strsplit(mzfiles,"/"), function (x) x[[3]]))
pd <- data.frame(sample_name = sub(basename(mzfiles), pattern = ".mzXML", replacement = "", fixed = TRUE), 
                 sample_group = group, 
                 stringsAsFactors = FALSE)
```

``` r
raw_data <- readMSData(files = mzfiles, 
                       pdata = new("NAnnotatedDataFrame", pd), 
                       mode = "onDisk", 
                       msLevel = 1, 
                       verbose = T, 
                       centroided = T)
```

    ## Reading 5229 spectra from file Biorec288729_UDN_POS.mzML

    ## Reading 5067 spectra from file Biorec288732_UDN_POS.mzML

## Peak picking

``` r
cwp <- CentWaveParam(peakwidth = c(9.5, 36),
                     ppm = 11.5,
                     noise = 0, 
                     snthresh = 10, 
                     mzdiff = -0.001, 
                     prefilter = c(3, 100), 
                     mzCenterFun = "wMean", 
                     integrate = 1, 
                     fitgauss = FALSE)
```

``` r
xset <- findChromPeaks(raw_data, param = cwp)
```

## Peak alignment

``` r
arp <- ObiwarpParam(distFun = "cor_opt", 
                    binSize = 1, 
                    response = 1, 
                    gapInit = 0.32, 
                    gapExtend = 2.688, 
                    factorDiag = 2, 
                    factorGap = 1,
                    localAlignment = FALSE)
```

``` r
xset <- adjustRtime(xset, param = arp)
```

## Peak grouping

``` r
pdp <- PeakDensityParam(sampleGroups = xset$sample_group, 
                        bw = 0.879999999999999, 
                        binSize = 0.02412, 
                        minFraction = 0.000001, 
                        minSamples = 1, 
                        maxFeatures = 50)
```

``` r
xset <- groupChromPeaks(xset, param = pdp)
```

    ## Processing 145180 mz slices ... OK

## Selection of parameters for peak picking, alignment, and grouping

## Imputation of missing values

``` r
xset <- fillChromPeaks(xset)
```

    ## Defining peak areas for filling-in .... OK
    ## Start integrating peak areas from original files

## Data export

``` r
pks <- chromPeaks(xset)
grs <- featureDefinitions(xset)
mtx <- featureValues(xset, method="maxint", value="into", filled=T) 
```

## Filtering of peaks

## Normalization

## Software used

``` r
sessionInfo()
```

    ## R version 4.0.3 (2020-10-10)
    ## Platform: x86_64-apple-darwin17.0 (64-bit)
    ## Running under: macOS Catalina 10.15.7
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRblas.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] ru_RU.UTF-8/ru_RU.UTF-8/ru_RU.UTF-8/C/ru_RU.UTF-8/ru_RU.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats4    parallel  stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] DT_0.18             ggplot2_3.3.5       xcms_3.12.0        
    ##  [4] MSnbase_2.16.1      ProtGenerics_1.22.0 S4Vectors_0.28.1   
    ##  [7] mzR_2.24.1          Rcpp_1.0.7          BiocParallel_1.24.1
    ## [10] Biobase_2.50.0      BiocGenerics_0.36.1
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] MatrixGenerics_1.2.1        vsn_3.58.0                 
    ##  [3] foreach_1.5.1               assertthat_0.2.1           
    ##  [5] BiocManager_1.30.16         affy_1.68.0                
    ##  [7] GenomeInfoDbData_1.2.4      yaml_2.2.1                 
    ##  [9] robustbase_0.93-8           impute_1.64.0              
    ## [11] pillar_1.6.2                lattice_0.20-44            
    ## [13] glue_1.4.2                  limma_3.46.0               
    ## [15] digest_0.6.27               GenomicRanges_1.42.0       
    ## [17] RColorBrewer_1.1-2          XVector_0.30.0             
    ## [19] colorspace_2.0-2            htmltools_0.5.1.1          
    ## [21] preprocessCore_1.52.1       Matrix_1.3-4               
    ## [23] plyr_1.8.6                  MALDIquant_1.20            
    ## [25] XML_3.99-0.6                pkgconfig_2.0.3            
    ## [27] zlibbioc_1.36.0             purrr_0.3.4                
    ## [29] scales_1.1.1                RANN_2.6.1                 
    ## [31] affyio_1.60.0               tibble_3.1.3               
    ## [33] generics_0.1.0              IRanges_2.24.1             
    ## [35] ellipsis_0.3.2              withr_2.4.2                
    ## [37] SummarizedExperiment_1.20.0 MassSpecWavelet_1.56.0     
    ## [39] magrittr_2.0.1              crayon_1.4.1               
    ## [41] evaluate_0.14               ncdf4_1.17                 
    ## [43] fansi_0.5.0                 doParallel_1.0.16          
    ## [45] MASS_7.3-54                 tools_4.0.3                
    ## [47] lifecycle_1.0.0             matrixStats_0.60.0         
    ## [49] stringr_1.4.0               munsell_0.5.0              
    ## [51] DelayedArray_0.16.3         pcaMethods_1.82.0          
    ## [53] compiler_4.0.3              GenomeInfoDb_1.26.7        
    ## [55] mzID_1.28.0                 rlang_0.4.11               
    ## [57] grid_4.0.3                  RCurl_1.98-1.3             
    ## [59] iterators_1.0.13            htmlwidgets_1.5.3          
    ## [61] MsCoreUtils_1.2.0           bitops_1.0-7               
    ## [63] rmarkdown_2.10              gtable_0.3.0               
    ## [65] codetools_0.2-18            DBI_1.1.1                  
    ## [67] R6_2.5.0                    knitr_1.33                 
    ## [69] dplyr_1.0.7                 utf8_1.2.2                 
    ## [71] stringi_1.7.3               vctrs_0.3.8                
    ## [73] DEoptimR_1.0-9              tidyselect_1.1.1           
    ## [75] xfun_0.25
