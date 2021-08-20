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
library(IPO)
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
plotAdjustedRtime(xset)
```

![](LipidomicAnalysis_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

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

``` r
peakpickingParameters <- getDefaultXcmsSetStartingParams('centWave')
peakpickingParameters$min_peakwidth = c(0,10)
peakpickingParameters$max_peakwidth = c(10,30)
peakpickingParameters$ppm = c(0,10)
```

``` r
#resultPeakpicking <- optimizeXcmsSet(files = mzfiles, 
#                                     params = peakpickingParameters, 
#                                     nSlaves = 0, 
#                                     subdir = NULL)
```

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
    ##  [1] IPO_1.16.0          CAMERA_1.46.0       rsm_2.10.2         
    ##  [4] DT_0.18             ggplot2_3.3.5       xcms_3.12.0        
    ##  [7] MSnbase_2.16.1      ProtGenerics_1.22.0 S4Vectors_0.28.1   
    ## [10] mzR_2.24.1          Rcpp_1.0.7          BiocParallel_1.24.1
    ## [13] Biobase_2.50.0      BiocGenerics_0.36.1
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] bitops_1.0-7                matrixStats_0.60.0         
    ##  [3] doParallel_1.0.16           RColorBrewer_1.1-2         
    ##  [5] GenomeInfoDb_1.26.7         backports_1.2.1            
    ##  [7] tools_4.0.3                 utf8_1.2.2                 
    ##  [9] R6_2.5.0                    affyio_1.60.0              
    ## [11] rpart_4.1-15                Hmisc_4.5-0                
    ## [13] DBI_1.1.1                   colorspace_2.0-2           
    ## [15] nnet_7.3-16                 withr_2.4.2                
    ## [17] gridExtra_2.3               tidyselect_1.1.1           
    ## [19] compiler_4.0.3              MassSpecWavelet_1.56.0     
    ## [21] preprocessCore_1.52.1       graph_1.68.0               
    ## [23] htmlTable_2.2.1             DelayedArray_0.16.3        
    ## [25] checkmate_2.0.0             scales_1.1.1               
    ## [27] DEoptimR_1.0-9              robustbase_0.93-8          
    ## [29] affy_1.68.0                 RBGL_1.66.0                
    ## [31] stringr_1.4.0               digest_0.6.27              
    ## [33] foreign_0.8-81              rmarkdown_2.10             
    ## [35] XVector_0.30.0              jpeg_0.1-9                 
    ## [37] base64enc_0.1-3             pkgconfig_2.0.3            
    ## [39] htmltools_0.5.1.1           MatrixGenerics_1.2.1       
    ## [41] highr_0.9                   limma_3.46.0               
    ## [43] htmlwidgets_1.5.3           rlang_0.4.11               
    ## [45] rstudioapi_0.13             impute_1.64.0              
    ## [47] generics_0.1.0              mzID_1.28.0                
    ## [49] dplyr_1.0.7                 RCurl_1.98-1.3             
    ## [51] magrittr_2.0.1              GenomeInfoDbData_1.2.4     
    ## [53] Formula_1.2-4               MALDIquant_1.20            
    ## [55] Matrix_1.3-4                munsell_0.5.0              
    ## [57] fansi_0.5.0                 MsCoreUtils_1.2.0          
    ## [59] lifecycle_1.0.0             vsn_3.58.0                 
    ## [61] stringi_1.7.3               yaml_2.2.1                 
    ## [63] MASS_7.3-54                 SummarizedExperiment_1.20.0
    ## [65] zlibbioc_1.36.0             plyr_1.8.6                 
    ## [67] grid_4.0.3                  crayon_1.4.1               
    ## [69] lattice_0.20-44             splines_4.0.3              
    ## [71] knitr_1.33                  pillar_1.6.2               
    ## [73] igraph_1.2.6                GenomicRanges_1.42.0       
    ## [75] codetools_0.2-18            XML_3.99-0.6               
    ## [77] glue_1.4.2                  evaluate_0.14              
    ## [79] latticeExtra_0.6-29         data.table_1.14.0          
    ## [81] pcaMethods_1.82.0           BiocManager_1.30.16        
    ## [83] png_0.1-7                   vctrs_0.3.8                
    ## [85] foreach_1.5.1               gtable_0.3.0               
    ## [87] RANN_2.6.1                  purrr_0.3.4                
    ## [89] assertthat_0.2.1            xfun_0.25                  
    ## [91] ncdf4_1.17                  survival_3.2-12            
    ## [93] tibble_3.1.3                iterators_1.0.13           
    ## [95] IRanges_2.24.1              cluster_2.1.2              
    ## [97] ellipsis_0.3.2
