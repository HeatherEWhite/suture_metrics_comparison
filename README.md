# A comparison of metrics for quantifying cranial suture complexity

__Authors:__
[Heather E White](mailto:heather.white.17@ucl.ac.uk), 
[Julien Clavel](https://github.com/JClavel),
Abigail S Tucker, 
Anjali Goswami

__To cite the paper__: 

>White HW, Clavel J, Abigail S Tucker, Goswami A. A comparison of complexity metrics for quantifying cranial suture complexity. Journal of the Royal Society Interface. 2020.

Available at: https://github.com/HeatherEWhite/suture_metrics_comparison

If using any of this code or data please cite the paper above and this repo

__To cite this repo__: 

> White HW, Clavel J, Abigail S Tucker, Goswami A. A comparison of complexity metrics for quantifying cranial suture complexity. Journal of the Royal Society Interface. 2020. Github: https://github.com/HeatherEWhite/suture_metrics_comparison plus the Zenodo DOI: 



## Data

In this folder you will find .csv files used for running code supplied in 'Code' folder

.csv documents have single tabs 

The data are provided in the `Raw data` folder
1. `CF_results` - the complexity factor results calculated for each specimen based on Saunders et al 1995 (.csv)
2. `PC_loadings` - PC loadings for each of the five methods on the PCA of complexity scores (.csv)
3. `SI_results` - sinuosity index (SI) complexity scores for each specimen, for later use in calculating SCI (.csv)
4. `STFT_results` - short-time Fourier transform array for PSD caluclation (.csv)
5. `Species_names_ordered` - specimen names for associating with landmarks and complexity method scores (.csv)
6. `all_complexity_results` - complexity scores for all five metrics for each specimen to be used in method comparison analysis (.csv)
7. `all_complexity_results_with_PC_scores` - complexity scores for all five metrics for each specimen to be used in method comparison analysis and PC scores from shape data to compare shape variation with complexity variation (.csv)
8. A sub-folder `Resampled_landmarks` containing 2D landmarks used for analysis

## Analysis
In this repository you will find raw data (.csv) and code for analyses (code supplied as .R files)

 :file_folder:
* **Data**

As above outlined above in 'Data' including a sub-folder 'Resampled_landmarks' containing 2D landmarks used for analysis.

 :file_folder:
* **Code_for_analyses**

`calculating_fractal_dimension.R`

`collecting_suture_curves.R`

`comparing_complexity_methods.R`

`opening_3D_skulls_in_R_and_taking_photos.R`

`PCA.R`

`SCI.R`

`sinuosity_index.R`

`stft_and_psd.R`



## License :page_with_curl:
This project is licensed under the MIT License - see the [LICENSE.md](https://github.com/HeatherEWhite/suture_complexity_metrics/blob/master/LICENSE) file for details

## Session Info :clipboard:
For reproducibility purposes, here is the output of `devtools::session_info()` used to perform the analyses in the publication. 

```{r}
─ Session info ───────────────────────────────────────────────────────────────────────────────
 setting  value                       
 version  R version 3.6.1 (2019-07-05)
 os       macOS Catalina 10.15.5      
 system   x86_64, darwin15.6.0        
 ui       RStudio                     
 language (EN)                        
 collate  en_GB.UTF-8                 
 ctype    en_GB.UTF-8                 
 tz       Europe/London               
 date     2020-08-03                  

─ Packages ───────────────────────────────────────────────────────────────────────────────────
 package           * version    date       lib source        
 ade4              * 1.7-13     2018-08-31 [1] CRAN (R 3.6.0)
 animation           2.6        2018-12-11 [1] CRAN (R 3.6.0)
 ape               * 5.3        2019-03-17 [1] CRAN (R 3.6.0)
 arrayhelpers      * 1.1-0      2020-02-04 [1] CRAN (R 3.6.0)
 assertthat          0.2.1      2019-03-21 [1] CRAN (R 3.6.0)
 backports           1.1.5      2019-10-02 [1] CRAN (R 3.6.0)
 bitops              1.0-6      2013-08-17 [1] CRAN (R 3.6.0)
 broom               0.5.3      2019-12-14 [1] CRAN (R 3.6.0)
 callr               3.4.3      2020-03-28 [1] CRAN (R 3.6.2)
 caTools             1.18.0     2020-01-17 [1] CRAN (R 3.6.0)
 cellranger          1.1.0      2016-07-27 [1] CRAN (R 3.6.0)
 cli                 2.0.2      2020-02-28 [1] CRAN (R 3.6.0)
 cluster             2.1.0      2019-06-19 [1] CRAN (R 3.6.1)
 clusterGeneration   1.3.4      2015-02-18 [1] CRAN (R 3.6.0)
 coda                0.19-3     2019-07-05 [1] CRAN (R 3.6.0)
 codetools           0.2-16     2018-12-24 [1] CRAN (R 3.6.1)
 colorRamps          2.3        2012-10-29 [1] CRAN (R 3.6.0)
 colorspace          1.4-1      2019-03-18 [1] CRAN (R 3.6.0)
 combinat            0.0-8      2012-10-29 [1] CRAN (R 3.6.0)
 corpcor           * 1.6.9      2017-04-01 [1] CRAN (R 3.6.0)
 corrgram          * 1.13       2018-07-09 [1] CRAN (R 3.6.0)
 corrplot          * 0.84       2017-10-16 [1] CRAN (R 3.6.0)
 crayon              1.3.4      2017-09-16 [1] CRAN (R 3.6.0)
 crosstalk           1.0.0      2016-12-21 [1] CRAN (R 3.6.0)
 data.table        * 1.12.8     2019-12-09 [1] CRAN (R 3.6.0)
 DBI                 1.1.0      2019-12-15 [1] CRAN (R 3.6.0)
 dbplyr              1.4.2      2019-06-17 [1] CRAN (R 3.6.0)
 dendextend          1.13.2     2019-12-02 [1] CRAN (R 3.6.0)
 desc                1.2.0      2018-05-01 [1] CRAN (R 3.6.0)
 devtools          * 2.3.1      2020-07-21 [1] CRAN (R 3.6.2)
 digest              0.6.22     2019-10-21 [1] CRAN (R 3.6.1)
 doParallel          1.0.15     2019-08-02 [1] CRAN (R 3.6.0)
 dotCall64           1.0-0      2018-07-30 [1] CRAN (R 3.6.0)
 dplyr             * 0.8.3      2019-07-04 [1] CRAN (R 3.6.0)
 ellipsis            0.3.0      2019-09-20 [1] CRAN (R 3.6.0)
 expm                0.999-4    2019-03-21 [1] CRAN (R 3.6.0)
 factoextra        * 1.0.6      2019-12-05 [1] CRAN (R 3.6.0)
 fansi               0.4.0      2018-10-05 [1] CRAN (R 3.6.0)
 fastmap             1.0.1      2019-10-08 [1] CRAN (R 3.6.0)
 fastmatch           1.1-0      2017-01-28 [1] CRAN (R 3.6.0)
 forcats           * 0.4.0      2019-02-17 [1] CRAN (R 3.6.0)
 foreach             1.4.7      2019-07-27 [1] CRAN (R 3.6.0)
 fs                  1.3.1      2019-05-06 [1] CRAN (R 3.6.0)
 gclus               1.3.2      2019-01-07 [1] CRAN (R 3.6.0)
 gdata               2.18.0     2017-06-06 [1] CRAN (R 3.6.0)
 generics            0.0.2      2018-11-29 [1] CRAN (R 3.6.0)
 geomorph          * 3.1.3      2019-10-09 [1] CRAN (R 3.6.0)
 ggplot2           * 3.2.1      2019-08-10 [1] CRAN (R 3.6.0)
 ggrepel             0.8.1      2019-05-07 [1] CRAN (R 3.6.0)
 glassoFast          1.0        2018-07-30 [1] CRAN (R 3.6.0)
 glue                1.3.1      2019-03-12 [1] CRAN (R 3.6.0)
 gplots              3.0.3      2020-02-25 [1] CRAN (R 3.6.0)
 gridExtra           2.3        2017-09-09 [1] CRAN (R 3.6.0)
 gtable              0.3.0      2019-03-25 [1] CRAN (R 3.6.0)
 gtools              3.8.1      2018-06-26 [1] CRAN (R 3.6.0)
 haven               2.2.0      2019-11-08 [1] CRAN (R 3.6.0)
 hms                 0.5.2      2019-10-30 [1] CRAN (R 3.6.0)
 htmltools           0.4.0      2019-10-04 [1] CRAN (R 3.6.0)
 htmlwidgets         1.5.1      2019-10-08 [1] CRAN (R 3.6.0)
 httpuv              1.5.2      2019-09-11 [1] CRAN (R 3.6.0)
 httr                1.4.1      2019-08-05 [1] CRAN (R 3.6.0)
 igraph              1.2.4.2    2019-11-27 [1] CRAN (R 3.6.0)
 iterators           1.0.12     2019-07-26 [1] CRAN (R 3.6.0)
 jpeg                0.1-8.1    2019-10-24 [1] CRAN (R 3.6.0)
 jsonlite            1.7.0      2020-06-25 [1] CRAN (R 3.6.2)
 KernSmooth          2.23-15    2015-06-29 [1] CRAN (R 3.6.1)
 knitr               1.26       2019-11-12 [1] CRAN (R 3.6.0)
 later               1.0.0      2019-10-04 [1] CRAN (R 3.6.0)
 lattice             0.20-38    2018-11-04 [1] CRAN (R 3.6.1)
 lazyeval            0.2.2      2019-03-15 [1] CRAN (R 3.6.0)
 lifecycle           0.1.0      2019-08-01 [1] CRAN (R 3.6.0)
 lubridate           1.7.4      2018-04-11 [1] CRAN (R 3.6.0)
 magrittr            1.5        2014-11-22 [1] CRAN (R 3.6.0)
 manipulateWidget    0.10.0     2018-06-11 [1] CRAN (R 3.6.0)
 maps              * 3.3.0      2018-04-03 [1] CRAN (R 3.6.0)
 MASS              * 7.3-51.4   2019-03-31 [1] CRAN (R 3.6.1)
 Matrix              1.2-17     2019-03-22 [1] CRAN (R 3.6.1)
 memoise             1.1.0      2017-04-21 [1] CRAN (R 3.6.0)
 mime                0.7        2019-06-11 [1] CRAN (R 3.6.0)
 miniUI              0.1.1.1    2018-05-18 [1] CRAN (R 3.6.0)
 mnormt              1.5-5      2016-10-15 [1] CRAN (R 3.6.0)
 modelr              0.1.5      2019-08-08 [1] CRAN (R 3.6.0)
 Morpho            * 2.7        2019-05-16 [1] CRAN (R 3.6.0)
 munsell             0.5.0      2018-06-12 [1] CRAN (R 3.6.0)
 mvMORPH           * 1.1.1      2020-01-08 [1] CRAN (R 3.6.0)
 nlme                3.1-140    2019-05-12 [1] CRAN (R 3.6.1)
 numDeriv            2016.8-1.1 2019-06-06 [1] CRAN (R 3.6.0)
 paleomorph        * 0.1.4      2017-04-19 [1] CRAN (R 3.6.0)
 pbmcapply           1.5.0      2019-07-10 [1] CRAN (R 3.6.0)
 phangorn            2.5.5      2019-06-19 [1] CRAN (R 3.6.0)
 phytools          * 0.7-47     2020-06-01 [1] CRAN (R 3.6.2)
 pillar              1.4.2      2019-06-29 [1] CRAN (R 3.6.0)
 pkgbuild            1.1.0      2020-07-13 [1] CRAN (R 3.6.2)
 pkgconfig           2.0.3      2019-09-22 [1] CRAN (R 3.6.0)
 pkgload             1.1.0      2020-05-29 [1] CRAN (R 3.6.2)
 plotrix             3.7-7      2019-12-05 [1] CRAN (R 3.6.0)
 prettyunits         1.0.2      2015-07-13 [1] CRAN (R 3.6.0)
 processx            3.4.1      2019-07-18 [1] CRAN (R 3.6.0)
 promises            1.1.0      2019-10-04 [1] CRAN (R 3.6.0)
 ps                  1.3.0      2018-12-21 [1] CRAN (R 3.6.0)
 purrr             * 0.3.3      2019-10-18 [1] CRAN (R 3.6.0)
 quadprog            1.5-8      2019-11-20 [1] CRAN (R 3.6.0)
 R6                  2.4.0      2019-02-14 [1] CRAN (R 3.6.0)
 Rcpp                1.0.2      2019-07-25 [1] CRAN (R 3.6.0)
 readr             * 1.3.1      2018-12-21 [1] CRAN (R 3.6.0)
 readxl              1.3.1      2019-03-13 [1] CRAN (R 3.6.0)
 registry            0.5-1      2019-03-05 [1] CRAN (R 3.6.0)
 remotes             2.2.0      2020-07-21 [1] CRAN (R 3.6.2)
 reprex              0.3.0      2019-05-16 [1] CRAN (R 3.6.0)
 rgl               * 0.100.30   2019-08-19 [1] CRAN (R 3.6.0)
 rlang               0.4.7      2020-07-09 [1] CRAN (R 3.6.2)
 rprojroot           1.3-2      2018-01-03 [1] CRAN (R 3.6.0)
 RRPP              * 0.4.3      2019-10-07 [1] CRAN (R 3.6.0)
 rstudioapi          0.11       2020-02-07 [1] CRAN (R 3.6.0)
 Rvcg                0.18       2018-09-28 [1] CRAN (R 3.6.0)
 rvest               0.3.5      2019-11-08 [1] CRAN (R 3.6.0)
 scales              1.0.0      2018-08-09 [1] CRAN (R 3.6.0)
 scatterplot3d       0.3-41     2018-03-14 [1] CRAN (R 3.6.0)
 seriation           1.2-8      2019-08-27 [1] CRAN (R 3.6.0)
 sessioninfo         1.1.1      2018-11-05 [1] CRAN (R 3.6.0)
 shiny               1.4.0      2019-10-10 [1] CRAN (R 3.6.0)
 spam                2.5-1      2019-12-12 [1] CRAN (R 3.6.0)
 stringi             1.4.3      2019-03-12 [1] CRAN (R 3.6.0)
 stringr           * 1.4.0      2019-02-10 [1] CRAN (R 3.6.0)
 subplex           * 1.5-4      2018-04-05 [1] CRAN (R 3.6.0)
 svd               * 0.5        2019-08-19 [1] CRAN (R 3.6.0)
 svUnit              1.0.3      2020-04-20 [1] CRAN (R 3.6.2)
 testthat            2.3.2      2020-03-02 [1] CRAN (R 3.6.0)
 tibble            * 2.1.3      2019-06-06 [1] CRAN (R 3.6.0)
 tidyr             * 1.0.0      2019-09-11 [1] CRAN (R 3.6.0)
 tidyselect          0.2.5      2018-10-11 [1] CRAN (R 3.6.0)
 tidyverse         * 1.3.0      2019-11-21 [1] CRAN (R 3.6.0)
 TSP                 1.1-10     2020-04-17 [1] CRAN (R 3.6.2)
 usethis           * 1.6.1      2020-04-29 [1] CRAN (R 3.6.2)
 vctrs               0.2.0      2019-07-05 [1] CRAN (R 3.6.0)
 viridis             0.5.1      2018-03-29 [1] CRAN (R 3.6.0)
 viridisLite         0.3.0      2018-02-01 [1] CRAN (R 3.6.0)
 webshot             0.5.2      2019-11-22 [1] CRAN (R 3.6.0)
 withr               2.1.2      2018-03-15 [1] CRAN (R 3.6.0)
 xfun                0.11       2019-11-12 [1] CRAN (R 3.6.0)
 xml2                1.2.2      2019-08-09 [1] CRAN (R 3.6.0)
 xtable              1.8-4      2019-04-21 [1] CRAN (R 3.6.0)
 zeallot             0.1.0      2018-01-28 [1] CRAN (R 3.6.0)

[1] /Library/Frameworks/R.framework/Versions/3.6/Resources/library
```
