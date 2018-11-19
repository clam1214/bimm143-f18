Untitled
================

PDB database comp stats
-----------------------

Download PDB stats data as CSV file and

``` r
pdbstats <- read.csv("~/Desktop/bimm143_fall18-master/Data Export Summary.csv")
```

Data Table
----------

Lets look at the table

``` r
library(knitr)
kable(pdbstats)
```

| Experimental.Method |  Proteins|  Nucleic.Acids|  Protein.NA.Complex|  Other|   Total|
|:--------------------|---------:|--------------:|-------------------:|------:|-------:|
| X-Ray               |    122263|           1960|                6333|     10|  130566|
| NMR                 |     10898|           1263|                 253|      8|   12422|
| Electron Microscopy |      1822|             31|                 657|      0|    2510|
| Other               |       244|              4|                   6|     13|     267|
| Multi Method        |       119|              5|                   2|      1|     127|

Q1
--

``` r
# total number of entries
nstru <- sum(pdbstats$Total)

precent <- round((pdbstats$Total / nstru)*100, 2)
```

A1: There are 89.49 % X-ray structures and 8.51 % EM structures solved.

``` r
nstats <- pdbstats
nstats$Precent <- precent
kable(nstats)
```

| Experimental.Method |  Proteins|  Nucleic.Acids|  Protein.NA.Complex|  Other|   Total|  Precent|
|:--------------------|---------:|--------------:|-------------------:|------:|-------:|--------:|
| X-Ray               |    122263|           1960|                6333|     10|  130566|    89.49|
| NMR                 |     10898|           1263|                 253|      8|   12422|     8.51|
| Electron Microscopy |      1822|             31|                 657|      0|    2510|     1.72|
| Other               |       244|              4|                   6|     13|     267|     0.18|
| Multi Method        |       119|              5|                   2|      1|     127|     0.09|

``` r
round((sum(nstats$Proteins) / nstru)*100, 2)
```

    ## [1] 92.77

``` r
library(bio3d)
```
