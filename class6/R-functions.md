R Functions
================
Charles Lam

Improving analysis code by writing functions
--------------------------------------------

A. We will improve the following code using functions.

``` r
# A. Can you improve this analysis code?
#df <- data.frame(a=1:10, b=seq(200,400,length=10),c=11:20,d=NA)
#df$a <- (df$a - min(df$a)) / (max(df$a) - min(df$a))
#df$b <- (df$b - min(df$a)) / (max(df$b) - min(df$b))
#df$c <- (df$c - min(df$c)) / (max(df$c) - min(df$c))
#df$d <- (df$d - min(df$d)) / (max(df$a) - min(df$d))
```

A. The simplified code is below

``` r
df <- data.frame(a=1:10, b=seq(200,400,length=10),c=11:20,d=NA)

# simplified function
df_simple <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

# simplified answers
df_simple(df$a)
```

    ##  [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667
    ##  [8] 0.7777778 0.8888889 1.0000000

``` r
df_simple(df$b)
```

    ##  [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667
    ##  [8] 0.7777778 0.8888889 1.0000000

``` r
df_simple(df$c)
```

    ##  [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667
    ##  [8] 0.7777778 0.8888889 1.0000000

``` r
df_simple(df$d)
```

    ##  [1] NA NA NA NA NA NA NA NA NA NA

B. Next improve the below example code for the analysis of protein drug interactions by abstracting the main activities in your own new function. Then answer questions 1 to 6 below.

``` r
# Can you improve this analysis code?
# s1 <- read.pdb("4AKE")  # kinase with drug
# s2 <- read.pdb("1AKE")  # kinase no drug
# s3 <- read.pdb("1E4Y")  # kinase with drug
# s1.chainA <- trim.pdb(s1, chain="A", elety="CA")
# s2.chainA <- trim.pdb(s2, chain="A", elety="CA")
# s3.chainA <- trim.pdb(s1, chain="A", elety="CA")
# s1.b <- s1.chainA$atom$b
# s2.b <- s2.chainA$atom$b
# s3.b <- s3.chainA$atom$b
# plotb3(s1.b, sse=s1.chainA, typ="l", ylab="Bfactor")
# plotb3(s2.b, sse=s2.chainA, typ="l", ylab="Bfactor")
# plotb3(s3.b, sse=s3.chainA, typ="l", ylab="Bfactor")
```

B. The simplified function is as follows.

``` r
#simplified function
protein_reader <- function(x) {
  read_file <- read.pdb(x)
  trim_read <- trim.pdb(read_file, chain="A", elety="CA")
  spec_prot <- trim_read$atom$b
  plot_prot <- plotb3(spec_prot, sse=trim_read, typ="l", ylab="Bfactor")
}

# simplified answers
protein_reader("4AKE") #kinase with drug
```

    ##   Note: Accessing on-line PDB file

![](R-functions_files/figure-markdown_github/unnamed-chunk-4-1.png)

``` r
protein_reader("1AKE") #kinase no drug
```

    ##   Note: Accessing on-line PDB file
    ##    PDB has ALT records, taking A only, rm.alt=TRUE

![](R-functions_files/figure-markdown_github/unnamed-chunk-4-2.png)

``` r
protein_reader("1E4Y") #kianse with drug
```

    ##   Note: Accessing on-line PDB file

![](R-functions_files/figure-markdown_github/unnamed-chunk-4-3.png)

Q1: What does the `read.pdb()` function do?
-------------------------------------------

The `read.pdb()` function reads the protein file taken from the Protein Data Bank (PDB).
