Class 14: Genome Informatics I
================

Asthma SNPs
-----------

Examine Asthma SNPs in the MXL (Mexican Ancestry in Los Angeles, California) 1000 Genomes sequencing data.

``` r
mxl <- read.csv("373531-SampleGenotypes-Homo_sapiens_Variation_Sample_rs8067378.csv")
```

Lets focus on 2nd column that contains genotype info

``` r
genotypes <- round(table( mxl[,2] ) / nrow(mxl) * 100, 2)
genotypes
```

    ## 
    ##   A|A   A|G   G|A   G|G 
    ## 34.38 32.81 18.75 14.06

There are 34.38 % AA genotype in this population.

Interpreting Base Qualities in R
--------------------------------

``` r
#install.packages("seqinr")
#install.packages("gtools")
```

``` r
library(seqinr)
library(gtools)
phred <- asc( s2c("DDDDCDEDCDDDDBBDDDCC@") ) - 33
phred
```

    ##  D  D  D  D  C  D  E  D  C  D  D  D  D  B  B  D  D  D  C  C  @ 
    ## 35 35 35 35 34 35 36 35 34 35 35 35 35 33 33 35 35 35 34 34 31

Section 4: Population Scale Analysis
------------------------------------

One sample is obviously not enough to know what is happening in a population. You are interested in assessing genetic differences on a population scale. So, you processed about ~230 samples and did the normalization on a genome level. Now, you want to find whether there is any association of the 4 asthma-associated SNPs (rs8067378...) on ORMDL3 expression.

This is the final file you got ( <https://bioboot.github.io/bimm143_S18/class-material/> rs8067378\_ENSG00000172057.6.txt ). The first column is sample name, the second column is genotype and the third column is the expression value.

``` r
## read expresion data file
expr <- read.table("rs8067378_ENSG00000172057.6.txt")
```

``` r
table(expr$geno)
```

    ## 
    ## A/A A/G G/G 
    ## 108 233 121

``` r
inds.aa <- expr$geno == "A/A"
summary(expr$exp[inds.aa])
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##   11.40   27.02   31.25   31.82   35.92   51.52

``` r
inds.ag <- expr$geno == "A/G"
summary(expr$exp[inds.ag])
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##   7.075  20.626  25.065  25.397  30.552  48.034

``` r
inds.gg <- expr$geno == "G/G"
summary(expr$exp[inds.gg])
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##   6.675  16.903  20.074  20.594  24.457  33.956

``` r
#boxplot(count ~ spray, data = InsectSprays, col = "lightgray")
boxplot(exp ~ geno, data=expr)
```

![](class14_files/figure-markdown_github/unnamed-chunk-10-1.png)
