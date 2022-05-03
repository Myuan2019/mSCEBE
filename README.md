Instruction
-----------

mSCEBE is a simultaneous correction method that provides bias correction for multiple covariate amalysis based on empirical Bayesian estimation in mixed-effects models for longitudinal data.

#### The source codes and examples are available [**Here**](https://github.com/Myuan2019/mSCEBE).

Requirements
------------

mSCEBE requires the following R packages:

-  `Matrix`, `lme4`

You could download them directly in CRAN through the following commands in your R console.

    install.packages(c('Matrix', 'lme4'))
    
Usage
-----

-  `scebe_m.R` is the main program to perform the method.
-  `sampledata.RData` is a sample longitudinal data, and its corresponding covariate data matrix

### Input

-groupde data with the five columns representing the individual ID, the measuring time points, the trait, and the covariates, respectively.

    load('sampledata.RData')
    head(myData,10)
    
    ##    Grouped Data: Conc ~ Time | ID
    ##    <environment: 0x000000001a8962e8>
    ##   ID Time        Conc          x1          x2
    ##1   1 0.05  0.17478530  0.05319654 -0.32911900
    ##2   1 0.15  0.30276705  0.05319654 -0.32911900
    ##3   1 0.30  0.18952321  0.05319654 -0.32911900
    ##4   1 0.60  0.40129862  0.05319654 -0.32911900
    ##5   1 1.00  0.59090389  0.05319654 -0.32911900
    ##6   2 0.05  0.30323297 -0.09637448 -0.03060553
    ##7   2 0.15 -0.14195085 -0.09637448 -0.03060553
    ##8   2 0.30 -0.04760426 -0.09637448 -0.03060553
    ##9   2 0.60 -0.29212213 -0.09637448 -0.03060553
    ##10  2 1.00  0.04344139 -0.09637448 -0.03060553
    
-covariate data matrix with each row representing an individual, each column representing a covariate.
    
### Output

-The results of NEBE, mSCEBE, rmSCEBE are output as a list named as follows:

    rst <- scebe_m(mydata,x)
    names(rst)
    
    ## [1] "NEBE"    "mSCEBE"  "rmSCEBE"
    
-For each method, the covariate effect estimation, standard deviation, t-test statistics, and p-value are reported. Take rmSCEBE as an example:

    rst$rmSCEBE
    
           est         sd  t-value      p-value
    ##[1,] 0.4570546 0.06129157 7.457055 8.847754e-14
    ##[2,] 0.5925337 0.07996014 7.410364 1.259530e-13
    ##[3,] 0.5034304 0.06033973 8.343266 7.226640e-17
    ##[4,] 0.4595916 0.07906311 5.812972 6.137349e-09
    
    
    
    
