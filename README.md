
<!-- README.md is generated from README.Rmd. Please edit that file -->

:information\_source: `RegressionCalibration` is still under active
development.

## Overview

`RegressionCalibration` is an R-package to perform three-sample
Mendelian Randomisation (MR) analyses using the Regression Calibration
framework. MR estimates can be subject to bias due the use of weak
instruments and winner’s curse. The regression calibration framework can
be used when three samples (discovery and replication for the exposure,
and outcome) to account for these biases and re-calibrate the causal
effect estimate.  
This package builds up on the
[`TwoSampleMR`](https://github.com/MRCIEU/TwoSampleMR/) R-package for
the implementation of the most common MR estimators. We also recommend
using `TwoSampleMR` functions to pre-process the data (identify
independent instruments, harmonise data, etc…).

There are several functions available:  
- **`RC()`**  
main function, to perform regression calibration for all estimators
(IVW, simple/weighted median, simple/weighted mode)  
- **`RC_ivw()`**  
IVW implementation in the regression calibration framework  
- **`RC_simple_median()`**  
simple median implementation in the regression calibration framework  
- **`RC_weighted_median()`**  
weigthed median implementation in the regression calibration framework  
- **`RC_simple_mode()`**  
simple mode implementation in the regression calibration framework  
- **`RC_weighted_mode()`**  
weigthed mode implementation in the regression calibration framework

## Installation

You can install the current version of `RegressionCalibration` with:

``` r
# Directly install the package from github
# install.packages("remotes")
remotes::install_github("n-mounier/RegressionCalibration")
library(RegressionCalibration)
```

<!--- Note: using remotes instead of devtools leads to re-build the package
and apparently, it may be a problem with R 3.4 and macOS, 
see https://stackoverflow.com/questions/43595457/alternate-compiler-for-installing-r-packages-clang-error-unsupported-option/43943631#43943631 --->

## Usage - Main Function

The `RegressionCalibration` R-Package has been designed to be compatible
with the `TwoSampleMR` R-Package: you can use it to extract data and
pre-process the data.  
For example, we can use it to obtain BMI data from 3 different samples
(UKB, GIANT-males, GIANT-females):

``` r
# BMI (exposure - discovery) UKB Neale Lab 2017 
disc <- TwoSampleMR::extract_instruments(outcomes = 'ukb-a-248') 
```

    ## API: public: http://gwas-api.mrcieu.ac.uk/

``` r
# BMI (exposure - replication) GIANT 2015 males + BMI (outcome) GIANT 2015 females 
rep_out <- TwoSampleMR::extract_outcome_data(snps = disc$SNP, outcomes =c('ieu-a-785', 'ieu-a-974'))
```

    ## Extracting data for 315 SNP(s) from 2 GWAS(s)

    ## Finding proxies for 180 SNPs in outcome ieu-a-785

    ## Extracting data for 180 SNP(s) from 1 GWAS(s)

    ## Finding proxies for 180 SNPs in outcome ieu-a-974

    ## Extracting data for 180 SNP(s) from 1 GWAS(s)

``` r
# harmonise
dat <- TwoSampleMR::harmonise_data(disc, rep_out)
```

    ## Harmonising Body mass index (BMI) || id:ukb-a-248 (ukb-a-248) and Body mass index || id:ieu-a-785 (ieu-a-785)

    ## Removing the following SNPs for being palindromic with intermediate allele frequencies:
    ## rs11208779, rs1454687, rs1840661, rs2253310, rs7568228, rs815715, rs9489620, rs9536449

    ## Harmonising Body mass index (BMI) || id:ukb-a-248 (ukb-a-248) and Body mass index || id:ieu-a-974 (ieu-a-974)

    ## Removing the following SNPs for being palindromic with intermediate allele frequencies:
    ## rs11208779, rs1454687, rs1840661, rs2253310, rs7568228, rs815715, rs9489620, rs9536449

Note that when extracted data for the replication and for the outcome
samples using `TwoSampleMR::extract_outcome_data()`. To distinguish the
two, an extra column is needed:

``` r
# add column for regression calibration: sample should be "replication" or "outcome"
dat %>%
  mutate(sample = case_when(
    id.outcome == 'ieu-a-785' ~ 'replication',
    id.outcome == 'ieu-a-974' ~ 'outcome')) -> dat
```

The `dat` object now contains all information needed to perform the
Regression Calibration analyses:

``` r
# columns needed are: 
# SNP, exposure, outcome, beta.exposure, beta.outcome, se.exposure, se.outcome, mr_keep, sample
resRC = RegressionCalibration::RC(dat)
```

    ## Exposure, discovery:  Body mass index (BMI) || id:ukb-a-248 
    ## Exposure, replication:  Body mass index || id:ieu-a-785 
    ## Outcome:  Body mass index || id:ieu-a-974

``` r
resRC
```

    ##                      exposure.discovery            exposure.replication
    ## 1 Body mass index (BMI) || id:ukb-a-248 Body mass index || id:ieu-a-974
    ## 2 Body mass index (BMI) || id:ukb-a-248 Body mass index || id:ieu-a-974
    ## 3 Body mass index (BMI) || id:ukb-a-248 Body mass index || id:ieu-a-974
    ## 4 Body mass index (BMI) || id:ukb-a-248 Body mass index || id:ieu-a-974
    ## 5 Body mass index (BMI) || id:ukb-a-248 Body mass index || id:ieu-a-974
    ##                           outcome                    method nsnps          pval
    ## 1 Body mass index || id:ieu-a-974 Inverse variance weighted   239 9.385406e-201
    ## 2 Body mass index || id:ieu-a-974             Simple median   239 3.339987e-111
    ## 3 Body mass index || id:ieu-a-974           Weighted median   239 5.011951e-130
    ## 4 Body mass index || id:ieu-a-974               Simple mode   239  1.630656e-08
    ## 5 Body mass index || id:ieu-a-974             Weighted mode   239  3.514769e-52
    ##          b         se
    ## 1 1.045197 0.03457412
    ## 2 1.046779 0.04671594
    ## 3 1.053022 0.04340317
    ## 4 1.086436 0.19238431
    ## 5 1.055437 0.06943484

The output of the `RegressionCalibration::RC` is similar to the one of
`TwoSampleMR::mr()`.  
For this example, we can compare the results using the Regression
Calibration approach to the ones obtained without correction for
Winner’s Curse and weak instrument bias.

``` r
# just discovery -> outcome
resMR = TwoSampleMR::mr(dat %>% filter(sample == "outcome"))
```

    ## Analysing 'ukb-a-248' on 'ieu-a-974'

``` r
# using full GIANT sample as outcome
# also possible to use UKB on (full GIANT)
out <- TwoSampleMR::extract_outcome_data(snps = disc$SNP, outcomes =c('ieu-a-2'))
```

    ## Extracting data for 315 SNP(s) from 1 GWAS(s)

    ## Finding proxies for 177 SNPs in outcome ieu-a-2

    ## Extracting data for 177 SNP(s) from 1 GWAS(s)

``` r
dat_fullGIANT <- TwoSampleMR::harmonise_data(disc, out)
```

    ## Harmonising Body mass index (BMI) || id:ukb-a-248 (ukb-a-248) and Body mass index || id:ieu-a-2 (ieu-a-2)

    ## Removing the following SNPs for being palindromic with intermediate allele frequencies:
    ## rs11208779, rs1454687, rs1840661, rs2253310, rs7568228, rs815715, rs9489620, rs9536449, rs9835772

``` r
resMR_fullGIANT = TwoSampleMR::mr(dat_fullGIANT)
```

    ## Analysing 'ukb-a-248' on 'ieu-a-2'

``` r
# compare IVW estimates
resMR %>% filter(method=="Inverse variance weighted")
```

    ##   id.exposure id.outcome                         outcome
    ## 1   ukb-a-248  ieu-a-974 Body mass index || id:ieu-a-974
    ##                                exposure                    method nsnp
    ## 1 Body mass index (BMI) || id:ukb-a-248 Inverse variance weighted  239
    ##           b         se          pval
    ## 1 0.7873282 0.02604408 9.385406e-201

``` r
resMR_fullGIANT %>% filter(method=="Inverse variance weighted")
```

    ##   id.exposure id.outcome                       outcome
    ## 1   ukb-a-248    ieu-a-2 Body mass index || id:ieu-a-2
    ##                                exposure                    method nsnp
    ## 1 Body mass index (BMI) || id:ukb-a-248 Inverse variance weighted  240
    ##           b        se          pval
    ## 1 0.7650597 0.0231447 1.296585e-239

``` r
resRC %>% filter(method=="Inverse variance weighted")
```

    ##                      exposure.discovery            exposure.replication
    ## 1 Body mass index (BMI) || id:ukb-a-248 Body mass index || id:ieu-a-974
    ##                           outcome                    method nsnps          pval
    ## 1 Body mass index || id:ieu-a-974 Inverse variance weighted   239 9.385406e-201
    ##          b         se
    ## 1 1.045197 0.03457412

As the causal effect of BMI on itself is expected to be 1, we can see
that the standard IVW estimate (regardless of the outcome used) is
biased towards the null. The IVW using Regression Calibration approach
recovers the correct value.

## Usage - Other Functions

Other functions (method-specific) are also available. As they require a
different type of input (vectors of genetic effects and standard errors)
they should be use with care.

``` r
dat_disc = dat %>%
  filter(sample=="replication") %>%
  transmute(SNP,
            EA = effect_allele.exposure,
            OA = other_allele.exposure,
            beta = beta.exposure,
            se = se.exposure)
dat_rep = dat %>%
  filter(sample=="replication") %>%
  transmute(SNP,
            EA = effect_allele.outcome,
            OA = other_allele.outcome,
            beta = beta.outcome,
            se = se.outcome)
dat_out = dat %>%
  filter(sample=="outcome") %>%
  transmute(SNP,
            EA = effect_allele.outcome,
            OA = other_allele.outcome,
            beta = beta.outcome,
            se = se.outcome)
# check that for each sample-specific data we have the same SNPs in the same order
table(dat_disc$OA == dat_rep$OA)
```

    ## 
    ## TRUE 
    ##  247

``` r
table(dat_disc$EA == dat_rep$EA)
```

    ## 
    ## TRUE 
    ##  247

``` r
table(dat_disc$SNP == dat_out$SNP)
```

    ## 
    ## TRUE 
    ##  247

``` r
table(dat_disc$OA == dat_out$OA)
```

    ## 
    ## TRUE 
    ##  247

``` r
table(dat_disc$EA == dat_out$EA)
```

    ## 
    ## TRUE 
    ##  247

``` r
RegressionCalibration::RC_ivw(b_disc = dat_disc$beta,
                              b_rep = dat_rep$beta,
                              b_out = dat_out$beta,
                              se_disc = dat_disc$se,
                              se_rep = dat_rep$se,
                              se_out = dat_out$se)
```

    ## $b
    ## [1] 1.0463
    ## 
    ## $se
    ## [1] 0.03410832
    ## 
    ## $pval
    ## [1] 1.197067e-206
    ## 
    ## $nsnp
    ## [1] 247

## Citation

## Contact

<mounier.ninon@gmail.com>