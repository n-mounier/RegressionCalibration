
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

We also implemented the MR-RAPS using ZP and MR-RAPS using UMVCUE
approaches:  
- **`raps_ZP()`**  
- **`raps_UMVCUE`()\`**

More details can be found in the
[manual](doc/RegressionCalibration-manual.pdf).

Note that scripts to reproduce simulation results from our paper are
available here:  
[launch\_simulations.R](inst/Scripts/launch_simulations.R)  
[make\_figures\_simulations.R](inst/Scripts/launch_simulations.R)

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
print(resRC)
```

    ##                      exposure.discovery            exposure.replication
    ## 1 Body mass index (BMI) || id:ukb-a-248 Body mass index || id:ieu-a-974
    ## 2 Body mass index (BMI) || id:ukb-a-248 Body mass index || id:ieu-a-974
    ## 3 Body mass index (BMI) || id:ukb-a-248 Body mass index || id:ieu-a-974
    ## 4 Body mass index (BMI) || id:ukb-a-248 Body mass index || id:ieu-a-974
    ## 5 Body mass index (BMI) || id:ukb-a-248 Body mass index || id:ieu-a-974
    ##                           outcome                    method nsnps          pval
    ## 1 Body mass index || id:ieu-a-974 Inverse variance weighted   239 1.030811e-201
    ## 2 Body mass index || id:ieu-a-974             Simple median   239 1.527915e-125
    ## 3 Body mass index || id:ieu-a-974           Weighted median   239 1.768643e-150
    ## 4 Body mass index || id:ieu-a-974               Simple mode   239  3.278644e-09
    ## 5 Body mass index || id:ieu-a-974             Weighted mode   239  6.736170e-52
    ##          b         se
    ## 1 1.047358 0.03456228
    ## 2 1.080210 0.04532456
    ## 3 1.074307 0.04111761
    ## 4 1.089359 0.18410666
    ## 5 1.053989 0.06953475

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
# compare IVW estimates
resMR %>% filter(method=="Inverse variance weighted")
```

    ##   id.exposure id.outcome                         outcome
    ## 1   ukb-a-248  ieu-a-974 Body mass index || id:ieu-a-974
    ##                                exposure                    method nsnp
    ## 1 Body mass index (BMI) || id:ukb-a-248 Inverse variance weighted  239
    ##           b         se          pval
    ## 1 0.7889564 0.02603515 1.030811e-201

``` r
resRC %>% filter(method=="Inverse variance weighted")
```

    ##                      exposure.discovery            exposure.replication
    ## 1 Body mass index (BMI) || id:ukb-a-248 Body mass index || id:ieu-a-974
    ##                           outcome                    method nsnps          pval
    ## 1 Body mass index || id:ieu-a-974 Inverse variance weighted   239 1.030811e-201
    ##          b         se
    ## 1 1.047358 0.03456228

The causal effect of BMI on itself is expected to be 1, but we can see
that the standard IVW estimate is biased towards the null (0.789,
SE=0.026). The IVW using Regression Calibration approach recovers the
correct value (1.047, SE=0.035).

## Usage - Other Functions

Other functions (method-specific) are also available. As they require a
different type of input (vectors of genetic effects and standard errors)
they should be use with care.

``` r
dat_disc = dat %>%
  filter(sample=="replication", mr_keep) %>%
  transmute(SNP,
            EA = effect_allele.exposure,
            OA = other_allele.exposure,
            beta = beta.exposure,
            se = se.exposure)
dat_rep = dat %>%
  filter(sample=="replication", mr_keep) %>%
  transmute(SNP,
            EA = effect_allele.outcome,
            OA = other_allele.outcome,
            beta = beta.outcome,
            se = se.outcome)
dat_out = dat %>%
  filter(sample=="outcome", mr_keep) %>%
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
    ##  239

``` r
table(dat_disc$EA == dat_rep$EA)
```

    ## 
    ## TRUE 
    ##  239

``` r
table(dat_disc$SNP == dat_out$SNP)
```

    ## 
    ## TRUE 
    ##  239

``` r
table(dat_disc$OA == dat_out$OA)
```

    ## 
    ## TRUE 
    ##  239

``` r
table(dat_disc$EA == dat_out$EA)
```

    ## 
    ## TRUE 
    ##  239

``` r
RegressionCalibration::RC_ivw(b_disc = dat_disc$beta,
                              b_rep = dat_rep$beta,
                              b_out = dat_out$beta,
                              se_disc = dat_disc$se,
                              se_rep = dat_rep$se,
                              se_out = dat_out$se)
```

    ## $b
    ## [1] 1.047358
    ## 
    ## $se
    ## [1] 0.03456228
    ## 
    ## $pval
    ## [1] 1.030811e-201
    ## 
    ## $nsnp
    ## [1] 239

``` r
# MR-RAPS using ZP correction (as presented in the paper) 
# just discovery -> outcome
# less biased than standard IVW, but still not recovering true causal effect
RegressionCalibration::raps_ZP(b_exp = dat_disc$beta,
                               b_out = dat_out$beta,
                               se_exp = dat_disc$se,
                               se_out = dat_out$se,
                               threshold = 5e-8)
```

    ## $b
    ## [1] 0.8659352
    ## 
    ## $se
    ## [1] 0.02807875
    ## 
    ## $pval
    ## [1] 7.743111e-209
    ## 
    ## $nsnp
    ## [1] 239

``` r
# MR-RAPS using UMVCUE correction (as presented in the paper)
# less biased than standard IVW, but still not recovering true causal effect
RegressionCalibration::raps_UMVCUE(b_disc = dat_disc$beta,
                               b_rep = dat_rep$beta,
                               b_out = dat_out$beta,
                               se_disc = dat_disc$se,
                               se_rep = dat_rep$se,
                               se_out = dat_out$se,
                               threshold = 5e-8)
```

    ## $b
    ## [1] 0.9289908
    ## 
    ## $se
    ## [1] 0.0229748
    ## 
    ## $pval
    ## [1] 0
    ## 
    ## $nsnp
    ## [1] 239

## Citation

## Contact

<mounier.ninon@gmail.com>
