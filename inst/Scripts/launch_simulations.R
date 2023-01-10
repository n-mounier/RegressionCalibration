#############################################
####                                     ####
####      Winner's curse                 ####
####      Simulations                    ####
####                                     ####
####                                     ####
####      Ninon Mounier                  ####
####      10/01/2023                     ####
#############################################

# assume than in the current directory, there is a "results/" folder
# where results will be saved
working_dir = "results/"
setwd(working_dir)

# will create 7 RDS files
# balanced pleiotropy, 5e-8 (res_sim_balancedpleio.RDS)
# strong correlated pleiotropy, 5e-8 (res_sim_correlatedpleio_strong.RDS)
# weak correlated pleiotropy, 5e-8 (res_sim_correlatedpleio_weak.RDS)
# balanced pleiotropy, different thresholds (res_varyT_balancedpleio.RDS)
# strong correlated pleiotropy, different thresholds (res_varyT_correlatedpleio_strong.RDS)
# weak correlated pleiotropy, different thresholds (res_varyT_correlatedpleio_weak.RDS)
# balanced pleiotropy, different discovery/replication ratio (res_varyN_balancedpleio.RDS)
# balanced pleiotropy, different total exposure sample size (res_varyNtotexp_balancedpleio.RDS)
# balanced pleiotropy, different outcome sample size (res_varyNy_balancedpleio.RDS)

# each file contains a data.frame, with one row / simulated dataset
# the causal effect estimate for the different methods (IVW, IVW_3S, RAPS_3S, RAPS_UMVCUE, RAPS_ZP,
# IVW_RC, median_RC, Wmedian_RC, mode_RC, Wmode_RC)
# the F-statistic in each sample (discovery, replication, discovery+replication) (F_disc, F_rep, F_combined)
# the number of IVs in the discovery sample (m, null_IVs_disc, pleiotropic_IVs_disc)
# (total number of IVs, number of IVs violating the relevance assumption, number of IVs violating the exclusion restriction assumption
# and in the combined sample (m_combined, null_IVs_combined, pleiotropic_IVs_combined)

rm(list=ls())

library(mr.raps)
library(TwoSampleMR)
library(truncnorm)
library(gwas.winners.curse)

library(tidyverse)




#### main function to launch simulations ####
# to speed up things ignore median-/mode-based regression calibration
# when there is no pleiotropy
launch_simulations <- function(Nx_disc,   # exposure discovery sample size
                               Nx_rep,    # exposure replication sample size
                               Ny,        # outcome sample size
                               M,         # number of genetic variants
                               # assume spike-and-slab genetic architecture
                               # with normally distributed effect size
                               # - direct effects on the exposure
                               pi_x,      # proportion of causal genetic variants for the exposure
                               h2_x,      # (direct) heritability for the exposure
                               # - direct effects on the outcome
                               pi_y,      # proportion of causal genetic variants for the outcome
                               h2_y,      # (direct) heritability for the outcome
                               # - effects through a confounder (correlated pleiotropy)
                               # if no pleiotropy keep the following parameters = 0
                               pi_u = 0,  # proportion of causal genetic variants for the confounder
                               h2_u = 0,  # heritability for the confounder
                               q_x = 0,   # effect of the confounder on the exposure
                               q_y =0 ,   # effect of the confounder on the outcome
                               beta,      # causal effect (exposure on outcome)
                               threshold, # selection threshold (p-value)
                               sn, # number of simulations
                               noMedian_noMode = F){

  # note: we are assuming normalised effects, so phenotypes should be normalised
  # and var(X) = var(Y) = 1
  # make sure parameters are compatible with this
  if(h2_x+q_x^2>=1) stop("X can not have a variance of 1 if h2_x + q_x^2 >= 1.")
  if(h2_y+beta^2+q_y^2+2*beta*q_x*q_y>=1) stop("Y can not have a variance of 1 if h2_y + beta^2 + q_y^2 + 2 * beta * q_y * q_x >= 1.")


  ##### TRUE EFFECTS #####
  # simulate true effects for causal genetic variants
  # (causal SNPs for X / Y could be randomly overlapping)

  Mcausal_x = round(pi_x*M)
  effects_x = rnorm(Mcausal_x, 0, sqrt(h2_x/Mcausal_x))
  # ensure variance explained = h2, particularly needed if low number of causal SNPs
  effects_x = sqrt(h2_x)*effects_x/(sqrt(sum(effects_x^2)))
  effects_x = c(effects_x, rep(0, M - Mcausal_x))

  Mcausal_y = round(pi_y*M)
  effects_y = rnorm(Mcausal_y, 0, sqrt(h2_y/Mcausal_y))
  # ensure variance explained = h2, particularly needed if low number of causal SNPs
  effects_y = sqrt(h2_y)*effects_y/(sqrt(sum(effects_y^2)))
  effects_y = sample(c(effects_y, rep(0, M - Mcausal_y)),
                     size = M, replace = F)


  # if pleiotropy
  if(!any(c(pi_u, h2_u, q_x, q_y) == 0)){
    M_pleiotropy = round(pi_u*M)
    effects_pleiotropy = rnorm(M_pleiotropy, 0, sqrt(h2_u/M_pleiotropy))
    # could be randomly overlapping with causal SNPs for X / Y
    effects_pleiotropy = sample(c(effects_pleiotropy, rep(0, M - M_pleiotropy)),
                                size = M, replace = F)
    SNPs_pleiotropy = which(effects_pleiotropy!=0)
  } else {
    # force all of them to zero
    pi_u = h2_u = q_x = q_y = 0
    effects_pleiotropy = rep(0)
    SNPs_pleiotropy = c()
  }

  # results needed for figures/tables afterwards
  # (non-exhaustive list, we might need to add other things)
  # IVW - combined disc+rep
  # IVW_3samples -to keep
  # RAPS_3samples -to keep
  # RAPS_UMVCUE - MR-RAPS(UMVCUE)
  # RAPS_ZP - MR-RAPS(ZP)
  # IVW_RC
  # MedianBased_RC
  # WMedianBased_RC
  # ModeBased_RC
  # WModeBased_RC
  # F_disc, F_rep, F_combined
  # m, null_IVs_disc, pleiotropic_IVs_disc
  # m_combined, null_IVs_combined, pleiotropic_IVs_disc

  results = data.frame(IVW = rep(NA_real_, sn), IVW_3samples = rep(NA_real_, sn),
                       RAPS_3samples = rep(NA_real_, sn),
                       RAPS_UMVCUE = rep(NA_real_, sn), RAPS_ZP = rep(NA_real_, sn),
                       IVW_RC = rep(NA_real_, sn), MedianBased_RC = rep(NA_real_, sn), WMedianBased_RC = rep(NA_real_, sn),
                       ModeBased_RC = rep(NA_real_, sn), WModeBased_RC = rep(NA_real_, sn),
                       F_disc = rep(NA_real_, sn), F_rep =  rep(NA_real_, sn), F_combined = rep(NA_real_, sn),
                       m = rep(NA_real_, sn), null_IVs_disc = rep(NA_real_, sn), pleiotropic_IVs_disc = rep(NA_real_, sn),
                       m_combined = rep(NA_real_, sn), null_IVs_combined = rep(NA_real_, sn), pleiotropic_IVs_combined = rep(NA_real_, sn))

  if(noMedian_noMode){
    results$MedianBased_RC = NULL
    results$ModeBased_RC = NULL
    results$WMedianBased_RC = NULL
    results$WModeBased_RC = NULL
  }

  for(i in 1:sn){
    if(i%%50==0) print(i)

    # get genetic variant effect size estimates in each sample
    #  we are working with standardized effects
    # i.e. se = sqrt(1/N)
    # linear regression results X (discovery)
    betaX_disc = effects_x + effects_pleiotropy * q_x +
      rnorm(M, 0, sqrt(1/Nx_disc))
    betaX_disc_se = rep(1/sqrt(Nx_disc), M)
    # linear regression results X (replication)
    betaX_rep = effects_x +  effects_pleiotropy * q_x +
      rnorm(M, 0, sqrt(1/Nx_rep))
    betaX_rep_se = rep(1/sqrt(Nx_rep), M)
    # linear regression results Y
    betaY = effects_y + (effects_x + effects_pleiotropy * q_x) * beta +
      effects_pleiotropy * q_y + rnorm(M, 0, sqrt(1/Ny))
    betaY_se = rep(1/sqrt(Ny), M)

    betaX_combined = (betaX_disc/betaX_disc_se^2 + betaX_rep/betaX_rep_se^2) / (1/betaX_disc_se^2 + 1/betaX_rep_se^2)
    betaX_combined_se = sqrt(1 / (1/betaX_disc_se^2 + 1/betaX_rep_se^2))

    # select IVs based on p-value in discovery or in discovery+replication combined
    Tr = qnorm(threshold/2, lower.tail = F)

    IVs_disc = which(abs(betaX_disc/betaX_disc_se)>Tr)
    m = length(IVs_disc)
    null_IVs_disc = sum(IVs_disc>Mcausal_x & !IVs_disc %in% SNPs_pleiotropy)
    pleiotropic_IVs_disc = sum(IVs_disc %in% SNPs_pleiotropy)
    F_disc = mean(betaX_disc[IVs_disc]^2/betaX_disc_se[IVs_disc]^2)
    F_rep = mean(betaX_rep[IVs_disc]^2/betaX_rep_se[IVs_disc]^2)

    IVs_combined = which(abs(betaX_combined/betaX_combined_se)>Tr)
    m_combined = length(IVs_combined)
    null_IVs_combined = sum(IVs_combined>Mcausal_x & !IVs_combined %in% SNPs_pleiotropy)
    pleiotropic_IVs_combined = sum(IVs_combined %in% SNPs_pleiotropy)
    F_combined = mean(betaX_combined[IVs_combined]^2/betaX_combined_se[IVs_combined]^2)

    ### IVW : combined -> outcome
    IVW_model = TwoSampleMR::mr_ivw(betaX_combined[IVs_combined], betaY[IVs_combined], betaX_combined_se[IVs_combined], betaY_se[IVs_combined])
    IVW = IVW_model$b

    ### IVW : 3-sample design disc for selection / rep for estimation  -> outcome
    IVW_3S_model = TwoSampleMR::mr_ivw(betaX_rep[IVs_disc], betaY[IVs_disc], betaX_rep_se[IVs_disc], betaY_se[IVs_disc])
    IVW_3S = IVW_3S_model$b

    ### MR-RAPS: 3-sample design disc for selection / rep for estimation  -> outcome
    RAPS_3S = mr.raps.overdispersed.robust(betaX_rep[IVs_disc], betaY[IVs_disc], betaX_rep_se[IVs_disc], betaY_se[IVs_disc], loss.function = "huber")$beta.hat
    #RAPS_3S = mr.raps.overdispersed(betaX_rep[IVs_disc], betaY[IVs_disc], betaX_rep_se[IVs_disc], betaY_se[IVs_disc])$beta.hat


    ### UMVCUE, MR-RAPS
    get_RAPS_UMVCUE <- function(b1, s1, b2, s2, out, out_se, m, Tr){

      s1sq = s1^2
      s2sq = s2^2

      MLE = (s2sq*b1 + s1sq*b2)/(s1sq+s2sq)

      BB1 = (sqrt(s1sq+s2sq)/s1sq)*(MLE-Tr*s1)
      BB2 = (sqrt(s1sq+s2sq)/s1sq)*(MLE+Tr*s1)
      betaX_UMVCUE = MLE - (s2sq/sqrt(s1sq+s2sq))*((dnorm(BB1)-dnorm(BB2))/(pnorm(BB1)+pnorm(-BB2)))

      EstBoot = matrix(nrow=200,ncol=m)
      for(k in 1:200){

        K = m

        a = Tr*sqrt(s1sq)
        b = rep(Inf, K)

        minus = pnorm(-Tr-betaX_UMVCUE/sqrt(s1sq))/(pnorm(-Tr+betaX_UMVCUE/sqrt(s1sq)) + pnorm(-Tr-betaX_UMVCUE/sqrt(s1sq)))
        minus = (runif(K) <= minus)

        a[minus] = -Inf
        b[minus] = -Tr*sqrt(s1sq[minus])

        BXGs = rtruncnorm(1, a, b, mean = betaX_UMVCUE, sd = sqrt(s1sq))
        BXG2s = rnorm(K, mean = betaX_UMVCUE, sd = sqrt(s2sq))

        mle      = (s2sq*BXGs +  s1sq*BXG2s)/(s1sq+s2sq)
        BB1      = (sqrt(s1sq+s2sq)/s1sq)*(mle-Tr*sqrt(s1sq))
        BB2      = (sqrt(s1sq+s2sq)/s1sq)*(mle+Tr*sqrt(s1sq))
        EstBoot[k,] = mle - (s2sq/(sqrt(s1sq+s2sq)))*(dnorm(BB1)-dnorm(BB2))/(pnorm(BB1)+pnorm(-BB2))

      }

      betaX_UMVCUE_se = apply(EstBoot,2,sd)
      res = mr.raps.overdispersed.robust(betaX_UMVCUE, out, betaX_UMVCUE_se, out_se, loss.function = "huber")$beta.hat
      # res = mr.raps.overdispersed(betaX_UMVCUE, out, betaX_UMVCUE_se, out_se)$beta.hat

      return(res)
    }
    RAPS_UMVCUE =  get_RAPS_UMVCUE(b1 = betaX_disc[IVs_disc],  s1 = betaX_disc_se[IVs_disc],
                                   b2 = betaX_rep[IVs_disc], s2 = betaX_rep_se[IVs_disc],
                                   out = betaY[IVs_disc], out_se = betaY_se[IVs_disc], m, Tr)

    ### ZP, MR-RAPS (ZP on combined data)
    get_RAPS_ZP <- function(exp, exp_se, out, out_se, threshold, N){
      tmpfile = paste0("tmp_", runif(1, 1, 50000))
      input.file = paste0(tmpfile, ".tsv")
      output.file = paste0(tmpfile, "_out.tsv")
      dat = data.frame(beta_disc = exp,
                       se_disc = exp_se,
                       N_disc = N,
                       af_disc = 0.5, # because standardized effects
                       beta_rep = NA_real_,
                       se_rep = NA_real_,
                       N_rep = NA_real_,
                       af_rep = NA_real_)
      write.table(dat, input.file, sep="\t", quote=F, row.names=F)


      trait.mean <- 0 # because standardized effects
      p.threshold <- threshold
      gwas.winners.curse::correct.winners.curse(input.file,
                                                output.file,
                                                trait.mean,
                                                p.threshold,
                                                header = TRUE,
                                                sep = "\t")
      corrected = data.table::fread(output.file)

      res = mr.raps.overdispersed.robust(corrected$debiased.beta.mle, out, exp_se, out_se, loss.function = "huber")$beta.hat
      #res = mr.raps.overdispersed(corrected$debiased.beta.mle, out, exp_se, out_se)$beta.hat


      file.remove(input.file)
      file.remove(output.file)

      return(res)

    }

    RAPS_ZP =  get_RAPS_ZP(exp = betaX_combined[IVs_combined],  exp_se = betaX_combined_se[IVs_combined],
                           out = betaY[IVs_combined], out_se = betaY_se[IVs_combined], threshold, N = (Nx_disc+Nx_rep))



    ### IVW, Regression-Calibration
    get_IVW_RC <- function(b1, s1, b2, s2, out, out_se){
      # get IVW estimate (X discovery, Y)
      IVW_model = TwoSampleMR::mr_ivw(b1, out, s1, out_se)

      IVW_C = TwoSampleMR::mr_ivw(b1, b2, s1, s2)$b
      res = IVW_model$b/IVW_C

      return(res)


    }

    IVW_RC = get_IVW_RC(b1 = betaX_disc[IVs_disc],  s1 = betaX_disc_se[IVs_disc],
                        b2 = betaX_rep[IVs_disc], s2 = betaX_rep_se[IVs_disc],
                        out = betaY[IVs_disc], out_se = betaY_se[IVs_disc])

    if(noMedian_noMode){
      results[i, ] = c(IVW, IVW_3S, RAPS_3S,
                       RAPS_UMVCUE, RAPS_ZP,
                       IVW_RC,
                       F_disc, F_rep, F_combined,
                       m, null_IVs_disc, pleiotropic_IVs_disc,
                       m_combined, null_IVs_combined, pleiotropic_IVs_combined)
    } else {

      get_median_RC <- function(b1, s1, b2, s2, out, out_se, weighted=T){

        # median-based estimator
        if(weighted){
          MedianBased_model = TwoSampleMR::mr_weighted_median(b1, out, s1, out_se, parameters = list(nboot=2))

          MedianBased_C <- TwoSampleMR::mr_weighted_median(b1, b2, s1, s2, parameters = list(nboot=2))$b
        } else {
          MedianBased_model = TwoSampleMR::mr_simple_median(b1, out, s1, out_se, parameters = list(nboot=2))

          MedianBased_C <- TwoSampleMR::mr_simple_median(b1, b2, s1, s2, parameters = list(nboot=2))$b

        }

        res = MedianBased_model$b/MedianBased_C

        return(res)

      }

      median_RC = get_median_RC(b1 = betaX_disc[IVs_disc],  s1 = betaX_disc_se[IVs_disc],
                                b2 = betaX_rep[IVs_disc], s2 = betaX_rep_se[IVs_disc],
                                out = betaY[IVs_disc], out_se = betaY_se[IVs_disc], weighted = F)
      Wmedian_RC = get_median_RC(b1 = betaX_disc[IVs_disc],  s1 = betaX_disc_se[IVs_disc],
                                 b2 = betaX_rep[IVs_disc], s2 = betaX_rep_se[IVs_disc],
                                 out = betaY[IVs_disc], out_se = betaY_se[IVs_disc], weighted = T)


      get_mode_RC <- function(b1, s1, b2, s2, out, out_se, weighted=T){
        # mode-based estimator
        param = TwoSampleMR::default_parameters()
        param$nboot = 2
        if(weighted){
          ModeBased_model = TwoSampleMR::mr_weighted_mode(b1, out, s1, out_se, parameters = param)

          ModeBased_C <- TwoSampleMR::mr_weighted_mode(b1, b2, s1, s2, parameters = param)$b
        } else {
          ModeBased_model = TwoSampleMR::mr_simple_mode(b1, out, s1, out_se, parameters = param)

          ModeBased_C <- TwoSampleMR::mr_simple_mode(b1, b2, s1, s2, parameters = param)$b

        }

        res = ModeBased_model$b/ModeBased_C

        return(res)

      }

      mode_RC = get_mode_RC(b1 = betaX_disc[IVs_disc],  s1 = betaX_disc_se[IVs_disc],
                            b2 = betaX_rep[IVs_disc], s2 = betaX_rep_se[IVs_disc],
                            out = betaY[IVs_disc], out_se = betaY_se[IVs_disc], weighted = F)

      Wmode_RC = get_mode_RC(b1 = betaX_disc[IVs_disc],  s1 = betaX_disc_se[IVs_disc],
                             b2 = betaX_rep[IVs_disc], s2 = betaX_rep_se[IVs_disc],
                             out = betaY[IVs_disc], out_se = betaY_se[IVs_disc], weighted = T)

      results[i, ] = c(IVW, IVW_3S, RAPS_3S,
                       RAPS_UMVCUE, RAPS_ZP,
                       IVW_RC, median_RC, Wmedian_RC, mode_RC, Wmode_RC,
                       F_disc, F_rep, F_combined,
                       m, null_IVs_disc, pleiotropic_IVs_disc,
                       m_combined, null_IVs_combined, pleiotropic_IVs_combined)
    }
  }

  return(results)
}



#### tests with / without pleiotropy ####
### without pleiotropy ####
set.seed(0)
balanced_pleio = launch_simulations(Nx_disc = 90000,
                                    Nx_rep = 50000,
                                    Ny = 50000,
                                    M = 30000,
                                    pi_x = 0.1,
                                    h2_x = 0.15,
                                    pi_y = 0.05,
                                    h2_y = 0.1,
                                    pi_u = 0,
                                    h2_u = 0,
                                    q_x = 0,
                                    q_y = 0,
                                    beta = 0.2,
                                    threshold = 5e-8,
                                    sn=1000, noMedian_noMode = T)

saveRDS(balanced_pleio, "res_sim_balancedpleio.RDS")


### with pleiotropy - strong ####
set.seed(0)
correlated_pleio_strong = launch_simulations(Nx_disc = 90000,
                                             Nx_rep = 50000,
                                             Ny = 50000,
                                             M = 30000,
                                             pi_x = 0.10,
                                             h2_x = 0.15,
                                             pi_y = 0.05,
                                             h2_y = 0.1,
                                             pi_u = 0.0005,
                                             h2_u = 0.2,
                                             q_x = 0.4,
                                             q_y = 0.3,
                                             beta = 0.2,
                                             threshold = 5e-8,
                                             sn=1000)


saveRDS(correlated_pleio_strong, "res_sim_correlatedpleio_strong.RDS")


### with pleiotropy - weak ####
set.seed(0)
correlated_pleio_weak = launch_simulations(Nx_disc = 90000,
                                           Nx_rep = 50000,
                                           Ny = 50000,
                                           M = 30000,
                                           pi_x = 0.10,
                                           h2_x = 0.15,
                                           pi_y = 0.05,
                                           h2_y = 0.1,
                                           pi_u = 0.01,
                                           h2_u = 0.2,
                                           q_x = 0.4,
                                           q_y = 0.3,
                                           beta = 0.2,
                                           threshold = 5e-8,
                                           sn=1000)


saveRDS(correlated_pleio_weak, "res_sim_correlatedpleio_weak.RDS")


#### vary threshold ####
### without pleiotropy ####
## run simulations for a given design and multiple thresholds
# here we need larger sample size(s) to be able to get IVs for the most stringent thresholds
params = list(Nx_disc = 90000,
              Nx_rep = 50000,
              Ny = 50000,
              M = 30000,
              pi_x = 0.10,
              h2_x = 0.15,
              pi_y = 0.05,
              h2_y = 0.1,
              pi_u = 0,
              h2_u = 0,
              q_x = 0,
              q_y = 0,
              beta = 0.2,
              sn=1000,
              noMedian_noMode = T)
if(exists("all_res_balanced_pleio")) rm(all_res_balanced_pleio)
# can't go above 7, no IVs in some datasets
for(my_threshold in seq(3, 7, 0.25)){
  print(my_threshold)
  my_pthreshold =  pnorm(-abs(my_threshold))*2
  set.seed(0)
  tmp = with(params, launch_simulations(Nx_disc, Nx_rep, Ny, M,
                                        pi_x, h2_x,
                                        pi_y, h2_y,
                                        pi_u, h2_u, q_x, q_y,
                                        beta,
                                        threshold = my_pthreshold,
                                        sn, noMedian_noMode))
  tmp$threshold = my_threshold
  tmp$pthreshold = my_pthreshold
  if(exists("all_res_balanced_pleio")){
    all_res_balanced_pleio = rbind(all_res_balanced_pleio, tmp)
  } else {
    all_res_balanced_pleio = tmp
  }
}

saveRDS(all_res_balanced_pleio, "res_varyT_balancedpleio.RDS")



### with pleiotropy - strong ####
params = list(Nx_disc = 90000,
              Nx_rep = 50000,
              Ny = 50000,
              M = 30000,
              pi_x = 0.10,
              h2_x = 0.15,
              pi_y = 0.05,
              h2_y = 0.1,
              pi_u = 0.0005,
              h2_u = 0.2,
              q_x = 0.4,
              q_y = 0.3,
              beta = 0.2,
              sn=1000,
              noMedian_noMode = F)
if(exists("all_res_correlated_pleio_strong")) rm(all_res_correlated_pleio_strong)
for(my_threshold in seq(3, 7, 0.25)){
  print(my_threshold)
  my_pthreshold =  pnorm(-abs(my_threshold))*2
  set.seed(1)
  tmp = with(params, launch_simulations(Nx_disc, Nx_rep, Ny, M,
                                        pi_x, h2_x,
                                        pi_y, h2_y,
                                        pi_u, h2_u, q_x, q_y,
                                        beta,
                                        threshold = my_pthreshold,
                                        sn, noMedian_noMode))
  tmp$threshold = my_threshold
  tmp$pthreshold = my_pthreshold
  if(exists("all_res_correlated_pleio_strong")){
    all_res_correlated_pleio_strong = rbind(all_res_correlated_pleio_strong, tmp)
  } else {
    all_res_correlated_pleio_strong = tmp
  }
}

saveRDS(all_res_correlated_pleio_strong, "res_varyT_correlatedpleio_strong.RDS")


### with pleiotropy - weak ####
params = list(Nx_disc = 90000,
              Nx_rep = 50000,
              Ny = 50000,
              M = 30000,
              pi_x = 0.10,
              h2_x = 0.15,
              pi_y = 0.05,
              h2_y = 0.1,
              pi_u = 0.01,
              h2_u = 0.2,
              q_x = 0.4,
              q_y = 0.3,
              beta = 0.2,
              sn=1000,
              noMedian_noMode = F)
if(exists("all_res_correlated_pleio_weak")) rm(all_res_correlated_pleio_weak)
for(my_threshold in seq(3, 7, 0.25)){
  print(my_threshold)
  my_pthreshold =  pnorm(-abs(my_threshold))*2
  set.seed(1)
  tmp = with(params, launch_simulations(Nx_disc, Nx_rep, Ny, M,
                                        pi_x, h2_x,
                                        pi_y, h2_y,
                                        pi_u, h2_u, q_x, q_y,
                                        beta,
                                        threshold = my_pthreshold,
                                        sn, noMedian_noMode))
  tmp$threshold = my_threshold
  tmp$pthreshold = my_pthreshold
  if(exists("all_res_correlated_pleio_weak")){
    all_res_correlated_pleio_weak = rbind(all_res_correlated_pleio_weak, tmp)
  } else {
    all_res_correlated_pleio_weak = tmp
  }
}

saveRDS(all_res_correlated_pleio_weak, "res_varyT_correlatedpleio_weak.RDS")

#### vary Nx_disc/Nx_rep ratio ####
# Nx_disc + Nx_rep = 140000 (constant)
params = list(Ny = 50000,
              M = 30000,
              pi_x = 0.10,
              h2_x = 0.15,
              pi_y = 0.05,
              h2_y = 0.1,
              pi_u = 0,
              h2_u = 0,
              q_x = 0,
              q_y = 0,
              beta = 0.2,
              threshold = 5e-8,
              sn=1000,
              noMedian_noMode = T)
Nx_tot = 140000
# start at 0.35 otherwise if Nx_disc < 45,000 no IVs at 5e-8
if(exists("all_res_N")) rm(all_res_N)
for(prop in seq(0.35, 0.9, 0.05)){
  print(prop)
  set.seed(1)
  tmp = with(params, launch_simulations(Nx_disc = prop * Nx_tot, Nx_rep = (1-prop) * Nx_tot,
                                        Ny, M,
                                        pi_x, h2_x,
                                        pi_y, h2_y,
                                        pi_u, h2_u, q_x, q_y,
                                        beta, threshold,
                                        sn, noMedian_noMode))
  tmp$prop = prop
  if(exists("all_res_N")){
    all_res_N = rbind(all_res_N, tmp)
  } else {
    all_res_N = tmp
  }
}

saveRDS(all_res_N, "res_varyN_balancedpleio.RDS")


# # larger N
# params = list(Ny = 50000,
#               M = 30000,
#               pi_x = 0.10,
#               h2_x = 0.15,
#               pi_y = 0.05,
#               h2_y = 0.1,
#               pi_u = 0,
#               h2_u = 0,
#               q_x = 0,
#               q_y = 0,
#               beta = 0.2,
#               threshold = 5e-8,
#               sn=1000,
#               noMedian_noMode = T)
# Nx_tot = 300000
# # start at 0.35 otherwise if Nx_disc < 45,000 no IVs at 5e-8
# if(exists("all_res_N_larger")) rm(all_res_N_larger)
# for(prop in seq(0.35, 0.9, 0.05)){
#   print(prop)
#   set.seed(1)
#   tmp = with(params, launch_simulations(Nx_disc = prop * Nx_tot, Nx_rep = (1-prop) * Nx_tot,
#                                         Ny, M,
#                                         pi_x, h2_x,
#                                         pi_y, h2_y,
#                                         pi_u, h2_u, q_x, q_y,
#                                         beta, threshold,
#                                         sn, noMedian_noMode))
#   tmp$prop = prop
#   if(exists("all_res_N_larger")){
#     all_res_N_larger = rbind(all_res_N_larger, tmp)
#   } else {
#     all_res_N_larger = tmp
#   }
# }
#
#
# saveRDS(all_res_N_larger, "res_varyN_balancedpleio_largerN.RDS")



#### vary total exposure (Nx_disc+Nx_rep) sample size ####
params = list(Ny = 50000,
              M = 30000,
              pi_x = 0.10,
              h2_x = 0.15,
              pi_y = 0.05,
              h2_y = 0.1,
              pi_u = 0,
              h2_u = 0,
              q_x = 0,
              q_y = 0,
              beta = 0.2,
              threshold = 5e-8,
              sn=1000,
              noMedian_noMode = T)
if(exists("all_res_Ntot")) rm(all_res_Ntot)
for(tot_N in seq(80000, 300000, 20000)){
  print(tot_N)
  set.seed(1)
  tmp = with(params, launch_simulations(Nx_disc = tot_N * 0.64, Nx_rep = tot_N * 0.36,
                                        Ny, M,
                                        pi_x, h2_x,
                                        pi_y, h2_y,
                                        pi_u, h2_u, q_x, q_y,
                                        beta, threshold,
                                        sn, noMedian_noMode))
  tmp$tot_N = tot_N
  if(exists("all_res_Ntot")){
    all_res_Ntot = rbind(all_res_Ntot, tmp)
  } else {
    all_res_Ntot = tmp
  }
}

saveRDS(all_res_Ntot, "res_varyNtotexp_balancedpleio.RDS")



#### vary outcome sample size ####
params = list(Nx_disc = 90000,
              Nx_rep = 50000,
              M = 30000,
              pi_x = 0.10,
              h2_x = 0.15,
              pi_y = 0.05,
              h2_y = 0.1,
              pi_u = 0,
              h2_u = 0,
              q_x = 0,
              q_y = 0,
              beta = 0.2,
              threshold = 5e-8,
              sn=1000,
              noMedian_noMode = T)
if(exists("all_res_Ny")) rm(all_res_Ny)
for(my_Ny in  seq(20000, 260000, 30000)){
  print(my_Ny)
  set.seed(1)
  tmp = with(params, launch_simulations(Nx_disc, Nx_rep,
                                        Ny=my_Ny, M,
                                        pi_x, h2_x,
                                        pi_y, h2_y,
                                        pi_u, h2_u, q_x, q_y,
                                        beta, threshold,
                                        sn, noMedian_noMode))
  tmp$Ny = my_Ny
  if(exists("all_res_Ny")){
    all_res_Ny = rbind(all_res_Ny, tmp)
  } else {
    all_res_Ny = tmp
  }
}
saveRDS(all_res_Ny, "res_varyNy_balancedpleio.RDS")
