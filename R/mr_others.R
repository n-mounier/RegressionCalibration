#' MR-RAPS using ZP correction
#'
#' @param b_exp SNP-exposure association estimate
#' @param b_out SNP-outcome association estimate
#' @param se_exp SNP-exposure association standard error
#' @param se_out SNP-outome association standard error
#' @param threshold p-value threshold used to select SNPs
#' @export
#' @return List with the following elements:
#' \describe{
#' \item{b}{Causal effect estimate}
#' \item{se}{Causal effect standard error}
#' \item{pval}{Causal effect p-value}
#' \item{nsnps}{Number of instruments used}
#' }
#' @importFrom rlang .data
#' @importFrom magrittr "%>%"

raps_ZP <- function(b_exp, b_out, se_exp, se_out, threshold){
  # how many instruments do we need as a minimum?
  if(sum(!is.na(b_exp) & !is.na(b_out) &
         !is.na(se_exp) & !is.na(se_out)) < 5)
    stop("Not enough instruments (minimum is 5).")

  tmpfile = paste0("tmp_", runif(1, 1, 50000))
  input.file = paste0(tmpfile, ".tsv")
  output.file = paste0(tmpfile, "_out.tsv")
  dat = data.frame(beta_disc = b_exp,
                   se_disc = se_exp,
                   N_disc = NA_real_, # does not matter
                   af_disc = 0.5, # does not matter?
                   beta_rep = NA_real_,
                   se_rep = NA_real_,
                   N_rep = NA_real_,
                   af_rep = NA_real_)
  write.table(dat, input.file, sep="\t", quote=F, row.names=F)


  trait.mean <- 0 # does not matter?
  p.threshold <- threshold
  gwas.winners.curse::correct.winners.curse(input.file,
                                            output.file,
                                            trait.mean,
                                            p.threshold,
                                            header = TRUE,
                                            sep = "\t")
  corrected = data.table::fread(output.file)

  res_mr = mr.raps::mr.raps.overdispersed.robust(corrected$debiased.beta.mle, b_out, se_exp, se_out, loss.function = "huber")

  file.remove(input.file)
  file.remove(output.file)

  b = res_mr$beta.hat
  se = res_mr$beta.se
  pval <- 2 * stats::pnorm(abs(b/se), lower.tail=FALSE)

  return(list(b = b, se = se, pval = pval, nsnp=length(b_exp)))

}



#' MR-RAPS using UMVCUE correction
#'
#' @param b_disc SNP-exposure association estimate in discovery sample
#' @param b_rep SNP-exposure association estimate in replication sample
#' @param b_out SNP-outcome association estimate
#' @param se_disc SNP-exposure association standard error in discovery sample
#' @param se_rep SNP-exposure association standard error in replication sample
#' @param se_out SNP-outome association standard error
#' @param threshold p-value threshold used to select SNPs
#' @export
#' @return List with the following elements:
#' \describe{
#' \item{b}{Causal effect estimate}
#' \item{se}{Causal effect standard error}
#' \item{pval}{Causal effect p-value}
#' \item{nsnps}{Number of instruments used}
#' }
#' @importFrom rlang .data
#' @importFrom magrittr "%>%"

raps_UMVCUE <- function(b_disc, b_rep, b_out, se_disc, se_rep, se_out, threshold){
  # how many instruments do we need as a minimum?
  if(sum(!is.na(b_disc) & !is.na(b_rep) & !is.na(b_out) &
         !is.na(se_disc) & !is.na(se_rep) & !is.na(se_out)) < 5)
    stop("Not enough instruments (minimum is 5).")

  Tr = qnorm(threshold/2, lower.tail = F)
  m = length(b_disc)

  s1sq = se_disc^2
  s2sq = se_rep^2

  MLE = (s2sq*b_disc + s1sq*b_rep)/(s1sq+s2sq)

  BB1 = (sqrt(s1sq+s2sq)/s1sq)*(MLE-Tr*se_disc)
  BB2 = (sqrt(s1sq+s2sq)/s1sq)*(MLE+Tr*se_disc)
  betaX_UMVCUE = MLE - (s2sq/sqrt(s1sq+s2sq))*((dnorm(BB1)-dnorm(BB2))/(pnorm(BB1)+pnorm(-BB2)))

  EstBoot = matrix(nrow=200,ncol=m)
  for(k in 1:200){

    K = m

    a = Tr*sqrt(s1sq)
    b = rep(Inf, K)

    minus = stats::pnorm(-Tr-betaX_UMVCUE/sqrt(s1sq))/(stats::pnorm(-Tr+betaX_UMVCUE/sqrt(s1sq)) + stats::pnorm(-Tr-betaX_UMVCUE/sqrt(s1sq)))
    minus = (stats::runif(K) <= minus)

    a[minus] = -Inf
    b[minus] = -Tr*sqrt(s1sq[minus])

    BXGs = truncnorm::rtruncnorm(1, a, b, mean = betaX_UMVCUE, sd = sqrt(s1sq))
    BXG2s = stats::rnorm(K, mean = betaX_UMVCUE, sd = sqrt(s2sq))

    mle      = (s2sq*BXGs +  s1sq*BXG2s)/(s1sq+s2sq)
    BB1      = (sqrt(s1sq+s2sq)/s1sq)*(mle-Tr*sqrt(s1sq))
    BB2      = (sqrt(s1sq+s2sq)/s1sq)*(mle+Tr*sqrt(s1sq))
    EstBoot[k,] = mle - (s2sq/(sqrt(s1sq+s2sq)))*(stats::dnorm(BB1)-stats::dnorm(BB2))/(stats::pnorm(BB1)+stats::pnorm(-BB2))

  }

  betaX_UMVCUE_se = apply(EstBoot,2,sd)
  res_mr = mr.raps::mr.raps.overdispersed.robust(betaX_UMVCUE, b_out, betaX_UMVCUE_se, se_out, loss.function = "huber")

  b = res_mr$beta.hat
  se = res_mr$beta.se
  pval <- 2 * stats::pnorm(abs(b/se), lower.tail=FALSE)


  return(list(b = b, se = se, pval = pval, nsnp=length(b_disc)))

}

