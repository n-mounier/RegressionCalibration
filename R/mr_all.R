#' Perform Regression Calibration for all estimators (IVW, simple/weighted median, simple/weighted mode)
#'
#' @param dat Harmonised exposure (discovery and replication) and outcome data. Output from \code{TwoSampleMR::harmonise_data}.
#' Note than one extra column is needed (\code{sample}) and than only the following columns are used:
#' \code{SNP}, \code{beta.exposure}, \code{beta.outcome}, \code{se.exposure}, \code{se.outcome} and \code{mr_keep}
#' @export
#' @return data.frame with Regression Calibration results
#' For each method (IVW, simple/weighted median, simple/weighted mode)
#' \describe{
#' \item{nsnps}{Number of instruments used}
#' \item{b}{Causal effect estimate}
#' \item{se}{Causal effect standard error}
#' \item{pval}{Causal effect p-value}
#' }
#' @importFrom rlang .data



RC <- function(dat){
  # assumes than when harmonisation is done for multiple outcome, the SNP exposure coding is kept the same
  # i.e. alignment is done in a way that effect/other allele for the exposure do not change
  disc = dat$exposure[1]
  rep =  dat$outcome[dat$sample=="replication"][1]
  out = dat$outcome[dat$sample=="outcome"][1]
  cat("Exposure, discovery: ", disc, "\n")
  cat("Exposure, replication: ", rep, "\n")
  cat("Outcome: ", out, "\n")

  # disc - rep
  dat %>%
    dplyr::filter(.data$sample=="replication",
           .data$mr_keep) %>%
    dplyr::transmute(.data$SNP,
              b_disc = .data$beta.exposure,
              b_rep = .data$beta.outcome,
              se_disc = .data$se.exposure,
              se_rep = .data$se.outcome) -> dat_exp

  # disc - out
  dat %>%
    dplyr::filter(.data$sample=="outcome",
           .data$mr_keep)  %>%
    dplyr::transmute(.data$SNP,
              b_out = .data$beta.outcome,
              se_out = .data$se.outcome) -> dat_out

  dat_clean = dplyr::inner_join(dat_exp, dat_out, by="SNP") %>% dplyr::mutate(SNP=NULL)

  res = data.frame(exposure.discovery = disc,
                   exposure.replication = out,
                   outcome = out,
                   method = c("Inverse variance weighted", "Simple median", "Weighted median",
                              "Simple mode", "Weighted mode"),
                   b = NA_real_,
                   se = NA_real_,
                   pval = NA_real_,
                   nsnps = NA_real_)
  res[1, 5:8] = unlist(with(dat_clean, RegressionCalibration::RC_ivw(b_disc, b_rep, b_out, se_disc, se_rep, se_out)))
  res[2, 5:8] = unlist(with(dat_clean, RegressionCalibration::RC_simple_median(b_disc, b_rep, b_out, se_disc, se_rep, se_out)))
  res[3, 5:8] = unlist(with(dat_clean, RegressionCalibration::RC_weighted_median(b_disc, b_rep, b_out, se_disc, se_rep, se_out)))
  res[4, 5:8] = unlist(with(dat_clean, RegressionCalibration::RC_simple_mode(b_disc, b_rep, b_out, se_disc, se_rep, se_out)))
  res[5, 5:8] = unlist(with(dat_clean, RegressionCalibration::RC_weighted_mode(b_disc, b_rep, b_out, se_disc, se_rep, se_out)))

  res %>%
    dplyr::relocate(.data$nsnps, .after=.data$method) %>%
    dplyr::relocate(.data$pval, .after=.data$nsnps)-> res

  return(res)
}


#' IVW Regression Calibration
#'
#' @param b_disc SNP-exposure association estimate in discovery sample
#' @param b_rep SNP-exposure association estimate in replication sample
#' @param b_out SNP-outcome association estimate
#' @param se_disc SNP-exposure association standard error in discovery sample
#' @param se_rep SNP-exposure association standard error in replication sample
#' @param se_out SNP-outome association standard error
#' @export
#' @return List with the following elements:
#' \describe{
#' \item{b}{Causal effect estimate}
#' \item{se}{Causal effect standard error}
#' \item{pval}{Causal effect p-value}
#' \item{nsnps}{Number of instruments used}
#' }

RC_ivw <- function(b_disc, b_rep, b_out, se_disc, se_rep, se_out){
  # how many instruments do we need as a minimum?
  if(sum(!is.na(b_disc) & !is.na(b_rep) & !is.na(b_out) &
         !is.na(se_disc) & !is.na(se_rep) & !is.na(se_out)) < 5)
    stop("Not enough instruments (minimum is 5).")

  IVW_model = TwoSampleMR::mr_ivw(b_disc, b_out, se_disc, se_out)

  IVW_C = TwoSampleMR::mr_ivw(b_disc, b_rep, se_disc, se_rep)$b

  b <- IVW_model$b/IVW_C
  se <- sqrt(IVW_model$se**2 / IVW_C**2)
  pval <- 2 * stats::pnorm(abs(b/se), lower.tail=FALSE)

  return(list(b = b, se = se, pval = pval, nsnp=length(b_disc)))

}


#' Simple Median Regression Calibration
#'
#' @param b_disc SNP-exposure association estimate in discovery sample
#' @param b_rep SNP-exposure association estimate in replication sample
#' @param b_out SNP-outcome association estimate
#' @param se_disc SNP-exposure association standard error in discovery sample
#' @param se_rep SNP-exposure association standard error in replication sample
#' @param se_out SNP-outome association standard error
#' @export
#' @return List with the following elements:
#' \describe{
#' \item{b}{Causal effect estimate}
#' \item{se}{Causal effect standard error}
#' \item{pval}{Causal effect p-value}
#' \item{nsnps}{Number of instruments used}
#' }

RC_simple_median <- function(b_disc, b_rep, b_out, se_disc, se_rep, se_out){
  # how many instruments do we need as a minimum?
  if(sum(!is.na(b_disc) & !is.na(b_rep) & !is.na(b_out) &
         !is.na(se_disc) & !is.na(se_rep) & !is.na(se_out)) < 5)
    stop("Not enough instruments (minimum is 5).")

  MedianBased_model = TwoSampleMR::mr_simple_median(b_disc, b_out, se_disc, se_out)

  MedianBased_C <- TwoSampleMR::mr_simple_median(b_disc, b_rep, se_disc, se_rep, parameters = list(nboot=2))$b

  b <- MedianBased_model$b/MedianBased_C
  se <- sqrt(MedianBased_model$se**2 / MedianBased_C**2)
  pval <- 2 * stats::pnorm(abs(b/se), lower.tail=FALSE)

  return(list(b = b, se = se, pval = pval, nsnp=length(b_disc)))
}


#' Weighted Median Regression Calibration
#'
#' @param b_disc SNP-exposure association estimate in discovery sample
#' @param b_rep SNP-exposure association estimate in replication sample
#' @param b_out SNP-outcome association estimate
#' @param se_disc SNP-exposure association standard error in discovery sample
#' @param se_rep SNP-exposure association standard error in replication sample
#' @param se_out SNP-outome association standard error
#' @export
#' @return List with the following elements:
#' \describe{
#' \item{b}{Causal effect estimate}
#' \item{se}{Causal effect standard error}
#' \item{pval}{Causal effect p-value}
#' \item{nsnps}{Number of instruments used}
#' }

RC_weighted_median <- function(b_disc, b_rep, b_out, se_disc, se_rep, se_out){
  # how many instruments do we need as a minimum?
  if(sum(!is.na(b_disc) & !is.na(b_rep) & !is.na(b_out) &
         !is.na(se_disc) & !is.na(se_rep) & !is.na(se_out)) < 5)
    stop("Not enough instruments (minimum is 5).")

  MedianBased_model = TwoSampleMR::mr_weighted_median(b_disc, b_out, se_disc, se_out)

  MedianBased_C <- TwoSampleMR::mr_weighted_median(b_disc, b_rep, se_disc, se_rep, parameters = list(nboot=2))$b

  b <- MedianBased_model$b/MedianBased_C
  se <- sqrt(MedianBased_model$se**2 / MedianBased_C**2)
  pval <- 2 * stats::pnorm(abs(b/se), lower.tail=FALSE)

  return(list(b = b, se = se, pval = pval, nsnp=length(b_disc)))
}

#' Simple Mode Regression Calibration
#'
#' @param b_disc SNP-exposure association estimate in discovery sample
#' @param b_rep SNP-exposure association estimate in replication sample
#' @param b_out SNP-outcome association estimate
#' @param se_disc SNP-exposure association standard error in discovery sample
#' @param se_rep SNP-exposure association standard error in replication sample
#' @param se_out SNP-outome association standard error
#' @export
#' @return List with the following elements:
#' \describe{
#' \item{b}{Causal effect estimate}
#' \item{se}{Causal effect standard error}
#' \item{pval}{Causal effect p-value}
#' \item{nsnps}{Number of instruments used}
#' }

RC_simple_mode <- function(b_disc, b_rep, b_out, se_disc, se_rep, se_out){
  # how many instruments do we need as a minimum?
  if(sum(!is.na(b_disc) & !is.na(b_rep) & !is.na(b_out) &
         !is.na(se_disc) & !is.na(se_rep) & !is.na(se_out)) < 5)
    stop("Not enough instruments (minimum is 5).")

  param = TwoSampleMR::default_parameters()
  param$nboot = 2

  ModeBased_model = TwoSampleMR::mr_simple_mode(b_disc, b_out, se_disc, se_out)

  ModeBased_C <- TwoSampleMR::mr_simple_mode(b_disc, b_rep, se_disc, se_rep, parameters = param)$b

  b <- ModeBased_model$b/ModeBased_C
  se <- sqrt(ModeBased_model$se**2 / ModeBased_C**2)
  pval <- 2 * stats::pnorm(abs(b/se), lower.tail=FALSE)

  return(list(b = b, se = se, pval = pval, nsnp=length(b_disc)))
}


#' Weihgted Mode Regression Calibration
#'
#' @param b_disc SNP-exposure association estimate in discovery sample
#' @param b_rep SNP-exposure association estimate in replication sample
#' @param b_out SNP-outcome association estimate
#' @param se_disc SNP-exposure association standard error in discovery sample
#' @param se_rep SNP-exposure association standard error in replication sample
#' @param se_out SNP-outome association standard error
#' @export
#' @return List with the following elements:
#' \describe{
#' \item{b}{Causal effect estimate}
#' \item{se}{Causal effect standard error}
#' \item{pval}{Causal effect p-value}
#' \item{nsnps}{Number of instruments used}
#' }

RC_weighted_mode <- function(b_disc, b_rep, b_out, se_disc, se_rep, se_out){
  # how many instruments do we need as a minimum?
  if(sum(!is.na(b_disc) & !is.na(b_rep) & !is.na(b_out) &
         !is.na(se_disc) & !is.na(se_rep) & !is.na(se_out)) < 5)
    stop("Not enough instruments (minimum is 5).")

  param = TwoSampleMR::default_parameters()
  param$nboot = 2

  ModeBased_model = TwoSampleMR::mr_weighted_mode(b_disc, b_out, se_disc, se_out)

  ModeBased_C <- TwoSampleMR::mr_weighted_mode(b_disc, b_rep, se_disc, se_rep, parameters = param)$b

  b <- ModeBased_model$b/ModeBased_C
  se <- sqrt(ModeBased_model$se**2 / ModeBased_C**2)
  pval <- 2 * stats::pnorm(abs(b/se), lower.tail=FALSE)

  return(list(b = b, se = se, pval = pval, nsnp=length(b_disc)))
}

