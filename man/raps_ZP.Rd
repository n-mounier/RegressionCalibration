% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mr_others.R
\name{raps_ZP}
\alias{raps_ZP}
\title{MR-RAPS using ZP correction}
\usage{
raps_ZP(b_exp, b_out, se_exp, se_out, threshold)
}
\arguments{
\item{b_exp}{SNP-exposure association estimate}

\item{b_out}{SNP-outcome association estimate}

\item{se_exp}{SNP-exposure association standard error}

\item{se_out}{SNP-outome association standard error}

\item{threshold}{p-value threshold used to select SNPs}
}
\value{
List with the following elements:
\describe{
\item{b}{Causal effect estimate}
\item{se}{Causal effect standard error}
\item{pval}{Causal effect p-value}
\item{nsnps}{Number of instruments used}
}
}
\description{
MR-RAPS using ZP correction
}
