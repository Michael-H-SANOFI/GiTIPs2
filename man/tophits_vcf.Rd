% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tophits_vcf.R
\name{tophits_vcf}
\alias{tophits_vcf}
\title{Function for obtaining top hits from a GWAS vcf file with and without clumping}
\usage{
tophits_vcf(
  path,
  clump_p = 5e-08,
  clump_r2 = 0.001,
  clump_kb = 10000,
  ldref,
  plink_bin = genetics.binaRies::get_plink_binary(),
  nThread
)
}
\arguments{
\item{path}{Specify the full path starting with /cloud-home/ to the vcf file, including the file name.}

\item{clump_p}{p value threshold to filter. Default = 5e-8}

\item{clump_r2}{Clumping r2 threshold. Default: 0.001}

\item{clump_kb}{Clumping kb window. Default: 10000. (very strict)}

\item{ldref}{LD reference data. Normally Europeans from the 1000 genomes reference panel. Default:"/cloud-home/U1006974/NeuroID/data/1kg/EUR"}

\item{plink_bin}{specify path to plink binary. Default = "genetics.binaRies::get_plink_binary()"}

\item{nThread}{Number of threads to use for parallel processes}
}
\value{
a data.frame for top hits from a GWAS vcf file with/without clumping
}
\description{
Obtain top hits from a GWAS vcf file with/without clumping
Note: gwasglue::clump_gwasvcf function clumping using 1kg reference genome or other reference genome UKbiobank,
the issue is that if a top variant is missing from the reference genome, it will it will "Removing top variants absence from LD reference panel",
 which is not ideal. As GWAS summary statistic are from imputed genome, for sure it will cover more SNPs than the unimputed reference genome.
The current solution is that get the top SNP in each chromosome (unclumped_top), then do the clumping using reference genome (clumped_df),
if (clumped_df) miss a chromosome in (unclumped_top), combine it into (clumped_df)
}
\references{
\url{https://mrcieu.github.io/gwasglue/reference/clump_gwasvcf.html}
}
