#' Summary Statistics Column Headers
#'
#' @description List of uncorrected column headers often found in GWAS Summary
#' Statistics column headers. Note the effect allele will always be the A2
#' allele, this is the approach done for
#' VCF(https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7805039/).
#' (https://github.com/neurogenomics/MungeSumstats/blob/HEAD/R/data.R)
#' This is enforced with the column header corrections here and also the check allele flipping
#' test.
#'
#' @format dataframe with 2 columns
#' @usage data("sumstatsColHeaders")
"sumstatsColHeaders"

#' UCSC Chain file hg38 to hg19
#'
#' @description UCSC Chain file hg38 to hg19, .chain.gz file, downloaded from
#' https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/ on 09/10/21
#'
#' @name hg38ToHg19
#' @section hg38ToHg19.over.chain.gz
#' @details UCSC Chain file hg38 to hg19, .chain.gz file, downloaded on 09/10/21
#' To be used as a back up if the download from UCSC fails.
#' @source The chain file was downloaded from
#' https://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/
#' \code{
#' utils::download.file('ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz',tempdir())
#' }
#' @format gunzipped chain file
#' @details NULL
NULL

#' UCSC Chain file hg19 to hg38
#'
#' @description UCSC Chain file hg19 to hg38, .chain.gz file, downloaded from
#' https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/ on 09/10/21
#'
#' @name hg19ToHg38
#' @section hg19ToHg38.over.chain.gz
#' @details UCSC Chain file hg19 to hg38, .chain.gz file, downloaded on 09/10/21
#' To be used as a back up if the download from UCSC fails.
#' @source The chain file was downloaded from
#' https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/
#' \code{
#' utils::download.file('ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz',tempdir())
#' }
#' @format gunzipped chain file
#' @details NULL
NULL
