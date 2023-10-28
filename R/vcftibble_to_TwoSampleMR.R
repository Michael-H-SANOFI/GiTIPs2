#' Function for creating exposure or outcome data format for TwoSampleMR from a tibble/data.frame
#'
#' creating exposure or outcome data format for TwoSampleMR from a tibble
#' which is created by vcf_to_granges(vcf) %>% dplyr::as_tibble()
#'
#' @param df a data.frame or a tibble
#' @param type "exposure" or "outcome"
#'
#' @return a data.frame for format of TwoSampleMR
#'
#' @export
#' @importFrom TwoSampleMR format_data
#'
vcftibble_to_TwoSampleMR <- function(df, type){
  if (!"ES" %in% names(df))
    df[["ES"]] <- NA_real_
  if (!"SE" %in% names(df))
    df[["SE"]] <- NA_real_
  if (!"LP" %in% names(df))
    df[["LP"]] <- NA_real_
  if (!"SS" %in% names(df))
    df[["SS"]] <- NA_real_
  if (!"NC" %in% names(df))
    df[["NC"]] <- NA_real_
  if (!"id" %in% names(df))
    df[["id"]] <- NA
  df[["LP"]] <- 10^-df[["LP"]]
  df[["NCONT"]] <- df[["SS"]] - df[["NC"]]
  TwoSampleMR::format_data(df, type = type, snp_col = "ID",
                           effect_allele_col = "ALT", other_allele_col = "REF",
                           eaf_col = "AF", chr_col = "seqnames", pos_col = "start",
                           beta_col = "ES", se_col = "SE", pval_col = "LP", samplesize_col = "SS",
                           ncase_col = "NC", ncontrol_col = "NCONT", phenotype_col = "id")
}
