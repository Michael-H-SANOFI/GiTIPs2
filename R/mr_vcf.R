#' Function for Mendelian randomization using vcf file
#'
#' Extract instruments from exposure vcf file by a pvalue threshold and do clumping
#' and then extract them from outcome vcf file, format both datasets into TwosampleMR
#' run Mendelian randomization using default methods, "mr_wald_ratio", "mr_egger_regression" and "mr_ivw"
#'
#' @param exposure_vcf Specify the full path starting with /cloud-home/ to the vcf file, including the file name, for exposure. Normally a molecular qtl file
#' @param outcome_vcf  Specify the full path starting with /cloud-home/ to the vcf file, including the file name, for outcome. Normally a gwas file
#' @param pval p value threshold to extract the instruments from exposure vcf file Default = 5e-8
#' @param clump_r2 Clumping r2 threshold. Default: 0.001
#' @param clump_kb Clumping kb window. Default: 10000. (very strict)
#' @param ldref LD reference data. Normally Europeans from the 1000 genomes reference panel. Default:"/cloud-home/U1006974/NeuroID/data/1kg/EUR"
#' @param plink_bin Specify path to plink binary. Default = "genetics.binaRies::get_plink_binary()"
#' @param exclude_HLA Exclude SNPs in HLA region or not. Default: TRUE. GRCh37; Location: chr6:28,477,797-33,448,354;
#' @param mr_method methods to use in mr analysis. See [mr_method_list()] for details. Default: subset(mr_method_list(), use_by_default)$obj[c(1,2,4)].
#' it will contain "mr_wald_ratio", "mr_egger_regression" and "mr_ivw"
#' @param nThread  Number of threads to use for parallel processes
#'
#' @return a data.frame for mr results
#'
#' @references \url{https://mrcieu.github.io/gwasglue/articles/mr.html}
#'
#' @export
#' @importFrom gwasvcf query_gwas
#' @importFrom gwasglue clump_gwasvcf gwasvcf_to_TwoSampleMR
#' @importFrom TwoSampleMR mr_method_list harmonise_data mr
#' @importFrom ieugwasr ld_clump
#' @importFrom genetics.binaRies get_plink_binary
#'
mr_vcf <- function(exposure_vcf,
                   outcome_vcf,
                   pval = 5e-08,
                   clump_r2 = 0.001,
                   clump_kb = 10000,
                   ldref,
                   plink_bin = genetics.binaRies::get_plink_binary(),
                   exclude_HLA = TRUE,
                   mr_method = subset(TwoSampleMR::mr_method_list(), use_by_default)$obj[c(1,2,4)],
                   nThread
                   ) {
  #Extract clumped instruments from exposure vcf file
  message(paste0("Extract clumped instruments from exposure vcf file ", exposure_vcf))
  exposure_iv = tophits_vcf(path = exposure_vcf,
                            clump_p = pval,
                            clump_r2 = clump_r2,
                            clump_kb = clump_kb,
                            ldref = ldref,
                            plink_bin = plink_bin,
                            nThread = nThread)
  exposure_iv = exposure_iv$clumped
  if(nrow(exposure_iv) != 0){
    if(exclude_HLA){
      #remove SNPs in HLA region
      snp_hla <- exposure_iv %>%
        dplyr::filter((seqnames == 6 & between(start, 28477797, 33448354)))
      message(paste0("Remove ", nrow(snp_hla), " SNPs in HLA region."))
      exposure_iv <- exposure_iv %>%
        dplyr::filter(!(seqnames == 6 & between(start, 28477797, 33448354)))
    }
    if(nrow(exposure_iv) != 0){
      #Extracting outcome data
      message("Extract SNPs from outcome vcf file and convert to TwoSampleMR format")
      outcome_iv_vcf <- gwasvcf::query_gwas(vcf     = outcome_vcf,
                                            rsid    = exposure_iv$rsid,
                                            proxies = "no",
                                            bfile   = ldref,
                                            threads = nThread)
      if(nrow(outcome_iv_vcf) != 0){
        exposure_dat = vcftibble_to_TwoSampleMR(exposure_iv, "exposure")
        outcome_dat <- gwasglue::gwasvcf_to_TwoSampleMR(outcome_iv_vcf, "outcome")
        #Proceed with harmonising
        message("Proceed with harmonising!")
        harmonised_dat <- TwoSampleMR::harmonise_data(exposure_dat = exposure_dat,
                                                      outcome_dat  = outcome_dat)
        if(nrow(harmonised_dat) != 0){
          message(paste0("Perform MR using methods: ",
                         paste0(mr_method, collapse = ", ")))
          mr_res <- TwoSampleMR::mr(harmonised_dat,
                                    method_list = mr_method)
          return(mr_res)
        }else{
          message("No SNPs left after harmonising!")
          return(NA)
        }
      }else{
        message("No SNPs found in the outcome vcf!")
        return(NA)
      }
    }else{
      message("No SNPs left after exclue HLA!")
      return(NA)
    }
  }else{
    message("No SNPs found after clumping in the exposure vcf file!")
    return(NA)
  }
}
