#' Function for obtaining top hits from a GWAS vcf file with and without clumping
#'
#' Obtain top hits from a GWAS vcf file with/without clumping
#' Note: gwasglue::clump_gwasvcf function clumping using 1kg reference genome or other reference genome UKbiobank,
#' the issue is that if a top variant is missing from the reference genome, it will it will "Removing top variants absence from LD reference panel",
#'  which is not ideal. As GWAS summary statistic are from imputed genome, for sure it will cover more SNPs than the unimputed reference genome.
#' The current solution is that get the top SNP in each chromosome (unclumped_top), then do the clumping using reference genome (clumped_df),
#' if (clumped_df) miss a chromosome in (unclumped_top), combine it into (clumped_df)
#'
#' @param path Specify the full path starting with /cloud-home/ to the vcf file, including the file name.
#' @param clump_p p value threshold to filter. Default = 5e-8
#' @param clump_r2 Clumping r2 threshold. Default: 0.001
#' @param clump_kb Clumping kb window. Default: 10000. (very strict)
#' @param ldref LD reference data. Normally Europeans from the 1000 genomes reference panel. Default:"/cloud-home/U1006974/NeuroID/data/1kg/EUR"
#' @param plink_bin  specify path to plink binary. Default = "genetics.binaRies::get_plink_binary()"
#' @param nThread  Number of threads to use for parallel processes
#'
#' @return a data.frame for top hits from a GWAS vcf file with/without clumping
#'
#' @references \url{https://mrcieu.github.io/gwasglue/reference/clump_gwasvcf.html}
#'
#' @export
#' @importFrom gwasvcf check_bcftools query_gwas vcf_to_tibble
#' @importFrom gwasglue clump_gwasvcf
#' @importFrom genetics.binaRies get_plink_binary
#'
tophits_vcf <-
  function(path,
           clump_p = 5e-08,
           clump_r2 = 0.001,
           clump_kb = 10000,
           ldref,
           plink_bin = genetics.binaRies::get_plink_binary(),
           nThread
           ) {
  unclumped <- gwasvcf::query_gwas(path, pval = clump_p, threads = nThread)
  message(paste0("Filtering ", path,  " by the threshold at ", clump_p, ", resulting in ", nrow(unclumped), " SNPs."))
  if(nrow(unclumped) != 0){
    unclumped_df <- gwasvcf::vcf_to_tibble(unclumped)
    unclumped_top <- unclumped_df %>%
      dplyr::group_by(seqnames) %>%
      dplyr::slice_max(LP) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(chrpos = paste0(seqnames, ":", start)) %>%
      dplyr::select(rsid, chrpos, seqnames)
    message(paste0("There is a top hit passing the threshold at ", clump_p, " in each Chromosome: ",
                   paste0(as.character(unique(unclumped_top$seqnames)), collapse = ", ")))
    #Perform clumping
    message("Perform clumping")
    clumped <- gwasglue::clump_gwasvcf(path,
                                       clump_p = clump_p,
                                       clump_r2 = clump_r2,
                                       clump_kb = clump_kb,
                                       bfile = ldref,
                                       plink_bin = plink_bin
                                       )
    clumped_df <- clumped %>%
      dplyr::mutate(seqnames = gsub("\\:.*", "", chrpos))
    message(paste0("Clumped SNPs based on reference genome are in chromosome: ",
                   paste0(as.character(unique(clumped_df$seqnames)), collapse = ", ")))
    missing_top = unclumped_top %>%
      dplyr::filter(seqnames %in% setdiff(unique(unclumped_top$seqnames),
                                          unique(clumped_df$seqnames)))
    message(paste0("Clumped SNPs based on reference genome are missing tophits from chromosome: ",
                   paste0(as.character(missing_top$seqnames), collapse = ", ")))
    clumped_df <- clumped_df %>%
      dplyr::bind_rows(missing_top) %>%
      dplyr::distinct() %>%
      dplyr::left_join(unclumped_df)
    message(paste0("Clumping SNPs based on reference genome and each chromosome tophits are from chromosome: ",
                   paste0(as.character(unique(clumped_df$seqnames)), collapse = ", ")))
    return(list(clumped = clumped_df,
           unclumped = unclumped_df))
  }else{
    return(list(clumped = NA,
           unclumped = NA))
  }
}
