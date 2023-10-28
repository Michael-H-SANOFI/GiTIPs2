#' Function for variant annotation from vcf file
#'
#' A wrapper for variant annotation by VEP via bioconductor ensemblVEP
#'
#' @param path Specify the full path starting with /cloud-home/ to the vcf file, including the file name.
#' @param pval p value threshold to filter. Default = 5e-8
#' @param dir_cache Specify the cache directory to use. Required. Default: "/cloud-data/its-cmo-darwin-magellan-workspaces-folders/WS_genetic_ts/PMCB_Genetics/data/VEP/110"
#' @param cache_version a numeric. Use a different cache version than the assumed default (the VEP version). Required. Default: 110
#' @param nThread  Number of threads to use for parallel processes
#'
#' @return a data.frame for variant annotation. ensemblVEP will return 82 columns from VEP.
#'
#' @references \url{https://grch37.ensembl.org/info/docs/tools/vep/script/vep_options.html}
#' @references \url{https://bioconductor.org/packages/release/bioc/vignettes/ensemblVEP/inst/doc/ensemblVEP.pdf}
#'
#' @export
#' @importFrom dplyr as_tibble
#' @importFrom gwasvcf check_bcftools query_gwas
#' @importFrom ensemblVEP VEPFlags ensemblVEP
#'
annotate_vcf_vep <- function(path,
                             pval = 5e-08,
                             dir_cache = "/cloud-data/its-cmo-darwin-magellan-workspaces-folders/WS_genetic_ts/PMCB_Genetics/data/VEP/110",
                             cache_version = 110,
                             nThread) {
  #check and set up bcftools
  stopifnot(gwasvcf::check_bcftools())
  bcftools <- options()[["tools_bcftools"]]
  #query first
  input_file_filtered <- gwasvcf::query_gwas(path, pval = pval)
  if(nrow(input_file_filtered) != 0){
    #filter the input vcf file based on pval
    #create a temp output vcf file
    tempvcf <- tempfile(pattern = "file", tmpdir = tempdir(), fileext = ".vcf")
    #bcftools command for filtering on pval and output to a tempvcf file
    message(paste0("Filtering ", path,  " by the threshold at ", pval, " with ", nrow(input_file_filtered), " SNPs written into a temp vcf file: ", tempvcf))
    cmd <- sprintf("%s filter -i 'FORMAT/LP > %s' -o %s %s", bcftools, -log10(pval), tempvcf, path)
    system(cmd)
    message(paste0("Annotate the temp vcf file: ", tempvcf, " with everything by VEP cache_version: ", cache_version))
    #param for ensemblVEP
    param <- ensemblVEP::VEPFlags(flags = list(
      cache = TRUE,
      dir_cache = dir_cache,
      cache_version = cache_version,
      everything = TRUE,
      force_overwrite = TRUE,
      check_existing = TRUE,
      vcf = FALSE, #results returned as GRanges
      offline = TRUE,
      fork  = nThread
    ))
    gr <- ensemblVEP::ensemblVEP(file = tempvcf,
                     param = param)
    annotation = gr %>%
      dplyr::as_tibble()
    return(annotation = annotation)
  }else{
    message(paste0("No hits passing the p-value threshold at ", pval))
    return(annotation = NULL)
  }
}
