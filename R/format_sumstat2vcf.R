#' General function for ingesting GWAS summary statistics into GWAS-VCF format
#'
#' Create a standardized VCF file
#'
#' @param path the file name with a path which the summary statistics is to be read from.This can be a compressed file as *.tsv.gz
#' @param SNP column name for SNP (rs id). If NULL, rs id inferred from reference genome by MungeSumstat
#' @param CHR column name for chromosome. Required.
#' @param BP  column name for base pair positions. Required.
#' @param A1  column name for reference allele. If NULL, inferred by MungeSumstat.
#' @param A2  column name for alternative allele. Required.  A note on the allele flipping check: MungeSumstats infers the effect allele will always be the A2 allele,
#' this is the approach done for IEU GWAS VCF and has such also been adopted here.
#' @param FRQ column name for alternate allele frequency of the SNP.If there is no column for FRQ, set as NULL, inferred from 1kg reference
#' @param BETA column name for effect size estimate relative to the alternative allele. If NULL, MungeSumstats sets impute_beta = TRUE
#' @param SE column name for standard error of effect size estimate. If NULL, then MungeSumstats sets impute_se = TRUE
#' @param P column name for unadjusted p-value of effect estimate. Required.
#' @param NCAS column name for number of cases, if studytype is "Binary",Required.If there is no column for NCAS, input an integer number. if studytype is "Continuous", NA.
#' @param N column name for total sample size. Required.
#' If there is no column for sample size, then input an integer, MungeSumstats sets compute_n to this integer value the N (sample size) for every SNP in the dataset
#' @param ref_genome 	name of the reference genome used for the GWAS ("GRCh37" or "GRCh38") Argument is case-insensitive. Required.
#' @param id name of the vcf output file and used in the gwas catalogue as ID identifier
#' @param studytype "Continuous" for quantitative trait and "Binary" for binary trait.
#' @param save_path  the path which the vcf will be written out. End with "/"
#' @param hg19_1kg_phase3 if FRQ is NULL, annotate the FRQ with the 1000 genomes phase 3 (hg19/GRCh37) allele frequencies and index
#' by "/cloud-home/U1006974/GiTIPs/data/reference/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.AF.tsv"
#' @param trait        trait name. Required.
#' @param group_name   group name. Optional. "public" or "private" or ...
#' @param year         Optional.
#' @param author       Optional.
#' @param consortium   Optional.
#' @param sex          Optional.
#' @param pmid         Optional.
#' @param population   Optional.
#' @param unit         Optional.
#' @param subcategory  Optional.
#' @param ontology     Optional.
#' @param note         Optional.
#' @param doi          Optional.
#' @param sd           Optional.
#' @param nThread  Number of threads to use for parallel processes
#' @param ... parameters from MungeSumstats. See the reference below
#'
#' @return one .vcf.gz file, one .vcf.gz.tbi file, one a log file .log, and one catalog file, .clg.
#' Note: The majority of VCF handling bioinformatics libraries use a .vcf.gz suffix, even for block gziped output.
#' However, VariantAnnotation::writeVcf() even with index = TRUE does not support this and forceably sets the suffix to .bgz.
#' See https://github.com/Bioconductor/VariantAnnotation/issues/35.
#' So the current solution is to rename from vcf.bgz to vcf.gz
#'
#' @references \url{https://github.com/MRCIEU/gwas-vcf-specification}
#' @references \url{https://neurogenomics.github.io/MungeSumstats/reference/index.html}
#'
#' @export
#' @importFrom MungeSumstats format_sumstats
#' @importFrom MungeSumstats read_header
#' @importFrom dplyr add_row filter rename mutate as_tibble left_join all_of if_else
#' @importFrom gwasvcf create_vcf
#' @importFrom data.table fread as.data.table
#' @importFrom data.table ":="
#' @importFrom vroom vroom
#' @importFrom magrittr "%>%"
#' @importFrom magrittr "%$%"
#' @importFrom VariantAnnotation meta header writeVcf
#' @importFrom utils write.table
#'
format_sumstat2vcf <-
  function(path,
           SNP,
           CHR,
           BP,
           A1,
           A2,
           FRQ,
           BETA,
           SE,
           P,
           NCAS,
           N,
           ref_genome,
           id,
           studytype,
           save_path,
           hg19_1kg_phase3,
           trait,
           group_name  = NA,
           year        = NA,
           author      = NA,
           consortium  = NA,
           sex         = NA,
           pmid        = NA,
           population  = NA,
           unit        = NA,
           subcategory = NA,
           ontology    = NA,
           note        = NA,
           doi         = NA,
           sd          = NA,
           nThread,
           ...) {
  #save message to log and show in the screen
  log.file = paste0(save_path, id, ".log")
  msg <- paste0("Saving output messages to:\n", log.file)
  message(msg)
  msgcon <- file(log.file, open = "a")
  sink(file = log.file, split = TRUE, append = FALSE)
  sink(msgcon, type = "message",  append = FALSE)
  #check required parameters
  if(is.null(CHR) ||
     is.null(BP)  ||
     is.null(A2)  ||
     is.null(P) ||
     is.null(N) ||
     is.null(ref_genome) ||
     is.null(id) ||
     is.null(studytype) ||
     is.null(trait)
     ){
    stp <- "Please check required parameters: CHR, BP, A2, P, N, ref_genome, id, studytype, trait!"
    stop(stp)
  }
  #check for NCAS
  if(studytype == "Binary" & is.na(NCAS)){
    stp <- "Studytype is Binary, NCAS is required!"
    stop(stp)
  }
  #for the required columns update the column mapping file
  #data(sumstatsColHeaders)
  sumstatsColHeaders_new <- sumstatsColHeaders %>%
    dplyr::add_row(Uncorrected  = CHR, Corrected = "CHR") %>%
    dplyr::add_row(Uncorrected  = BP,  Corrected = "BP") %>%
    dplyr::add_row(Uncorrected  = P,   Corrected = "P")
  #SNP
  if(!is.null(SNP)){
    sumstatsColHeaders_new <- sumstatsColHeaders_new %>%
      dplyr::add_row(Uncorrected  = SNP, Corrected = "SNP")
  }
  #when both A1 and A2 not NULL
  if(!is.null(A1) & !is.null(A2)){
    sumstatsColHeaders_new <- sumstatsColHeaders_new %>%
      dplyr::add_row(Uncorrected  = A1, Corrected = "A1") %>%
      dplyr::add_row(Uncorrected  = A2, Corrected = "A2")
  }
  #FRQ
  if(!is.null(FRQ)){
    sumstatsColHeaders_new <- sumstatsColHeaders_new %>%
      dplyr::add_row(Uncorrected  = FRQ, Corrected = "FRQ")
  }
  #BETA
  if(!is.null(BETA)){
    sumstatsColHeaders_new <- sumstatsColHeaders_new %>%
      dplyr::add_row(Uncorrected  = BETA, Corrected = "BETA")
    impute_beta = FALSE
  }else{
    impute_beta = TRUE
  }
  #SE
  if(!is.null(SE)){
    sumstatsColHeaders_new <- sumstatsColHeaders_new %>%
      dplyr::add_row(Uncorrected  = SE, Corrected = "SE")
    impute_se = FALSE
  }else{
    impute_se = TRUE
  }
  #N
  if(is.character(N)){
    sumstatsColHeaders_new <- sumstatsColHeaders_new %>%
      dplyr::add_row(Uncorrected  = N, Corrected = "N")
    compute_n = 0
  }else{
    compute_n = N
  }
  #### Check if tabular 1: infer from file name ####
  tab_suffixes <- supported_suffixes(vcf = FALSE,
                                     vcf_compressed = FALSE)
  tab_suffix_regexes <- gsub("\\.", "\\.",  paste0(tab_suffixes, "$"))
  is_tabular <- grepl(paste(tab_suffix_regexes, collapse = "|"), path)
  #### Check if tabular 2: infer from data ####
  if(isFALSE(is_tabular)){
    header <- MungeSumstats::read_header(path = path)
    is_tabular <- check_tabular(header = header)
  }
  #### Process tabular ####
  if (isTRUE(is_tabular)) {
    if(endsWith(path,".bgz")){
      message("Importing tabular bgz file: ", path)
      sumstats_file <- data.table::fread(
        text = readLines(con = path),
        nThread = nThread,
        nrows = Inf)
    }else {
      message("Importing tabular file: ", path, " using ", nThread, " threads")
      sumstats_file <- vroom::vroom(
        path,
        num_threads = nThread,
        progress = FALSE)
    }

  } else {
    suffixes <- supported_suffixes()
    stop(
      "Unrecognized file format.\n",
      "Must be one of: \n   ",
      paste(suffixes, collapse = "\n   ")
    )
  }
  #rename CHR
  sumstats_file <- sumstats_file %>%
    data.table::as.data.table(.) %>%
    dplyr::rename("CHR" = all_of(CHR)) %>%
    dplyr::mutate(CHR = as.character(CHR))
  message(paste0("chromosome style is in ", paste0(unique(sumstats_file$CHR),collapse = ",")))
  #rename 23 to X
  if("23" %in% unique(sumstats_file$CHR)){
    message("Rename 23 to X")
    sumstats_file[, CHR := gsub("23", "X", CHR)]
  }
  #rename M to MT
  if("M" %in% unique(sumstats_file$CHR)){
    message("Rename M to MT")
    sumstats_file[, CHR := gsub("M", "MT", CHR, ignore.case = TRUE)]
  }
  #when A1 is NULL and A2 is not NULL
  if(is.null(A1) & !is.null(A2)){ # remove mapping of A1/A2 to avoid 'ref' match 'alt'
    message("Remove mapping of A2 to avoid REF match ALT")
    #rename column A2 if column A2 is not NULL
    sumstats_file <- sumstats_file %>%
      dplyr::rename(UNCERTAIN_A2 = all_of(A2)) %>%
      dplyr::mutate(UNCERTAIN_A2 = toupper(UNCERTAIN_A2))
  }
  #NCAS
  if(studytype == "Binary"){
    if(is.numeric(NCAS)){
      sumstats_file <- sumstats_file %>%
        dplyr::mutate(NCAS = NCAS)
    }else{
      sumstats_file <- sumstats_file %>%
        dplyr::rename("NCAS" = all_of(NCAS))
    }
  }
  #format using MungeStat
  if(ref_genome == "GRCh37"){
    requireNamespace(SNPlocs.Hsapiens.dbSNP155.GRCh37)
    sumstats_dt <- MungeSumstats::format_sumstats(path = sumstats_file,
                                                  impute_beta = impute_beta,
                                                  impute_se = impute_se,
                                                  compute_n = compute_n,
                                                  rmv_chr = NULL,
                                                  dbSNP = "155",
                                                  ref_genome = "GRCh37",
                                                  drop_indels = TRUE,
                                                  return_data = TRUE,
                                                  return_format = "data.table",
                                                  nThread = nThread,
                                                  mapping_file = sumstatsColHeaders_new,
                                                  ...)
  }
  #format using MungeStat
  if(ref_genome == "GRCh38"){#liftover from hg38 to hg19
    requireNamespace(SNPlocs.Hsapiens.dbSNP155.GRCh38)
    sumstats_dt <- MungeSumstats::format_sumstats(path = sumstats_file,
                                                  impute_beta = impute_beta,
                                                  impute_se = impute_se,
                                                  compute_n = compute_n,
                                                  rmv_chr = NULL,
                                                  dbSNP = "155",
                                                  ref_genome = "GRCh38",
                                                  convert_ref_genome = "GRCh37",
                                                  drop_indels = TRUE,
                                                  return_data = TRUE,
                                                  return_format = "data.table",
                                                  nThread = nThread,
                                                  mapping_file = sumstatsColHeaders_new,
                                                  ...)
  }
  message(paste0("chromosome style is in ", paste0(unique(sumstats_dt$CHR),collapse = ",")))
  #convert to vcf
  if(is.null(FRQ)){
    message("Annotate vcf file with 1kg allele frequencies")
    AF <- vroom::vroom(hg19_1kg_phase3, col_names = c("CHR", "BP", "A1", "A2", "FRQ"), num_threads = nThread, progress = FALSE)
    out <- sumstats_dt %>%
      dplyr::as_tibble() %>%
      dplyr::left_join(AF)
    message(paste0(sum(is.na(out$FRQ)), " SNP's FRQ is missing"))
    #when A1 is NULL and A2 is not NULL
    if(is.null(A1) & !is.null(A2)){
      message(paste0("A1 and A2 inferred from reference genome, do another round of allele_flip_check to match A2 with ", A2))
      #A1 and A2 both derived from reference genome in the sumstats_dt
      #drop rows UNCERTAIN_A2 not match either A1 or A2
      if(nrow(sumstats_dt[A1 != UNCERTAIN_A2 &
                          A2 != UNCERTAIN_A2, ]) > 0) {
        print_msg0 <-
          paste0(
            "There are ",
            formatC(nrow(sumstats_dt[A1 != UNCERTAIN_A2 &
                                       A2 != UNCERTAIN_A2, ]), big.mark = ","),
            " SNPs where ", A2, " doesn't match neither A1 nor A2 in the reference genome.",
            " These will be removed."
          )
        message(print_msg0)
      }
      if(nrow(sumstats_dt[A2 != UNCERTAIN_A2, ]) > 0) {
        print_msg0 <-
          paste0(
            "There are ",
            formatC(nrow(sumstats_dt[A2 != UNCERTAIN_A2, ]), big.mark = ","),
            paste0(" SNPs where ", A2, " doesn't match A2 in the reference genome.",
                   " These will be flipped.")
          )
        message(print_msg0)
        #if UNCERTAIN_A2 matches A2, leave it
        #if UNCERTAIN_A2 matches A1(ref allele), flip it
        sumstats_dt <- sumstats_dt %>%
          dplyr::mutate(BETA = if_else(UNCERTAIN_A2 == A1, -BETA, BETA))
      }
    }
    #NCAS
    if(!is.null(NCAS)){
      out <- out %$%
        gwasvcf::create_vcf(chrom = CHR, pos = BP, nea = A1, ea = A2, snp = SNP, effect = BETA, ea_af = FRQ, se = SE, pval = P, ncase = NCAS, n = N, name = id)
    }else{
      out <- out %$%
        gwasvcf::create_vcf(chrom = CHR, pos = BP, nea = A1, ea = A2, snp = SNP, effect = BETA, ea_af = FRQ, se = SE, pval = P, n = N, name = id)
    }
  }else{
    #when A1 is NULL and A2 is not NULL
    if(is.null(A1) & !is.null(A2)){
      message(paste0("A1 and A2 inferred from reference genome, do another round of allele_flip_check to match A2 with ", A2))
      #A1 and A2 both derived from reference genome in the sumstats_dt
      #drop rows UNCERTAIN_A2 not match either A1 or A2
      if(nrow(sumstats_dt[A1 != UNCERTAIN_A2 &
                          A2 != UNCERTAIN_A2, ]) > 0) {
        print_msg0 <-
          paste0(
            "There are ",
            formatC(nrow(sumstats_dt[A1 != UNCERTAIN_A2 &
                                       A2 != UNCERTAIN_A2, ]), big.mark = ","),
            " SNPs where ", A2, " doesn't match neither A1 nor A2 in the reference genome.",
            " These will be removed."
          )
        message(print_msg0)
      }
      if(nrow(sumstats_dt[A2 != UNCERTAIN_A2, ]) > 0) {
        print_msg0 <-
          paste0(
            "There are ",
            formatC(nrow(sumstats_dt[A2 != UNCERTAIN_A2, ]), big.mark = ","),
            paste0(" SNPs where ", A2, " doesn't match A2 in the reference genome.",
                   "\nThese will be flipped with their effect column and FRQ column ")
          )
        message(print_msg0)
        #if UNCERTAIN_A2 matches A2, leave it
        #if UNCERTAIN_A2 matches A1(ref allele), flip it
        sumstats_dt <- sumstats_dt %>%
          dplyr::mutate(BETA = if_else(UNCERTAIN_A2 == A1, -BETA, BETA)) %>%
          dplyr::mutate(FRQ =  if_else(UNCERTAIN_A2 == A1, 1-FRQ, FRQ))
      }
    }
    #NCAS
    if(!is.null(NCAS)){
      out <- sumstats_dt %>%
        dplyr::as_tibble() %$%
        gwasvcf::create_vcf(chrom = CHR, pos = BP, nea = A1, ea = A2, snp = SNP, effect = BETA, ea_af = FRQ, se = SE, pval = P, ncase = NCAS, n = N, name = id)
    }else{
      out <- sumstats_dt %>%
        dplyr::as_tibble() %$%
        gwasvcf::create_vcf(chrom = CHR, pos = BP, nea = A1, ea = A2, snp = SNP, effect = BETA, ea_af = FRQ, se = SE, pval = P, n = N, name = id)
    }
  }
    message(paste0("studytype is ", studytype))
    VariantAnnotation::meta(VariantAnnotation::header(out))[["SAMPLE"]][["StudyType"]] = studytype
    rownames(VariantAnnotation::meta(header(out))[["SAMPLE"]]) = id
    message(paste0("Writing into a vcf file ", paste0(save_path, id, ".vcf.gz")))
    VariantAnnotation::writeVcf(out, file = paste0(save_path, id, ".vcf"), index = TRUE)
    #see note
    #rename from vcf.bgz to vcf.gz
    file.rename(from = paste0(save_path, id, ".vcf.bgz"),
                to   = paste0(save_path, id, ".vcf.gz"))
    #rename from vcf.bgz.tbi to vcf.gz.tbi
    file.rename(from = paste0(save_path, id, ".vcf.bgz.tbi"),
                to   = paste0(save_path, id, ".vcf.gz.tbi"))
    ##create a catalog file matching with opengwas
    #mr, coverage, study_design, priority in gwasinfo are removed
    clg  <-  data.frame(
      id,
      trait,
      ncase = NA,
      group_name,
      year,
      author,
      consortium,
      sex,
      pmid,
      population,
      unit,
      sample_size = max(sumstats_dt$N),
      build = "HG19/GRCh37",
      ncontrol = NA,
      category = studytype,
      subcategory,
      ontology,
      note,
      nsnp = nrow(sumstats_dt),
      doi,
      sd
    )
    #ncase
    if(("NCAS" %in% colnames(sumstats_dt))){
      clg$ncase = max(sumstats_dt$NCAS)
    }
    #ncontrol
    if(!is.na(clg$ncase)){
      clg$ncontrol = clg$sample_size - clg$ncase
    }
    message(paste0("Writing into a catalog file ", paste0(save_path, id, ".clg")))
    write.table(clg, file = paste0(save_path, id, ".clg"), sep = "\t",
                col.names = TRUE, row.names = FALSE, quote = FALSE)
    #close the external connection
    sink(NULL, type = "message")
    sink(NULL, type = "output")
}
