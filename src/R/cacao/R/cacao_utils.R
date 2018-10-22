#' Function that initiates the cacao reporting object
#' @param mode mode for reporting (somatic, hereditary, any)
#' @param sample_name name of sample
#' @param version cacao version
#' @param genome_assembly human genome assembly (grch37/grch38)
#' @param mapq_threshold alignment quality (MAPQ) filter
#' @param callability_levels_germline array of three integers denoting thresholds for LOW_COVERAGE, CALLABLE, and HIGH_COVERAGE (germline setting)
#' @param callability_levels_somatic array of three integers denoting thresholds for LOW_COVERAGE, CALLABLE, and HIGH_COVERAGE (somatic setting)
#'
#' @return coverage_df
#'
init_report <- function(mode, sample_name, version, genome_assembly, mapq_threshold, callability_levels_germline, callability_levels_somatic){

  cacao_report <- list()
  cacao_report[['sample_name']] <- sample_name
  cacao_report[['version']] <- version
  cacao_report[['mode']] <- mode
  cacao_report[['genome_assembly']] <- genome_assembly
  cacao_report[['mapq']] <- mapq_threshold
  cacao_report[['loci']] <- list()
  cacao_report[['global_distribution']] <- list()
  cacao_report[['eval']] <- list()

  cacao_report[['callability']] <- list()
  cacao_report[['callability']][['levels']] <- list()
  cacao_report[['callability']][['verbose']] <- list()
  cacao_report[['callability']][['levels']][['germline']] <- callability_levels_germline
  cacao_report[['callability']][['levels']][['somatic']] <- callability_levels_somatic

  for(c in c('somatic','germline')){
    cacao_report[['callability']][['verbose']][[c]] <- list()
    for(level in c('no_coverage','low','callable','high')){
      if(level == 'no_coverage'){
        cacao_report[['callability']][['verbose']][[c]][['no_coverage']] <- 'No coverage (zero sequencing depth)'
      }
      if(level == 'low'){
        cacao_report[['callability']][['verbose']][[c]][['low']] <- paste0('Sequencing depth from 1 to ',cacao_report[['callability']][['levels']][[c]][2] - 1)
      }
      if(level == 'callable'){
        cacao_report[['callability']][['verbose']][[c]][['callable']] <- paste0('Sequencing depth from ',cacao_report[['callability']][['levels']][[c]][2],' to ',cacao_report[['callability']][['levels']][[c]][3] - 1)
      }
      if(level == 'high'){
        cacao_report[['callability']][['verbose']][[c]][['high']] <- paste0('Sequencing depth above ',cacao_report[['callability']][['levels']][[c]][3])
      }
    }
  }

  for(c in c('hereditary','somatic_actionable','somatic_hotspot')){
    cacao_report[['loci']][[c]] <- list()
    for(class in c('all','callable','no_coverage','low','high')){
      cacao_report[['loci']][[c]][[class]] <- data.frame()
    }
    cacao_report[['global_distribution']][[c]] <- data.frame()
    cacao_report[['eval']][[c]] <- FALSE
  }

  if(mode == 'hereditary' || mode == 'any'){
    cacao_report[['eval']][['hereditary']] <- TRUE
  }
  if(mode == 'somatic' || mode == 'any'){
    cacao_report[['eval']][['somatic_actionable']] <- TRUE
    cacao_report[['eval']][['somatic_hotspot']] <- TRUE
  }

  return(cacao_report)
}

#' Function that assigns a categorical callability level to a numeric coverage estimate
#' @param coverage_df data frame with mosdepth-annotated loci (coverage at pathogenic/actionable cancer loci)
#' @param coverage_levels array of three integers denoting thresholds for LOW_COVERAGE, CALLABLE, and HIGH_COVERAGE
#'
#' @return coverage_df
#'
assign_callability <- function(coverage_df, coverage_levels){

  coverage_df$CALLABILITY <- NA

  if("COVERAGE" %in% colnames(coverage_df)){
    if(nrow(coverage_df[coverage_df$COVERAGE == 0,]) > 0){
      coverage_df[coverage_df$COVERAGE == 0,]$CALLABILITY <- 'NO_COVERAGE'
    }
    if(nrow(coverage_df[coverage_df$COVERAGE > 0 & coverage_df$COVERAGE < coverage_levels[2],]) > 0){
      coverage_df[coverage_df$COVERAGE > 0 & coverage_df$COVERAGE < coverage_levels[2],]$CALLABILITY <- 'LOW_COVERAGE'
    }
    if(nrow(coverage_df[coverage_df$COVERAGE >= coverage_levels[2] & coverage_df$COVERAGE < coverage_levels[3],]) > 0){
      coverage_df[coverage_df$COVERAGE >= coverage_levels[2] & coverage_df$COVERAGE < coverage_levels[3],]$CALLABILITY <- 'CALLABLE'
    }
    if(nrow(coverage_df[coverage_df$COVERAGE >= coverage_levels[3],]) > 0){
      coverage_df[coverage_df$COVERAGE >= coverage_levels[3],]$CALLABILITY <- 'HIGH_COVERAGE'
    }
  }
  return(coverage_df)

}
#' Function that appends variant/region annotations to coverage tracks
#' @param coverage_df data frame with mosdepth-annotated loci (coverage at pathogenic/actionable cancernloci)
#' @param annotation_df data frame with variant/region annotations
#'
#' @return coverage_df
#'
append_annotations <- function(coverage_df = NULL, annotation_df = NULL, mode = "hereditary"){

  if(mode == "hereditary" && !is.null(coverage_df) && !is.null(annotation_df)){
    annotation_df_status <- 'OK'
    for(var in c('name','symbol','class','refseq_mrna','ensembl_transcript_id','hgvs_c', 'codon',
                 'pathogenic_loci_trait', 'pathogenic_loci_trait_link','genomic_change')){
      if(!(var %in% colnames(annotation_df))){
        annotation_df_status <- 'MISSING_DATA'
      }
    }
    if(annotation_df_status == 'OK'){
      coverage_df <- coverage_df %>%
        dplyr::left_join(dplyr::select(annotation_df, name, symbol, class,refseq_mrna, ensembl_transcript_id,
                                       hgvs_c, codon, pathogenic_loci_trait, pathogenic_loci_trait_link, genomic_change),by=c("name")) %>%
        dplyr::rename(PHENOTYPE = pathogenic_loci_trait, CLINVAR = pathogenic_loci_trait_link, SYMBOL = symbol, LOCUSTYPE = class,
                      COVERAGE = coverage, HGVSc = hgvs_c, AMINO_ACID_POSITION = codon, NAME = name, REFSEQ_TRANSCRIPT_ID = refseq_mrna,
                      ENSEMBL_TRANSCRIPT_ID = ensembl_transcript_id, GENOMIC_CHANGE = genomic_change) %>%
        dplyr::mutate(REGION = paste(paste(chrom,start,sep=":"),end,sep="-"), COVERAGE = floor(COVERAGE)) %>%
        dplyr::select(-c(chrom,start,end)) %>%
        dplyr::filter(!is.na(COVERAGE))
    }
  }

  if(mode == "somatic_actionable" && !is.null(coverage_df) && !is.null(annotation_df)){
    annotation_df_status <- 'OK'
    for(var in c('name','symbol','ensembl_transcript_id','hgvsc', 'codon','pubmed_html_link', 'evidence_id','disease_ontology_id',
                 'evidence_type','evidence_level','therapeutic_context','clinical_significance','cancer_type','genomic_change')){
      if(!(var %in% colnames(annotation_df))){
        annotation_df_status <- 'MISSING_DATA'
      }
    }
    if(annotation_df_status == 'OK'){
      coverage_df <- as.data.frame(coverage_df %>%
        dplyr::left_join(dplyr::select(annotation_df, name, symbol, ensembl_transcript_id,pubmed_html_link,disease_ontology_id,clinical_significance,
                                       hgvsc, codon, evidence_level, evidence_id, evidence_type, cancer_type, therapeutic_context, genomic_change),by=c("name")) %>%
        dplyr::rename(EVIDENCE_LEVEL = evidence_level, EVIDENCE_TYPE = evidence_type, THERAPEUTIC_CONTEXT = therapeutic_context, SYMBOL = symbol,
                      COVERAGE = coverage, HGVSc = hgvsc, AMINO_ACID_POSITION = codon, NAME = name, CLINICAL_SIGNIFICANCE = clinical_significance,
                      CANCERTYPE = cancer_type, ENSEMBL_TRANSCRIPT_ID = ensembl_transcript_id, CITATION = pubmed_html_link, GENOMIC_CHANGE = genomic_change,
                      EVIDENCE_ID = evidence_id) %>%
        dplyr::mutate(REGION = paste(paste(chrom,start,sep=":"),end,sep="-"), COVERAGE = floor(COVERAGE)) %>%
        dplyr::select(-c(chrom,start,end)) %>%
        dplyr::filter(!is.na(COVERAGE)) %>%
        dplyr::mutate(NAME = stringr::str_replace(NAME,"^EID[0-9]{1,}:","")) %>%
        dplyr::group_by(NAME, SYMBOL, REGION, ENSEMBL_TRANSCRIPT_ID, AMINO_ACID_POSITION, COVERAGE, EVIDENCE_TYPE, EVIDENCE_LEVEL, CANCERTYPE) %>%
        dplyr::summarise(THERAPEUTIC_CONTEXT = paste(unique(THERAPEUTIC_CONTEXT), collapse=", "),
                         CLINICAL_SIGNIFICANCE = paste(unique(CLINICAL_SIGNIFICANCE),collapse=", "),
                         GENOMIC_CHANGE = paste(unique(GENOMIC_CHANGE),collapse=", "),
                         CITATION = paste(unique(CITATION), collapse=", ")))


        therapeutic_contexts <- dplyr::select(coverage_df, NAME, REGION, GENOMIC_CHANGE, CANCERTYPE, THERAPEUTIC_CONTEXT) %>%
          tidyr::separate_rows(THERAPEUTIC_CONTEXT,sep=", ") %>%
          dplyr::distinct()

    }
  }

  if(mode == "somatic_hotspot" && !is.null(coverage_df) && !is.null(annotation_df)){
    annotation_df_status <- 'OK'
    for(var in c('name','symbol','hgvsp', 'pvalue',
                 'cancer_type','genomic_change')){
      if(!(var %in% colnames(annotation_df))){
        annotation_df_status <- 'MISSING_DATA'
      }
    }
    if(annotation_df_status == 'OK'){
      coverage_df <- coverage_df %>%
        dplyr::left_join(dplyr::select(annotation_df, name, symbol,
                                       hgvsp, cancer_type, pvalue, genomic_change),by=c("name")) %>%
        dplyr::rename(SYMBOL = symbol, COVERAGE = coverage, NAME = name, P_VALUE = pvalue,
                      CANCERTYPE = cancer_type, GENOMIC_CHANGE = genomic_change, HGVSp = hgvsp) %>%
        dplyr::mutate(REGION = paste(paste(chrom,start,sep=":"),end,sep="-"), COVERAGE = floor(COVERAGE)) %>%
        dplyr::select(-c(chrom,start,end)) %>%
        dplyr::filter(!is.na(COVERAGE))
    }
  }


  return(coverage_df)
}

#' Function that makes a summary of callability levels for a given track (hereditary, somatic_actionable, somatic_hotspot)
#' @param coverage_df data frame with mosdepth-annotated loci (coverage at pathogenic/actionable cancernloci)
#'
#' @return coverage_df
#'
get_global_distribution <- function(coverage_df = NULL){

  coverage_dist <- NULL
  if(!is.null(coverage_df)){
    coverage_df <- coverage_df %>%
      dplyr::select(REGION, CALLABILITY) %>%
      dplyr::distinct()

    coverage_df$n_total <- nrow(coverage_df)
    coverage_dist <- as.data.frame(
      coverage_df %>%
        dplyr::group_by(CALLABILITY, n_total) %>%
        dplyr::summarise(n = n()) %>%
        dplyr::mutate(PERCENT = round(as.numeric((n / n_total) * 100),1))
    )
    coverage_dist$CALLABILITY <-  factor(coverage_dist$CALLABILITY, levels = c('NO_COVERAGE','LOW_COVERAGE','CALLABLE','HIGH_COVERAGE'), ordered = T)
  }
  return(coverage_dist)
}


