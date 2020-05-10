#' Function that initiates the cacao reporting object
#' @param mode mode for reporting (somatic, hereditary, any)
#' @param fname_aln host filename for the alignment (BAM/CRAM)
#' @param fname_target host filename for the target regions (BED)
#' @param sample_name name of sample
#' @param version cacao version
#' @param genome_assembly human genome assembly (grch37/grch38)
#' @param mapq_threshold alignment quality (MAPQ) filter
#' @param callability_levels_germline array of three integers denoting thresholds for LOW_COVERAGE, CALLABLE, and HIGH_COVERAGE (germline setting)
#' @param callability_levels_somatic array of three integers denoting thresholds for LOW_COVERAGE, CALLABLE, and HIGH_COVERAGE (somatic setting)
#'
#' @return coverage_df
#'
init_report <- function(mode, fname_aln, fname_target, sample_name, version, genome_assembly, mapq_threshold, callability_levels_germline, callability_levels_somatic){

  cacao_report <- list()
  cacao_report[['sample_name']] <- sample_name
  cacao_report[['version']] <- version
  cacao_report[['mode']] <- mode
  cacao_report[['genome_assembly']] <- genome_assembly
  cacao_report[['mapq']] <- mapq_threshold
  cacao_report[['host_alignment_fname']] <- fname_aln
  cacao_report[['host_target_fname']] <- fname_target
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
    cacao_report[['coverage_distribution']][[c]] <- list()
    for(p in c('global','gene','chromosome')){
      cacao_report[['coverage_distribution']][[c]][[p]] <- list()
      cacao_report[['coverage_distribution']][[c]][[p]][['data']] <- data.frame()
      cacao_report[['coverage_distribution']][[c]][[p]][['plot']] <- NULL
    }

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
    coverage_df <- coverage_df %>%
      dplyr::mutate(CALLABILITY = dplyr::if_else(TARGET == F,"NON_TARGET",as.character(CALLABILITY))) %>%
      dplyr::mutate(CALLABILITY = dplyr::if_else(COVERAGE == 0 & TARGET == T,"NO_COVERAGE",as.character(CALLABILITY))) %>%
      dplyr::mutate(CALLABILITY = dplyr::if_else(COVERAGE > 0 & TARGET == T & COVERAGE < coverage_levels[2],"LOW_COVERAGE",as.character(CALLABILITY))) %>%
      dplyr::mutate(CALLABILITY = dplyr::if_else(TARGET == T & COVERAGE >= coverage_levels[2] & COVERAGE < coverage_levels[3],"CALLABLE",as.character(CALLABILITY))) %>%
      dplyr::mutate(CALLABILITY = dplyr::if_else(TARGET == T & COVERAGE >= coverage_levels[3],"HIGH_COVERAGE",as.character(CALLABILITY)))
  }
  return(coverage_df)

}

#' Function that reads raw mosdepth coverage file
#' @param bed_coverage_fname coverage file (BED) from mosdepth + bedtools annotate
#' #'
#' @return coverage_pr_loci
#'
read_raw_coverage <- function(bed_coverage_fname){
  coverage_all_loci <- read.table(file=bed_coverage_fname,header=F,quote="",comment.char="",sep="\t",stringsAsFactors = F)
  if(ncol(coverage_all_loci) == 5){
    colnames(coverage_all_loci) <- c('CHROM','START','END','NAME','COVERAGE')
    coverage_all_loci$TARGET <- TRUE
  }
  else if(ncol(coverage_all_loci) == 6){
    colnames(coverage_all_loci) <- c('CHROM','START','END','NAME','COVERAGE','TARGET_OVERLAP')
    coverage_all_loci <- coverage_all_loci %>%
      dplyr::mutate(TARGET = dplyr::if_else(TARGET_OVERLAP == 0,FALSE,TRUE)) %>%
      dplyr::select(-TARGET_OVERLAP)
  }

  coverage_pr_loci <- list()
  coverage_pr_loci[['all']] <- coverage_all_loci
  coverage_pr_loci[['hereditary']] <- dplyr::filter(coverage_all_loci, stringr::str_detect(NAME,"^clinvar_path:"))
  coverage_pr_loci[['somatic_actionable']] <- dplyr::filter(coverage_all_loci, stringr::str_detect(NAME,"^civic:"))
  coverage_pr_loci[['somatic_hotspot']] <- dplyr::filter(coverage_all_loci, stringr::str_detect(NAME,"^hotspot:"))
  return(coverage_pr_loci)
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
        dplyr::left_join(dplyr::select(annotation_df, name, symbol, class, refseq_mrna, ensembl_transcript_id,
                                       hgvs_c, codon, pathogenic_loci_trait, pathogenic_loci_trait_link, genomic_change),by = c("NAME" = "name")) %>%
        dplyr::rename(PHENOTYPE = pathogenic_loci_trait, CLINVAR = pathogenic_loci_trait_link, SYMBOL = symbol, LOCUSTYPE = class,
                      HGVSc = hgvs_c, AMINO_ACID_POSITION = codon, REFSEQ_TRANSCRIPT_ID = refseq_mrna,
                      ENSEMBL_TRANSCRIPT_ID = ensembl_transcript_id, GENOMIC_CHANGE = genomic_change) %>%
        dplyr::mutate(REGION = paste(paste(CHROM,START,sep=":"),END,sep="-"), COVERAGE = floor(COVERAGE)) %>%
        dplyr::mutate(NAME = stringr::str_replace_all(NAME,"clinvar_path:","")) %>%
        dplyr::mutate(ENSEMBL_TRANSCRIPT_ID = paste0('<a href=\'http://www.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;t=',ENSEMBL_TRANSCRIPT_ID,'\' target=\'_blank\'>',ENSEMBL_TRANSCRIPT_ID,'</a>')) %>%
        dplyr::mutate(REFSEQ_TRANSCRIPT_ID = paste0('<a href=\'http://www.ncbi.nlm.nih.gov/nuccore/',REFSEQ_TRANSCRIPT_ID,'\' target=\'_blank\'>',REFSEQ_TRANSCRIPT_ID,'</a>')) %>%
        dplyr::select(-c(CHROM,START,END)) %>%
        dplyr::filter(!is.na(COVERAGE))
    }
  }

  if(mode == "somatic_actionable" && !is.null(coverage_df) && !is.null(annotation_df)){
    annotation_df_status <- 'OK'
    for(var in c('name','symbol','ensembl_transcript_id','refseq_mrna','hgvsc', 'codon','citation', 'evidence_id','disease_ontology_id',
                 'evidence_type','evidence_level','therapeutic_context','clinical_significance','cancer_type','genomic_change')){
      if(!(var %in% colnames(annotation_df))){
        annotation_df_status <- 'MISSING_DATA'
      }
    }
    if(annotation_df_status == 'OK'){
      coverage_df <- as.data.frame(coverage_df %>%
          dplyr::left_join(dplyr::select(annotation_df, name, symbol, ensembl_transcript_id,
                                         refseq_mrna, citation, disease_ontology_id,
                                         clinical_significance, hgvsc, codon, evidence_level,
                                         evidence_id, evidence_type, cancer_type, therapeutic_context,
                                         genomic_change), by=c("NAME" = "name")) %>%
          dplyr::rename(EVIDENCE_LEVEL = evidence_level, EVIDENCE_TYPE = evidence_type, THERAPEUTIC_CONTEXT = therapeutic_context, SYMBOL = symbol,
                      HGVSc = hgvsc, AMINO_ACID_POSITION = codon, CLINICAL_SIGNIFICANCE = clinical_significance, REFSEQ_TRANSCRIPT_ID = refseq_mrna,
                      CANCERTYPE = cancer_type, ENSEMBL_TRANSCRIPT_ID = ensembl_transcript_id, CITATION = citation, GENOMIC_CHANGE = genomic_change,
                      EVIDENCE_ID = evidence_id) %>%
          dplyr::mutate(REGION = paste(paste(CHROM,START,sep=":"),END,sep="-"), COVERAGE = floor(COVERAGE)) %>%
          dplyr::mutate(ENSEMBL_TRANSCRIPT_ID = paste0('<a href=\'http://www.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;t=',ENSEMBL_TRANSCRIPT_ID,'\' target=\'_blank\'>',ENSEMBL_TRANSCRIPT_ID,'</a>')) %>%
          dplyr::mutate(REFSEQ_TRANSCRIPT_ID = paste0('<a href=\'http://www.ncbi.nlm.nih.gov/nuccore/',REFSEQ_TRANSCRIPT_ID,'\' target=\'_blank\'>',REFSEQ_TRANSCRIPT_ID,'</a>')) %>%

          dplyr::select(-c(CHROM,START,END)) %>%
          dplyr::filter(!is.na(COVERAGE)) %>%
          dplyr::mutate(NAME = stringr::str_replace_all(NAME,"civic:EID[0-9]{1,}:","")) %>%
          dplyr::group_by(NAME, SYMBOL, REGION, ENSEMBL_TRANSCRIPT_ID, REFSEQ_TRANSCRIPT_ID, AMINO_ACID_POSITION, COVERAGE, TARGET, EVIDENCE_TYPE, EVIDENCE_LEVEL, CANCERTYPE) %>%
          dplyr::summarise(THERAPEUTIC_CONTEXT = paste(unique(THERAPEUTIC_CONTEXT), collapse=", "),
                         CLINICAL_SIGNIFICANCE = paste(unique(CLINICAL_SIGNIFICANCE),collapse=", "),
                         GENOMIC_CHANGE = paste(unique(GENOMIC_CHANGE),collapse=", "),
                         CITATION = paste(unique(CITATION), collapse=", ")))


        therapeutic_contexts <- dplyr::select(coverage_df, NAME, REGION, TARGET, GENOMIC_CHANGE, CANCERTYPE, THERAPEUTIC_CONTEXT) %>%
          tidyr::separate_rows(THERAPEUTIC_CONTEXT,sep=", ") %>%
          dplyr::distinct()

    }
  }

  if(mode == "somatic_hotspot" && !is.null(coverage_df) && !is.null(annotation_df)){
    annotation_df_status <- 'OK'
    for(var in c('name','symbol','hgvsp', 'pvalue','cancer_type','genomic_change')){
      if(!(var %in% colnames(annotation_df))){
        annotation_df_status <- 'MISSING_DATA'
      }
    }
    if(annotation_df_status == 'OK'){
      coverage_df <- coverage_df %>%
        dplyr::left_join(dplyr::select(annotation_df, name, symbol,
                                       hgvsp, cancer_type, pvalue, genomic_change),by=c("NAME" = "name")) %>%
        dplyr::rename(SYMBOL = symbol, P_VALUE = pvalue,
                      CANCERTYPE = cancer_type, GENOMIC_CHANGE = genomic_change, HGVSp = hgvsp) %>%
        dplyr::mutate(REGION = paste(paste(CHROM,START,sep=":"),END,sep="-"), COVERAGE = floor(COVERAGE)) %>%
        dplyr::mutate(NAME = stringr::str_replace_all(NAME,"hotspot:","")) %>%
        dplyr::select(-c(CHROM,START,END)) %>%
        dplyr::filter(!is.na(COVERAGE))
    }
  }


  return(coverage_df)
}

#' Function that makes a summary of callability levels for a given track (hereditary, somatic_actionable, somatic_hotspot)
#' @param coverage_pr_loci data frame with mosdepth-annotated loci (coverage at pathogenic/actionable cancernloci)
#'
#' @return coverage_pr_loci
#'
global_coverage_distribution <- function(coverage_pr_loci = NULL){

  coverage_dist <- data.frame()
  if(!is.null(coverage_pr_loci) & nrow(coverage_pr_loci) > 0){
    coverage_pr_loci <- coverage_pr_loci %>%
      dplyr::select(REGION, CALLABILITY) %>%
      dplyr::distinct()

    coverage_pr_loci$n_total <- nrow(coverage_pr_loci)
    coverage_dist <- as.data.frame(
      coverage_pr_loci %>%
        dplyr::group_by(CALLABILITY, n_total) %>%
        dplyr::summarise(n = n()) %>%
        dplyr::mutate(PERCENT = round(as.numeric((n / n_total) * 100),1))
    )
    coverage_dist$CALLABILITY <-  factor(coverage_dist$CALLABILITY, levels = c('NON_TARGET','NO_COVERAGE','LOW_COVERAGE','CALLABLE','HIGH_COVERAGE'), ordered = T)
  }
  return(coverage_dist)
}


#' Function that makes a summary of callability levels for a given track (hereditary, somatic_actionable, somatic_hotspot) pr. gene
#' @param coverage_pr_loci data frame with mosdepth-annotated loci (coverage at pathogenic/actionable cancer loci)
#'
#' @return c
#'
coverage_distribution_pr_gene <- function(coverage_pr_loci = NULL){

  callability_pr_gene <- as.data.frame(
    coverage_pr_loci %>%
      dplyr::select(SYMBOL,REGION,CALLABILITY) %>%
      dplyr::distinct() %>%
      dplyr::group_by(SYMBOL) %>%
      dplyr::summarise(n_loci = n()))

  b <- as.data.frame(
    coverage_pr_loci %>%
      dplyr::select(SYMBOL,REGION,CALLABILITY) %>%
      dplyr::distinct() %>%
      dplyr::group_by(SYMBOL, CALLABILITY) %>%
      dplyr::summarise(n = n()))

  c <- dplyr::left_join(b,callability_pr_gene,by=c("SYMBOL")) %>%
    dplyr::mutate(frac = round((n / n_loci) * 100,3))
  d <- as.data.frame(c %>% dplyr::group_by(SYMBOL) %>% dplyr::summarise(tot_loci = sum(n), tot_frac = sum(frac)))
  c <- dplyr::left_join(c,d,by=c("SYMBOL"))

  ## RANK GENES FROM FRACTION OF LOCI WITH NO COVERAGE TO GENES WITH ALL LOCI OF HIGH COVERAGE
  c$rank <- 0
  if(nrow(c[c$CALLABILITY == "NO_COVERAGE",]) > 0){
    c[c$CALLABILITY == "NO_COVERAGE",]$rank <- c[c$CALLABILITY == "NO_COVERAGE",]$frac * -2
  }
  if(nrow(c[c$CALLABILITY == "LOW_COVERAGE",]) > 0){
    c[c$CALLABILITY == "LOW_COVERAGE",]$rank <- c[c$CALLABILITY == "LOW_COVERAGE",]$frac * -1
  }
  if(nrow(c[c$CALLABILITY == "CALLABLE",]) > 0){
    c[c$CALLABILITY == "CALLABLE",]$rank <- c[c$CALLABILITY == "CALLABLE",]$frac * 1
  }
  if(nrow(c[c$CALLABILITY == "HIGH_COVERAGE",]) > 0){
    c[c$CALLABILITY == "HIGH_COVERAGE",]$rank <- c[c$CALLABILITY == "HIGH_COVERAGE",]$frac * 2
  }
  if(nrow(c[c$CALLABILITY == "NON_TARGET",]) > 0){
    c[c$CALLABILITY == "NON_TARGET",]$rank <- c[c$CALLABILITY == "NON_TARGET",]$frac * -3
  }

  m <- c %>% dplyr::group_by(SYMBOL) %>% dplyr::summarise(rank_order = sum(rank)) %>% dplyr::arrange(rank_order)
  c <- dplyr::left_join(c, m,by=c("SYMBOL"))
  c <- c %>% dplyr::mutate(SYMBOL = forcats::fct_reorder(SYMBOL, rank_order)) %>% dplyr::arrange(rank_order)

  c$CALLABILITY <-  factor(c$CALLABILITY, levels = c('NON_TARGET','NO_COVERAGE','LOW_COVERAGE','CALLABLE','HIGH_COVERAGE'), ordered = T)

  return(c)
}


#' Function that generates a ggplot of the variant loci coverage pr. gene
#' @param gene_coverage_distribution data frame with mosdepth-annotated loci (coverage at pathogenic/actionable cancer loci)
#' @param color_callability_map named vector mapping colors to callability levels
#' @param plotly logical indicating if ggplot should be makde with plotly
#'
#' @return p
#'
plot_gene_distribution <- function(gene_coverage_distribution, color_callability_map = NULL, plotly = FALSE){

  if(is.null(color_callability_map)){
    color_callability_map <- c("NON_TARGET" = "#E0DFDF", "NO_COVERAGE" = "#FC4E2A", "LOW_COVERAGE" = "#FD8D3C", "CALLABLE" = "#78C679", "HIGH_COVERAGE" = "#207733")
  }
  p <- ggplot2::ggplot(gene_coverage_distribution,ggplot2::aes(x=SYMBOL, y=frac, fill=CALLABILITY)) +
    ggplot2::geom_bar(position = ggplot2::position_stack(), stat="identity") +
    ggplot2::geom_text(ggplot2::aes(x = SYMBOL, y = tot_frac + 5, label = paste0("N = ",tot_loci)), size = 4) +
    ggplot2::scale_y_continuous("Percent of all variant loci",breaks=seq(0,100,by=10),labels=seq(0,100,by=10)) +
    ggplot2::coord_flip() +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(family = "Helvetica", size = 12, vjust = -0.1),
                   axis.title.x = ggplot2::element_text(family = "Helvetica", size = 13, vjust = -2.5),
                   axis.title.y = ggplot2::element_blank(),
                   legend.title = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_text(family = "Helvetica", size = 12, hjust = 0.5),
                   plot.margin = (grid::unit(c(0.5, 1, 1, 0.5), "cm")),
                   legend.text = ggplot2::element_text(family = "Helvetica", size = 12))
  p <- p + ggplot2::scale_fill_manual(values = color_callability_map)
  if(plotly == TRUE){
    p <- plotly::ggplotly(p)
  }
  return(p)

}

#' Function that generates a ggplot of the variant loci coverage pr. gene
#' @param gene_coverage_distribution data frame with mosdepth-annotated loci (coverage at pathogenic/actionable cancer loci)
#' @param color_callability_map named vector mapping colors to callability levels
#' @param plotly logical indicating if ggplot should be makde with plotly
#'
#' @return p
#'
plot_global_distribution <- function(global_coverage_distribution, color_callability_map = NULL, plotly = FALSE){

  if(is.null(color_callability_map)){
    color_callability_map <- c("NON_TARGET" = "#E0DFDF", "NO_COVERAGE" = "#FC4E2A", "LOW_COVERAGE" = "#FD8D3C", "CALLABLE" = "#78C679", "HIGH_COVERAGE" = "#207733")
  }
  tot_loci <- sum(global_coverage_distribution$n)

  p <- ggplot2::ggplot(ggplot2::aes(x="CALLABILITY", y=PERCENT, fill = CALLABILITY), data = global_coverage_distribution) +
    ggplot2::geom_bar(stat = 'identity', position = ggplot2::position_stack()) +
    ggplot2::coord_flip() +
    ggplot2::theme_classic() +
    ggplot2::geom_text(ggplot2::aes(y = 106), label = paste0("N = ",tot_loci), size = 4) +
    ggplot2::scale_y_continuous("Percent of all variant loci",breaks=seq(0,100,by=10),labels=seq(0,100,by=10)) +
    ggplot2::theme(legend.title = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(family = "Helvetica", size = 12, vjust = -0.1),
                   axis.title.x = ggplot2::element_text(family = "Helvetica", size = 13, vjust = -2.5),
                   axis.title.y = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_text(family = "Helvetica", size = 12, angle = -90, hjust = 0.5),
                   plot.margin = (grid::unit(c(0.5, 1, 1, 0.5), "cm")),
                   legend.text = ggplot2::element_text(family = "Helvetica", size = 12))

  p <- p + ggplot2::scale_fill_manual(values = color_callability_map)
  if(plotly == TRUE){
    p <- plotly::ggplotly(p)
  }
  return(p)

}

#' Function that strips HTML links (keeps the link text) from a data frame
#' @param coverage_pr_loci data frame with mosdepth-annotated loci (coverage at pathogenic/actionable cancer loci)
#'
#' @return coverage_pr_loci
#'

strip_html_links <- function(coverage_pr_loci){
  if('ENSEMBL_TRANSCRIPT_ID' %in% colnames(coverage_pr_loci)){
    coverage_pr_loci$ENSEMBL_TRANSCRIPT_ID <- unlist(lapply(stringr::str_match_all(coverage_pr_loci$ENSEMBL_TRANSCRIPT_ID,">.+<"),paste,collapse=","))
    coverage_pr_loci$ENSEMBL_TRANSCRIPT_ID <- stringr::str_replace_all(coverage_pr_loci$ENSEMBL_TRANSCRIPT_ID,">|<", "")
  }
  if('REFSEQ_TRANSCRIPT_ID' %in% colnames(coverage_pr_loci)){
    coverage_pr_loci$REFSEQ_TRANSCRIPT_ID <- unlist(lapply(stringr::str_match_all(coverage_pr_loci$REFSEQ_TRANSCRIPT_ID,">.+<"),paste,collapse=","))
    coverage_pr_loci$REFSEQ_TRANSCRIPT_ID <- stringr::str_replace_all(coverage_pr_loci$REFSEQ_TRANSCRIPT_ID,">|<", "")
  }
  return(coverage_pr_loci)
}
