#!/usr/bin/env Rscript

options(warn=-1)
args = commandArgs(trailingOnly=TRUE)

suppressWarnings(suppressPackageStartupMessages(library(magrittr)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(rmarkdown)))
suppressWarnings(suppressPackageStartupMessages(library(tidyr)))
suppressWarnings(suppressPackageStartupMessages(library(rlogging)))
suppressWarnings(suppressPackageStartupMessages(library(cacao)))

tsv_annotation_fname <- list()

bed_coverage_fname <- as.character(args[1])
host_name_alignment <- as.character(args[2])
host_name_target <- as.character(args[3])
tsv_annotation_fname[['hereditary']] <- as.character(args[4])
tsv_annotation_fname[['somatic_actionable']] <- as.character(args[5])
tsv_annotation_fname[['somatic_hotspot']] <- as.character(args[6])
sample_name <- as.character(args[7])
cacao_mode <- as.character(args[8])
coverage_levels_germline <- as.integer(stringr::str_split(as.character(args[9]),":")[[1]])
coverage_levels_somatic <- as.integer(stringr::str_split(as.character(args[10]),":")[[1]])
mapq_threshold <- as.integer(args[11])
genome_assembly <- as.character(args[12])
cacao_version <- as.character(args[13])
out_dir <- as.character(args[14])

fname_json <- paste0(sample_name,"_", genome_assembly,"_coverage_cacao.json")
fname_html <- paste0(sample_name,"_", genome_assembly,"_coverage_cacao.html")

rlogging::SetTimeStampFormat(ts.format="%Y-%m-%d %H:%M:%S ")
rlogging::SetLogFile(NULL)

coverage_pr_loci <- list()
loci_annotations <- list()
cacao_report <- cacao::init_report(cacao_mode, host_name_alignment, host_name_target, sample_name, cacao_version, genome_assembly, mapq_threshold, coverage_levels_germline, coverage_levels_somatic)

rlogging::message('------')
rlogging::message("Reading raw coverage output file from mosdepth")
if(bed_coverage_fname != "NA" && file.exists(bed_coverage_fname)){
  coverage_pr_loci <- cacao::read_raw_coverage(bed_coverage_fname)
}
rlogging::message('------')

for(c in c('hereditary','somatic_actionable','somatic_hotspot')){
  if(tsv_annotation_fname[[c]] != "NA" && file.exists(tsv_annotation_fname[[c]])){
    rlogging::message("Assigning callability levels across loci and appending annotations - track type:",c)
    rlogging::message('------')
    
    loci_annotations[[c]] <- read.table(file=tsv_annotation_fname[[c]],header=T,quote="",comment.char="",sep="\t",stringsAsFactors = F)
    
    coverage_pr_loci[[c]] <- cacao::append_annotations(coverage_df = coverage_pr_loci[[c]], annotation_df = loci_annotations[[c]], mode = c)
    if(c == 'hereditary'){
      if(nrow(coverage_pr_loci[[c]]) > 0){
        coverage_pr_loci[[c]] <- cacao::assign_callability(coverage_pr_loci[[c]], coverage_levels_germline)
        coverage_pr_loci[[c]] <- coverage_pr_loci[[c]] %>% dplyr::select(-TARGET) %>% dplyr::select(SYMBOL, HGVSc, CALLABILITY, COVERAGE, LOCUSTYPE, AMINO_ACID_POSITION, dplyr::everything())
      }
    }
    if(c == 'somatic_actionable'){
      if(nrow(coverage_pr_loci[[c]]) > 0){
        coverage_pr_loci[[c]] <- cacao::assign_callability(coverage_pr_loci[[c]], coverage_levels_somatic)
        coverage_pr_loci[[c]] <- coverage_pr_loci[[c]] %>% 
          dplyr::select(-TARGET) %>% 
          dplyr::select(NAME, CALLABILITY, CANCERTYPE, CLINICAL_SIGNIFICANCE, EVIDENCE_LEVEL, dplyr::everything())
      }
    }
    
    if(c == 'somatic_hotspot'){
      if(nrow(coverage_pr_loci[[c]]) > 0){
        coverage_pr_loci[[c]] <- cacao::assign_callability(coverage_pr_loci[[c]], coverage_levels_somatic)
        coverage_pr_loci[[c]] <- coverage_pr_loci[[c]]  %>% dplyr::select(-TARGET) %>% dplyr::select(SYMBOL, NAME, CALLABILITY, CANCERTYPE, COVERAGE, dplyr::everything())
      }
    }
    cacao_report[['loci']][[c]][['all']] <- coverage_pr_loci[[c]]
    if(nrow(coverage_pr_loci[[c]]) > 0 && "CALLABILITY" %in% colnames(coverage_pr_loci[[c]])){
      cacao_report[['loci']][[c]][['no_coverage']] <- dplyr::filter(coverage_pr_loci[[c]], CALLABILITY == "NO_COVERAGE")
      cacao_report[['loci']][[c]][['low']] <- dplyr::filter(coverage_pr_loci[[c]], CALLABILITY == "LOW_COVERAGE")
      cacao_report[['loci']][[c]][['callable']] <- dplyr::filter(coverage_pr_loci[[c]], CALLABILITY == "CALLABLE") 
      cacao_report[['loci']][[c]][['high']] <- dplyr::filter(coverage_pr_loci[[c]], CALLABILITY == "HIGH_COVERAGE")
      cacao_report[['coverage_distribution']][[c]][['global']][['data']] <- cacao::global_coverage_distribution(cacao_report[['loci']][[c]][['all']])
      cacao_report[['coverage_distribution']][[c]][['gene']][['data']] <- cacao::coverage_distribution_pr_gene(cacao_report[['loci']][[c]][['all']])
      cacao_report[['coverage_distribution']][[c]][['gene']][['plot']] <- cacao::plot_gene_distribution(cacao_report[['coverage_distribution']][[c]][['gene']][['data']])
      cacao_report[['coverage_distribution']][[c]][['global']][['plot']] <- cacao::plot_global_distribution(cacao_report[['coverage_distribution']][[c]][['global']][['data']])
    }
  }
}
rlogging::message("Writing JSON file with CACAO report contents")
cacao_content_json <- cacao_report
for(c in c('hereditary','somatic_actionable','somatic_hotspot')){
  cacao_content_json[['coverage_distribution']][[c]][['global']][['plot']] <- NULL
  cacao_content_json[['coverage_distribution']][[c]][['gene']][['plot']] <- NULL
}
cacao_json <- jsonlite::toJSON(cacao_content_json, pretty=T,na='string',null = 'null',force=T)
write(cacao_json, fname_json)
system(paste0('rm -f ',sample_name,'.mosdepth*dist* '), intern = F)

rlogging::message('------')
rlogging::message("Rendering HTML report with rmarkdown")

rmarkdown::render(system.file("templates","cacao_report.Rmd", package="cacao"), output_format = rmarkdown::html_document(theme = "default", toc = T, toc_depth = 3, toc_float = T, number_sections = F), output_file = fname_html, output_dir = out_dir, clean = T, intermediates_dir = out_dir, quiet = T)


# m <- cacao::coverage_distribution_pr_gene(coverage_pr_loci, mode = 'hereditary')
# n_loci <- m %>% dplyr::select(SYMBOL,tot_loci,tot_frac) %>% dplyr::distinct()
# 
# cols <- c("NO_COVERAGE" = "#FC4E2A", "LOW_COVERAGE" = "#FD8D3C", "CALLABLE" = "#78C679", "HIGH_COVERAGE" = "#41AB5D")
# cols <- c("NO_COVERAGE" = "#FC4E2A", "LOW_COVERAGE" = "#FD8D3C", "CALLABLE" = "#78C679", "HIGH_COVERAGE" = "#207733")
# 
# p <- ggplot(tail(m,60),aes(x=SYMBOL, y=frac, fill=CALLABILITY)) +
#   labs (y = "Percentage of pathogenic loci") +
#   geom_bar(position = position_stack(), stat="identity") +
#   geom_text(aes(x = SYMBOL, y = tot_frac + 4, label = paste0("N = ",tot_loci)), size = 4) +
#   coord_flip() +
#   ggplot2::theme_classic() +
#   ggplot2::theme(axis.text.x = ggplot2::element_text(family = "Helvetica", size = 12, vjust = -0.1),
#                  axis.title.x = ggplot2::element_text(family = "Helvetica", size = 12, vjust = -2),
#                  axis.title.y = ggplot2::element_blank(),
#                  legend.title = ggplot2::element_blank(),
#                  axis.text.y = ggplot2::element_text(family = "Helvetica", size = 12, hjust = 0.5),
#                  plot.margin = (grid::unit(c(0.5, 1, 1, 0.5), "cm")),
#                  legend.text = ggplot2::element_text(family = "Helvetica", size = 12))
# p <- p + ggplot2::scale_fill_manual(values = cols)
# p <- plotly::ggplotly(p)
# p


# tot_loci <- sum(cacao_report[['global_distribution']][['somatic_actionable']]$n)
# 
# p <- ggplot2::ggplot(ggplot2::aes(x="", y=PERCENT, fill = CALLABILITY), data = cacao_report[['global_distribution']][['somatic_actionable']]) +
#   ggplot2::geom_bar(stat = 'identity', position = ggplot2::position_stack()) +
#   ggplot2::coord_flip() +
#   ggplot2::theme_classic() +
#   ggplot2::geom_text(ggplot2::aes(y = 104), label = paste0("N = ",tot_loci), size = 4) +
#   ggplot2::scale_y_continuous("Percent of all variant loci",breaks=seq(0,100,by=10),labels=seq(0,100,by=10)) + 
#   ggplot2::theme(legend.title = ggplot2::element_blank(),
#                  axis.text.x = ggplot2::element_text(family = "Helvetica", size = 12, vjust = -0.1),
#                  axis.title.x = ggplot2::element_text(family = "Helvetica", size = 12, vjust = -2.5),
#                  axis.title.y = ggplot2::element_blank(),
#                  #axis.text.y = ggplot2::element_text(family = "Helvetica", size = 12, angle = -90, hjust = 0.5),
#                  plot.margin = (grid::unit(c(0.5, 1, 1, 0.5), "cm")),
#                  legend.text = ggplot2::element_text(family = "Helvetica", size = 12))
# 
# p <- p + ggplot2::scale_fill_manual(values = color_callability_map)
# p
# #plotly::ggplotly(p)
# htmltools::br()


