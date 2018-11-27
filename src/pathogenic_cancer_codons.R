
library(magrittr)
medgen_map <- readRDS(file="data-raw/medgen_map.rds")

datestamp <- '20181105'

all_tracks <- list()
for(ttype in c('bed','tsv')){
  all_tracks[[ttype]] <- list()
  for(b in c('grch37','grch38')){
    all_tracks[[ttype]][[b]] <- list()
    all_tracks[[ttype]][[b]][['clinvar_path']] <- data.frame()
    all_tracks[[ttype]][[b]][['civic']] <- data.frame()
    all_tracks[[ttype]][[b]][['hotspot']] <- data.frame()
  }
}

sort_bed_regions <- function(unsorted_regions){
  sorted_regions <- NULL
  if("chrom" %in% colnames(unsorted_regions) & "start" %in% colnames(unsorted_regions) & "end" %in% colnames(unsorted_regions)){
    chrOrder <- c(as.character(c(1:22)),"X","Y")
    unsorted_regions$chrom <- factor(unsorted_regions$chrom, levels=chrOrder)
    unsorted_regions <- unsorted_regions[order(unsorted_regions$chrom),]
    
    sorted_regions <- NULL
    for(chrom in chrOrder){
      if(nrow(unsorted_regions[unsorted_regions$chrom == chrom,]) > 0){
        chrom_regions <- unsorted_regions[unsorted_regions$chrom == chrom,]
        chrom_regions_sorted <- chrom_regions[with(chrom_regions, order(start, end)),]
        sorted_regions <- rbind(sorted_regions, chrom_regions_sorted)
      }
    }
    
    sorted_regions$start <- as.integer(sorted_regions$start)
    sorted_regions$end <- as.integer(sorted_regions$end)
  }
  return(sorted_regions)
  
}

## Hereditary cancer

for(build in c('grch37','grch38')){
  gencode_xref <- read.table(gzfile(paste0("data-raw/pcgr_onco_xref.",build,".tsv.gz")),header=T,stringsAsFactors = F,quote="",sep="\t",na.strings=c("","NA"),comment.char="#")
  gencode_xref <- dplyr::select(gencode_xref, -c(refseq_mrna, strand,transcript_start,transcript_end)) %>% dplyr::distinct()
  gencode_xref <- dplyr::filter(gencode_xref, !is.na(chrom) & !is.na(start) & !is.na(end))
  gencode_xref <- dplyr::rename(gencode_xref, transcript_start = start, transcript_end = end, biotype = gencode_transcript_type)
  gencode_xref <- dplyr::filter(gencode_xref, !stringr::str_detect(chrom,"^GL"))
  
  refseq_gencode_trans <- read.table(gzfile(paste0("data-raw/pcgr_onco_xref.",build,".tsv.gz")),header=T,stringsAsFactors = F,quote="",sep="\t",na.strings=c("","NA"),comment.char="#") %>%
    dplyr::select(ensembl_transcript_id, refseq_mrna) %>%
    dplyr::filter(!is.na(refseq_mrna))
  
  cancer_predisposition <- gencode_xref %>% 
    dplyr::select(ensembl_gene_id, symbol, entrezgene, predisp_source, predisp_cancer_cui, predisp_syndrome_cui, predisp_cancer_moi) %>%
    dplyr::filter(!is.na(predisp_source)) %>%
    dplyr::distinct()

  clinvar <- read.table(gzfile(paste0("data-raw/clinvar.",build,".tsv.gz")), header=T, stringsAsFactors = F, quote="", sep="\t",na.strings=c("","-",NA),comment.char="")
  clinvar$var_id <- paste(paste0('chr',clinvar$chrom),clinvar$pos,clinvar$ref,clinvar$alt,sep="_")
  clinvar <- dplyr::left_join(clinvar, refseq_gencode_trans)
  
  ## PATHOGENIC PROTEIN LOCI IN CANCER PREDISPOSITION GENES
  clinvar_pathogenic_cpg <- dplyr::semi_join(dplyr::filter(clinvar, (origin == 'germline' | origin == 'germline,somatic') & pathogenic == 1),cancer_predisposition, by=c("symbol"))
  
  
  cpg_pathogenic_insertions <- as.data.frame(
      clinvar_pathogenic_cpg %>%
      dplyr::select(symbol, chrom,pos, ref, alt, refseq_mrna, variation_type, molecular_consequence, ensembl_transcript_id, hgvs_p, hgvs_c, variation_id, cui) %>%
      dplyr::filter(variation_type == 'Duplication' | variation_type == 'Insertion') %>%
      dplyr::filter(!stringr::str_detect(molecular_consequence,"splice_")) %>%
      dplyr::filter(!is.na(hgvs_p)) %>%
      dplyr::mutate(codon = stringr::str_replace_all(hgvs_p,"\\*$|fs\\*[0-9]{1,}$","")) %>%
      dplyr::mutate(codon = stringr::str_replace_all(codon,"[A-Z]{1,}|[a-z]{1,}|\\.","")) %>%
      dplyr::mutate(name = paste0("clinvar_path:",symbol,":",stringr::str_replace(hgvs_p,"[A-Z]fs$",""))) %>%
      dplyr::mutate(genomic_change = paste(chrom,pos,ref,alt,sep="_")) %>%
      tidyr::separate_rows(cui) %>%
      dplyr::left_join(medgen_map) %>%
      dplyr::group_by(symbol,hgvs_c,hgvs_p,codon,alt,name,chrom,pos,variation_id,genomic_change,refseq_mrna, ensembl_transcript_id) %>%
      dplyr::summarise(trait = paste(unique(cui_name),collapse=";")) %>%
      dplyr::mutate(link = paste0("<a href='https://www.ncbi.nlm.nih.gov/clinvar/variation/",variation_id,"' target='_blank'>",symbol,"|",hgvs_c,"|",hgvs_p," - ",trait,"</a>")) %>%
      dplyr::group_by(symbol,chrom,codon,refseq_mrna,ensembl_transcript_id,name) %>%
      dplyr::summarise(pathogenic_loci_trait_link = paste(link, collapse=", "), pathogenic_loci_trait = paste0(unique(trait),collapse=";"), variation_id = paste(variation_id, collapse=","), genomic_change = paste(unique(genomic_change),collapse=","), hgvs_c = paste(unique(hgvs_c), collapse=","), hgvs_p = paste(unique(hgvs_p), collapse=","), max_length_alt = max(nchar(alt)), min_pos = min(pos), max_pos = max(pos)) %>%
      dplyr::mutate(start = as.integer(min_pos) - 6, end = max_pos + nchar(max_length_alt) + 4) %>%
      dplyr::select(-c(max_pos,max_length_alt,min_pos)))
  
  cpg_pathogenic_deletions <- as.data.frame(
    clinvar_pathogenic_cpg %>%
      dplyr::select(symbol, chrom,pos, ref, alt, refseq_mrna, variation_type, molecular_consequence, ensembl_transcript_id, hgvs_p, hgvs_c, variation_id, cui) %>%
      dplyr::filter(variation_type == 'Deletion') %>%
      dplyr::filter(!stringr::str_detect(molecular_consequence,"splice_")) %>%
      dplyr::filter(!is.na(hgvs_p)) %>%
      dplyr::mutate(codon = stringr::str_replace_all(hgvs_p,"\\*$|fs\\*[0-9]{1,}$","")) %>%
      dplyr::mutate(codon = stringr::str_replace_all(codon,"[A-Z]{1,}|[a-z]{1,}|\\.","")) %>%
      dplyr::mutate(name = paste0("clinvar_path:",symbol,":",stringr::str_replace(hgvs_p,"[A-Z]fs$",""))) %>%
      dplyr::mutate(genomic_change = paste(chrom,pos,ref,alt,sep="_")) %>%
      tidyr::separate_rows(cui) %>%
      dplyr::left_join(medgen_map) %>%
      dplyr::group_by(symbol,hgvs_c,hgvs_p,codon,ref,name,chrom,pos,variation_id,genomic_change,refseq_mrna, ensembl_transcript_id) %>%
      dplyr::summarise(trait = paste(unique(cui_name),collapse=";")) %>%
      dplyr::mutate(link = paste0("<a href='https://www.ncbi.nlm.nih.gov/clinvar/variation/",variation_id,"' target='_blank'>",symbol,"|",hgvs_c,"|",hgvs_p," - ",trait,"</a>")) %>%
      dplyr::group_by(symbol,chrom,codon,refseq_mrna,ensembl_transcript_id,name) %>%
      dplyr::summarise(pathogenic_loci_trait_link = paste(link, collapse=", "), pathogenic_loci_trait = paste0(unique(trait),collapse=";"), variation_id = paste(variation_id, collapse=","), genomic_change = paste(unique(genomic_change),collapse=", "), hgvs_c = paste(unique(hgvs_c), collapse=", "), hgvs_p = paste(unique(hgvs_p), collapse=", "), max_length_ref = max(nchar(ref)), min_pos = min(pos), max_pos = max(pos)) %>%
      dplyr::mutate(start = as.integer(min_pos) - 6, end = max_pos + max_length_ref + 4) %>%
      dplyr::select(-c(max_pos,max_length_ref,min_pos)))
  
  clinvar_pathogenic_cpg$variation_id <- as.character(clinvar_pathogenic_cpg$variation_id)
  
  cpg_pathogenic_insertions_varid <- cpg_pathogenic_insertions %>% dplyr::select(variation_id) %>% tidyr::separate_rows(variation_id)
  cpg_pathogenic_deletions_varid <- cpg_pathogenic_deletions %>% dplyr::select(variation_id) %>% tidyr::separate_rows(variation_id)
  clinvar_pathogenic_cpg <- dplyr::anti_join(clinvar_pathogenic_cpg,dplyr::bind_rows(cpg_pathogenic_insertions_varid,cpg_pathogenic_deletions_varid))
  
  cpg_pathogenic_splicesites <- as.data.frame(
    clinvar_pathogenic_cpg %>%
      dplyr::select(symbol, chrom,pos, ref,alt, refseq_mrna, ref,molecular_consequence, ensembl_transcript_id, hgvs_p, hgvs_c, variation_id, cui) %>%
      dplyr::filter(is.na(hgvs_p)) %>%
      dplyr::filter(stringr::str_detect(molecular_consequence,"splice_")) %>%
      dplyr::mutate(coding_nucleotide_number = stringr::str_match(hgvs_c,"c\\.[0-9]{1,}(\\+|-)[0-9]{1,}(del|ins|dup|A|G|C|T)")[,1]) %>%
      dplyr::mutate(genomic_change = paste(chrom,pos,ref,alt,sep="_")))
  
    cpg_pathogenic_splicesites[is.na(cpg_pathogenic_splicesites$coding_nucleotide_number),]$coding_nucleotide_number <- cpg_pathogenic_splicesites[is.na(cpg_pathogenic_splicesites$coding_nucleotide_number),]$hgvs_c
    
    cpg_pathogenic_splicesites <- as.data.frame(
      cpg_pathogenic_splicesites %>%
        tidyr::separate_rows(cui) %>%
        dplyr::left_join(medgen_map) %>%
        dplyr::group_by(chrom,pos,ref,symbol,hgvs_c,coding_nucleotide_number,variation_id,refseq_mrna, ensembl_transcript_id) %>%
        dplyr::summarise(trait = paste(unique(cui_name),collapse=";"), genomic_change = paste(unique(genomic_change),collapse=";")) %>%
        dplyr::mutate(link = paste0("<a href='https://www.ncbi.nlm.nih.gov/clinvar/variation/",variation_id,"' target='_blank'>",symbol,"|",hgvs_c," - ",trait,"</a>")) %>%
        dplyr::group_by(symbol,chrom,coding_nucleotide_number,refseq_mrna,ensembl_transcript_id) %>%
        dplyr::summarise(pathogenic_loci_trait_link = paste0(link, collapse=", "), pathogenic_loci_trait = paste0(unique(trait),collapse=";"), variation_id = paste(variation_id, collapse=","), genomic_change = paste(unique(genomic_change),collapse=", "), hgvs_c = paste(unique(hgvs_c), collapse=", "), max_length_ref = max(nchar(ref)), min_pos = min(pos), max_pos = max(pos)) %>%
        dplyr::mutate(name = paste0("clinvar_path:",symbol,":",coding_nucleotide_number), start = as.integer(min_pos) - 6, end = as.integer(max_pos) + max_length_ref + 4) %>%
        dplyr::select(-c(max_pos,max_length_ref,min_pos)))
  
  cpg_pathogenic_splicesites_varid <- cpg_pathogenic_splicesites %>% dplyr::select(variation_id) %>% tidyr::separate_rows(variation_id)
  clinvar_pathogenic_cpg <- dplyr::anti_join(clinvar_pathogenic_cpg,cpg_pathogenic_splicesites_varid)
  
  cpg_pathogenic_protein <- as.data.frame(dplyr::select(clinvar_pathogenic_cpg, chrom, pos, ref,alt,symbol, refseq_mrna, molecular_consequence, ensembl_transcript_id, hgvs_c, hgvs_p, variation_id, cui) %>%
    dplyr::filter(!is.na(hgvs_p)) %>%
    tidyr::separate_rows(cui) %>%
    dplyr::filter(stringr::str_detect(molecular_consequence, "nonsense|missense|synonymous|frameshift")) %>%
    dplyr::left_join(medgen_map) %>%
    dplyr::mutate(codon = stringr::str_replace_all(hgvs_p,"\\*$|fs\\*[0-9]{1,}$","")) %>%
    dplyr::mutate(codon = stringr::str_replace_all(codon,"[A-Z]{1,}|[a-z]{1,}|\\.","")) %>%
    dplyr::mutate(genomic_change = paste(chrom,pos,ref,alt,sep="_")) %>%
    dplyr::mutate(hgvs_p = stringr::str_replace_all(hgvs_p,"\\*$","X")) %>%
    dplyr::mutate(name = paste0("clinvar_path:",symbol,":",stringr::str_replace(hgvs_p,"[A-Z]$|delins[A-Z]{1,}$|del[A-Z]{1,}$|[A-Z]{1,}fs(\\*[0-9]{1,}){0,}$",""))) %>%
    dplyr::distinct() %>%
    dplyr::group_by(chrom,pos,genomic_change,name,symbol,hgvs_p,hgvs_c,variation_id,refseq_mrna, ensembl_transcript_id,codon) %>%
    dplyr::summarise(trait = paste(unique(cui_name),collapse=";")) %>%
    dplyr::mutate(link = paste0("<a href='https://www.ncbi.nlm.nih.gov/clinvar/variation/",variation_id,"' target='_blank'>",symbol,"|",hgvs_c,"|",hgvs_p," - ",trait,"</a>")) %>%
    dplyr::group_by(chrom,symbol,codon,refseq_mrna,ensembl_transcript_id,name) %>%
    dplyr::summarise(pathogenic_loci_trait_link = paste0(link, collapse=", "), pathogenic_loci_trait = paste0(unique(trait),collapse=";"), min_pos = min(pos), max_pos = max(pos), variation_id = paste(unique(variation_id), collapse=","),genomic_change = paste(unique(genomic_change),collapse=", "), hgvs_c = paste(unique(hgvs_c), collapse=", "), hgvs_p = paste(unique(hgvs_p), collapse=", "))) 
  
  cpg_pathogenic_protein_single_exon <- cpg_pathogenic_protein %>%
    dplyr::filter(max_pos - min_pos < 50) %>%
    dplyr::mutate(start = min_pos - 6, end = max_pos + 5) %>%
    dplyr::select(-c(min_pos,max_pos))
  
  cpg_pathogenic_protein_exon_spanning <- cpg_pathogenic_protein %>%
    dplyr::filter(max_pos - min_pos >= 50)
  
  all_spanning_entries <- data.frame()
  i <- 1
  while(i <= nrow(cpg_pathogenic_protein_exon_spanning)){
    entry <- as.data.frame(cpg_pathogenic_protein_exon_spanning[i,])
    entry_initial_exon <- entry %>%
      dplyr::mutate(start = min_pos - 6, end = min_pos + 5) %>%
      dplyr::select(-c(min_pos,max_pos))
    entry_subseqent_exon <- entry %>%
      dplyr::mutate(start = max_pos - 6, end = max_pos + 5) %>%
      dplyr::select(-c(min_pos,max_pos))
    all_spanning_entries <- dplyr::bind_rows(entry_initial_exon, entry_subseqent_exon)
    i <- i + 1
    
  }
  
  cpg_pathogenic_protein_regions <- dplyr::bind_rows(cpg_pathogenic_protein_single_exon, all_spanning_entries)
  cpg_pathogenic_deletions$class <- 'Frameshift deletion'
  cpg_pathogenic_insertions$class <- 'Frameshift insertion'
  cpg_pathogenic_splicesites$class <- 'Splice site'
  cpg_pathogenic_protein_regions$class <- 'Codon'
  
  cpg_pathogenic_tsv <- dplyr::bind_rows(cpg_pathogenic_protein_regions,cpg_pathogenic_splicesites,cpg_pathogenic_deletions,cpg_pathogenic_insertions) %>%
    dplyr::select(chrom,start,end,symbol,name,class,variation_id,hgvs_c,hgvs_p,refseq_mrna,ensembl_transcript_id,dplyr::everything())
  
  cpg_pathogenic_tsv$cacao_id <- seq(1:nrow(cpg_pathogenic_tsv))
  
  cpg_pathogenic_traits <- dplyr::select(cpg_pathogenic_tsv, cacao_id, pathogenic_loci_trait) %>%
    tidyr::separate_rows(pathogenic_loci_trait,sep=";") %>%
    dplyr::distinct() %>%
    dplyr::group_by(cacao_id) %>%
    dplyr::summarise(pathogenic_loci_trait = paste(unique(pathogenic_loci_trait),collapse=";"))
  
  cpg_pathogenic_tsv <- as.data.frame(
    cpg_pathogenic_tsv %>%
      dplyr::select(-pathogenic_loci_trait) %>%
      dplyr::left_join(cpg_pathogenic_traits)
  )
  
  cpg_pathogenic_bed <- dplyr::bind_rows(dplyr::select(cpg_pathogenic_protein_regions,chrom,start, end,name), 
                                         dplyr::select(cpg_pathogenic_insertions,chrom,start,end,name), 
                                         dplyr::select(cpg_pathogenic_deletions,chrom,start,end,name), 
                                         dplyr::select(cpg_pathogenic_splicesites,chrom,start,end,name)) %>%
    dplyr::filter(!is.na(chrom) & !is.na(start) & !is.na(end))
  
  
  all_tracks[['bed']][[build]][['clinvar_path']] <- cpg_pathogenic_bed
  all_tracks[['tsv']][[build]][['clinvar_path']] <- cpg_pathogenic_tsv
  
  # cpg_pathogenic_bed_sorted <- sort_bed_regions(cpg_pathogenic_bed)
  # 
  # 
  # write.table(cpg_pathogenic_bed_sorted,file=paste0("data/cancer_hereditary_pathogenic_loci.",datestamp,".",build,".bed"),sep="\t",row.names=F,quote=F,col.names=F)
  # cpg_pathogenic_bed_sorted$chrom <- paste0('chr',cpg_pathogenic_bed_sorted$chrom)
  # write.table(cpg_pathogenic_bed_sorted,file=paste0("data/cancer_hereditary_pathogenic_loci.",datestamp,".",build,".chr.bed"),sep="\t",row.names=F,quote=F,col.names=F)
  # write.table(cpg_pathogenic_tsv,file=paste0("data/cancer_hereditary_pathogenic_loci.",datestamp,".",build,".tsv"),sep="\t",row.names=F,quote=F,col.names=T)
  # system(paste0("cp data/cancer_hereditary_pathogenic_loci.",datestamp,".",build,".tsv data/cancer_hereditary_pathogenic_loci.",build,".tsv"))
  # system(paste0("cp data/cancer_hereditary_pathogenic_loci.",datestamp,".",build,".bed data/cancer_hereditary_pathogenic_loci.",build,".bed"))
  # system(paste0("cp data/cancer_hereditary_pathogenic_loci.",datestamp,".",build,".chr.bed data/cancer_hereditary_pathogenic_loci.",build,".chr.bed"))
  
}



## Sporadic cancer (somatic mutations associated with prognosis, diagnosis or drug sensitivity/resistance)

for(build in c('grch37','grch38')){
  
  civic_biomarkers <- read.table(file=paste0("data-raw/civic.",build,".biomarkers.tsv"),sep="\t",header=T,stringsAsFactors = F,quote="",na.strings=c("","NA"),comment.char="") %>%
    dplyr::select(genesymbol, evidence_id, evidence_type, clinical_significance, variant_name, variant_origin, pubmed_html_link, cancer_type, disease_ontology_id, eitem_codon, eitem_exon, alteration_type, mapping_category, evidence_level, therapeutic_context) %>%
    dplyr::filter(variant_origin == "Somatic Mutation") %>%
    dplyr::filter(mapping_category == "exact" | mapping_category == "codon" | mapping_category == "exon") %>%
    dplyr::filter(evidence_type != "Predisposing") %>%
    dplyr::filter(!is.na(clinical_significance)) %>%
    dplyr::rename(codon = eitem_codon, exon = eitem_exon, symbol = genesymbol) %>%
    dplyr::select(evidence_id, symbol, variant_name, clinical_significance,pubmed_html_link,cancer_type,evidence_type, evidence_level,disease_ontology_id,codon,exon, therapeutic_context) 
  
  civic_regions <- read.table(file=paste0("data-raw/civic.",build,".bed"),skip=1,sep="\t",header=F,stringsAsFactors = F,quote="",na.strings=c("","NA"),comment.char="#")
  colnames(civic_regions) <- c('chrom','start','end','evidence_id')
  civic_regions <- dplyr::inner_join(civic_regions, civic_biomarkers)
  civic_regions$name <- NA
  civic_regions <- dplyr::mutate(civic_regions, variant_name_2 = stringr::str_replace(variant_name,"FS$|FS\\*[0-9]{1,}",""))
  civic_regions[!is.na(civic_regions$exon),]$name <- paste0("civic:",civic_regions[!is.na(civic_regions$exon),]$evidence_id,":",civic_regions[!is.na(civic_regions$exon),]$symbol,":exon",civic_regions[!is.na(civic_regions$exon),]$exon)
  civic_regions[!is.na(civic_regions$codon) & stringr::str_detect(civic_regions$variant_name_2,"[A-Z]{1}[0-9]{1,}"),]$name <- paste0("civic:",civic_regions[!is.na(civic_regions$codon) & stringr::str_detect(civic_regions$variant_name_2,"[A-Z]{1}[0-9]{1,}"),]$evidence_id, ":", civic_regions[!is.na(civic_regions$codon) & stringr::str_detect(civic_regions$variant_name_2,"[A-Z]{1}[0-9]{1,}"),]$symbol,":p.",civic_regions[!is.na(civic_regions$codon) & stringr::str_detect(civic_regions$variant_name_2,"[A-Z]{1}[0-9]{1,}"),]$variant_name_2)
    
  civic_variants <- read.table(file=gzfile(paste0("data-raw/civic.pcgr_acmg.",build,".pass.tsv.gz")),skip=1,sep="\t",header=T,stringsAsFactors = F,quote="",na.strings=c("","NA"),comment.char="") %>%
    dplyr::select(Consequence,SYMBOL,HGVSp_short,CDS_CHANGE,CHROM,POS,REF,ALT,CIVIC_ID) %>%
    tidyr::separate_rows(CIVIC_ID) %>%
    dplyr::rename(consequence = Consequence, evidence_id = CIVIC_ID, hgvsp = HGVSp_short,symbol = SYMBOL, chrom = CHROM, pos = POS, ref = REF, alt = ALT) %>%
    dplyr::mutate(genomic_change = paste(chrom,pos,ref,alt,sep="_")) %>%
    dplyr::mutate(hgvsc = stringr::str_match(CDS_CHANGE,"c\\.[0-9]{1}.+:|c\\.[0-9]{1,}(\\+|-)[0-9]{1,}.+")[,1]) %>%
    dplyr::mutate(ensembl_transcript_id = stringr::str_match(CDS_CHANGE,"ENST[0-9]{1,}")[,1]) %>%
    dplyr::mutate(hgvsc = stringr::str_replace(hgvsc,":exon[0-9]{1,}:","")) %>%
    dplyr::select(chrom, pos, ref, alt, symbol, genomic_change, hgvsp, hgvsc, ensembl_transcript_id, consequence, evidence_id) %>%
    dplyr::distinct()
  
 
  civic_variants <- dplyr::inner_join(civic_variants, civic_biomarkers)
  
  civic_variants$hgvsp <- stringr::str_replace(civic_variants$hgvsp,"p\\.","")
  civic_variants <- dplyr::mutate(civic_variants, variant_name_2 = stringr::str_replace(variant_name,"\\*$","X"))
  civic_variants <- dplyr::mutate(civic_variants, variant_name_2 = stringr::str_replace(variant_name_2,"DUP","dup"))
  civic_variants <- dplyr::mutate(civic_variants, variant_name_2 = stringr::str_replace(variant_name_2,"FS","fs"))
  matching_isoform <- civic_variants[stringr::str_detect(string = civic_variants$variant_name_2, pattern = civic_variants$hgvsp),] %>% dplyr::select(evidence_id,variant_name)
  civic_variants_matching <- dplyr::semi_join(civic_variants,matching_isoform)
  non_matching_isoform <- civic_variants[!stringr::str_detect(string = civic_variants$variant_name_2, pattern = civic_variants$hgvsp),]
  non_matching_isoform <- dplyr::anti_join(non_matching_isoform,matching_isoform)
  civic_variants <- dplyr::bind_rows(civic_variants_matching, non_matching_isoform)
  civic_variants$name <- paste0("civic:",civic_variants$evidence_id, ":",civic_variants$symbol,":p.",stringr::str_replace(civic_variants$hgvsp,"[A-Z]$|dup$|delins[A-Z]{1,}$|fs$|ins[A-Z]{1,}$",""))
  civic_variants[civic_variants$hgvsp == '.' & !is.na(civic_variants$hgvsc),]$name <- paste0("civic:",civic_variants[civic_variants$hgvsp == '.' & !is.na(civic_variants$hgvsc),]$evidence_id,":",civic_variants[civic_variants$hgvsp == '.' & !is.na(civic_variants$hgvsc),]$symbol,":",civic_variants[civic_variants$hgvsp == '.' & !is.na(civic_variants$hgvsc),]$hgvsc)
  
  civic_variants <- as.data.frame(
    civic_variants %>%
    dplyr::group_by(symbol, chrom, hgvsp, consequence, evidence_id, variant_name, clinical_significance, pubmed_html_link, cancer_type, evidence_type, evidence_level,
                                     disease_ontology_id, codon, exon, therapeutic_context, variant_name_2, name) %>%
    dplyr::summarise(hgvsc = paste(unique(hgvsc),collapse=", "), genomic_change = paste(unique(genomic_change),collapse=", "), 
                     ensembl_transcript_id = paste(unique(ensembl_transcript_id), collapse=", "), 
                     start = min(pos) - 6, 
                     end = min(pos) + max(nchar(ref),nchar(alt)) + 4) %>%
      dplyr::arrange(chrom,start,end)
  )
  
  if(nrow(civic_variants[civic_variants$symbol == 'TERT' & civic_variants$consequence == 'upstream_gene_variant',]) > 0){
    civic_variants[civic_variants$symbol == 'TERT' & civic_variants$consequence == 'upstream_gene_variant',]$name <- "TERT:upstream_gene_variant_C228T"
  }
  if(nrow(civic_variants[civic_variants$symbol == 'CDKN2A' & civic_variants$consequence == 'upstream_gene_variant',]) > 0){
    civic_variants[civic_variants$symbol == 'CDKN2A' & civic_variants$consequence == 'upstream_gene_variant',]$name <- "CDKN2A:upstream_gene_variant_CA12951936"
  }
  
  civic_all <- dplyr::bind_rows(civic_variants,civic_regions) %>% dplyr::arrange(chrom,start,end) %>% dplyr::distinct()
  civic_all_bed <- dplyr::select(civic_all,chrom,start,end,name) %>% dplyr::arrange(chrom,start,end)
  
  cancer_somatic_actionable_loci_tsv <- civic_all
  
  
  all_tracks[['bed']][[build]][['civic']] <- civic_all_bed
  all_tracks[['tsv']][[build]][['civic']] <- cancer_somatic_actionable_loci_tsv
  
  # cancer_somatic_actionable_loci_bed_sorted <- sort_bed_regions(civic_all_bed)
  # 
  # write.table(cancer_somatic_actionable_loci_bed_sorted,file=paste0("data/cancer_somatic_actionable_loci.",datestamp,".",build,".bed"),sep="\t",row.names=F,quote=F,col.names=F)
  # cancer_somatic_actionable_loci_bed_sorted$chrom <- paste0('chr',cancer_somatic_actionable_loci_bed_sorted$chrom)
  # write.table(cancer_somatic_actionable_loci_bed_sorted,file=paste0("data/cancer_somatic_actionable_loci.",datestamp,".",build,".chr.bed"),sep="\t",row.names=F,quote=F,col.names=F)
  # write.table(cancer_somatic_actionable_loci_tsv,file=paste0("data/cancer_somatic_actionable_loci.",datestamp,".",build,".tsv"),sep="\t",row.names=F,quote=F,col.names=T)
  # system(paste0("cp data/cancer_somatic_actionable_loci.",datestamp,".",build,".tsv data/cancer_somatic_actionable_loci.",build,".tsv"))
  # system(paste0("cp data/cancer_somatic_actionable_loci.",datestamp,".",build,".bed data/cancer_somatic_actionable_loci.",build,".bed"))
  # system(paste0("cp data/cancer_somatic_actionable_loci.",datestamp,".",build,".chr.bed data/cancer_somatic_actionable_loci.",build,".chr.bed"))
  
}


## Somatic mutation hotspots (cancer_hotspots.org)

for(build in c('grch37','grch38')){
  
    hotspot_variants <- read.table(file=paste0("data-raw/cancer_hotspots.",build,".tsv"),skip=1,sep="\t",header=T,stringsAsFactors = F,quote="",na.strings=c("","NA"),comment.char="") %>%
      dplyr::select(MUTATION_HOTSPOT,MUTATION_HOTSPOT_CANCERTYPE,MUTATION_HOTSPOT_TRANSCRIPT,CHROM,POS,REF,ALT) %>%
      tidyr::separate_rows(MUTATION_HOTSPOT_CANCERTYPE,sep=",") %>%
      dplyr::rename(hotspot = MUTATION_HOTSPOT, hotspot_cancer_type = MUTATION_HOTSPOT_CANCERTYPE, chrom = CHROM, pos = POS, ref = REF, alt = ALT) %>%
      dplyr::mutate(genomic_change = paste(chrom,pos,ref,alt,sep="_")) %>%
      dplyr::select(chrom, pos, ref, alt, genomic_change, hotspot, hotspot_cancer_type) %>%
      dplyr::distinct()
    
    hotspot_variants_any <- as.data.frame(
      hotspot_variants %>%
      dplyr::mutate(hotspot_cancer_type = stringr::str_replace(hotspot_cancer_type,"esophagusstomach","Esophagus/Stomach")) %>%
        dplyr::mutate(hotspot_cancer_type = stringr::str_replace(hotspot_cancer_type,"lymph","Lymph")) %>%
        dplyr::mutate(hotspot_cancer_type = stringr::str_replace(hotspot_cancer_type,"blood","Blood")) %>%
        dplyr::mutate(hotspot_cancer_type = stringr::str_replace(hotspot_cancer_type,"peritoneum","Peritoneum")) %>%
        dplyr::mutate(hotspot_cancer_type = stringr::str_replace(hotspot_cancer_type,"skin","Skin")) %>%
        dplyr::mutate(hotspot_cancer_type = stringr::str_replace(hotspot_cancer_type,"breast","Breast")) %>%
        dplyr::mutate(hotspot_cancer_type = stringr::str_replace(hotspot_cancer_type,"pancreas","Pancreas")) %>%
        dplyr::mutate(hotspot_cancer_type = stringr::str_replace(hotspot_cancer_type,"bowel","Bowel")) %>%
        dplyr::mutate(hotspot_cancer_type = stringr::str_replace(hotspot_cancer_type,"uterus","Uterus")) %>%
        dplyr::mutate(hotspot_cancer_type = stringr::str_replace(hotspot_cancer_type,"prostate","Prostate")) %>%
        dplyr::mutate(hotspot_cancer_type = stringr::str_replace(hotspot_cancer_type,"cervix","Cervix")) %>%
        dplyr::mutate(hotspot_cancer_type = stringr::str_replace(hotspot_cancer_type,"lung","Lung")) %>%
        dplyr::mutate(hotspot_cancer_type = stringr::str_replace(hotspot_cancer_type,"cnsbrain","CNS/brain")) %>%
        dplyr::mutate(hotspot_cancer_type = stringr::str_replace(hotspot_cancer_type,"bladder","Bladder/Urinary Tract")) %>%
        dplyr::mutate(hotspot_cancer_type = stringr::str_replace(hotspot_cancer_type,"kidney","Kidney")) %>%
        dplyr::mutate(hotspot_cancer_type = stringr::str_replace(hotspot_cancer_type,"testis","Testis")) %>%
        dplyr::mutate(hotspot_cancer_type = stringr::str_replace(hotspot_cancer_type,"liver","Liver")) %>%
        dplyr::mutate(hotspot_cancer_type = stringr::str_replace(hotspot_cancer_type,"headandneck","Head and Neck")) %>%
        dplyr::mutate(hotspot_cancer_type = stringr::str_replace(hotspot_cancer_type,"ovaryfallopiantube","Ovary/Fallopian Tube")) %>%
        dplyr::mutate(hotspot_cancer_type = stringr::str_replace(hotspot_cancer_type,"thyroid","Thyroid")) %>%
        dplyr::mutate(hotspot_cancer_type = stringr::str_replace(hotspot_cancer_type,"biliarytract","Biliary Tract")) %>%
        dplyr::mutate(hotspot_cancer_type = stringr::str_replace(hotspot_cancer_type,"other","Other")) %>%
        dplyr::mutate(hotspot_cancer_type = stringr::str_replace(hotspot_cancer_type,"vulvavagina","Vulva/Vagina")) %>%
        dplyr::mutate(hotspot_cancer_type = stringr::str_replace(hotspot_cancer_type,"penis","Penis")) %>%
        dplyr::mutate(hotspot_cancer_type = stringr::str_replace(hotspot_cancer_type,"softtissue","Soft Tissue")) %>%
        dplyr::mutate(hotspot_cancer_type = stringr::str_replace(hotspot_cancer_type,"thymus","Thymus")) %>%
        dplyr::mutate(hotspot_cancer_type = stringr::str_replace(hotspot_cancer_type,"eye","Eye")) %>%
        dplyr::mutate(hotspot_cancer_type = stringr::str_replace(hotspot_cancer_type,"pleura","Pleura")) %>%
        dplyr::mutate(hotspot_cancer_type = stringr::str_replace(hotspot_cancer_type,"ampullaofvater","Ampulla of Vater")) %>%
        dplyr::mutate(hotspot_cancer_type = stringr::str_replace(hotspot_cancer_type,":[0-9]{1,}","")) %>%
        dplyr::mutate(hotspot_cancer_type = stringr::str_replace(hotspot_cancer_type,"bone","Bone")) %>%
        dplyr::mutate(hotspot_cancer_type = stringr::str_replace(hotspot_cancer_type,"adrenalgland","Adrenal Gland")) %>%
        
      #dplyr::select(hotspot_variants, -hotspot_cancer_types) %>%
        dplyr::group_by(chrom, hotspot) %>%
        dplyr::summarise(genomic_change = paste(unique(genomic_change),collapse=", "), 
                         cancer_type = paste(unique(hotspot_cancer_type), collapse=", "), 
                         min_pos = min(pos), 
                         max_pos = max(pos) + max(nchar(ref),nchar(alt))) %>%
        dplyr::mutate(name = paste0("hotspot:",hotspot)) %>%
        tidyr::separate(hotspot,c('symbol','hgvsp','pvalue'),sep="\\|",remove=F) %>%
        dplyr::mutate(hgvsp = paste0('p.',hgvsp)) %>%
        dplyr::mutate(start = min_pos - 6, end = max_pos + 4) %>%
        dplyr::select(-c(max_pos,min_pos)))
    
    cancer_somatic_hotspot_loci_tsv <- hotspot_variants_any %>% dplyr::arrange(chrom,start,end)
    cancer_somatic_hotspot_loci_bed_sorted <- sort_bed_regions(dplyr::select(hotspot_variants_any,chrom,start,end,name))
    
    all_tracks[['bed']][[build]][['hotspot']] <- cancer_somatic_hotspot_loci_bed_sorted
    all_tracks[['tsv']][[build]][['hotspot']] <- cancer_somatic_hotspot_loci_tsv
    
    
    # 
    # write.table(cancer_somatic_hotspot_loci_bed_sorted,file=paste0("data/cancer_somatic_hotspot_loci.",datestamp,".",build,".bed"),sep="\t",row.names=F,quote=F,col.names=F)
    # cancer_somatic_hotspot_loci_bed_sorted$chrom <- paste0('chr',cancer_somatic_hotspot_loci_bed_sorted$chrom)
    # write.table(cancer_somatic_hotspot_loci_bed_sorted,file=paste0("data/cancer_somatic_hotspot_loci.",datestamp,".",build,".chr.bed"),sep="\t",row.names=F,quote=F,col.names=F)
    # write.table(cancer_somatic_hotspot_loci_tsv,file=paste0("data/cancer_somatic_hotspot_loci.",datestamp,".",build,".tsv"),sep="\t",row.names=F,quote=F,col.names=T)
    # system(paste0("cp data/cancer_somatic_hotspot_loci.",datestamp,".",build,".tsv data/cancer_somatic_hotspot_loci.",build,".tsv"))
    # system(paste0("cp data/cancer_somatic_hotspot_loci.",datestamp,".",build,".bed data/cancer_somatic_hotspot_loci.",build,".bed"))
    # system(paste0("cp data/cancer_somatic_hotspot_loci.",datestamp,".",build,".chr.bed data/cancer_somatic_hotspot_loci.",build,".chr.bed"))
    
    # for(ctype in unique(stringr::str_replace(hotspot_variants$hotspot_cancer_types,":[0-9]{1,}",""))){
    #   vars <- hotspot_variants[stringr::str_detect(hotspot_variants$hotspot_cancer_types,ctype),] %>%
    #     dplyr::select(-hotspot_cancer_types) %>%
    #     dplyr::group_by(chrom, hotspot) %>%
    #     dplyr::summarise(genomic_change = paste(unique(genomic_change),collapse=","), min_pos = min(pos), max_pos = max(pos) + max(nchar(ref),nchar(alt))) %>%
    #     dplyr::mutate(name = hotspot) %>%
    #     dplyr::mutate(start = min_pos - 6, end = max_pos + 5) %>%
    #     dplyr::select(-c(max_pos,min_pos))
    #   
    #   cancer_somatic_hotspot_loci_bed_sorted <- sort_bed_regions(dplyr::select(vars,chrom,start,end,name))
    #   write.table(cancer_somatic_hotspot_loci_bed_sorted,file=paste0("data/cancer_somatic_hotspot_loci.",ctype,".",datestamp,".",build,".bed"),sep="\t",row.names=F,quote=F,col.names=F)
    #   cancer_somatic_hotspot_loci_bed_sorted$chrom <- paste0('chr',cancer_somatic_hotspot_loci_bed_sorted$chrom)
    #   write.table(cancer_somatic_hotspot_loci_bed_sorted,file=paste0("data/cancer_somatic_hotspot_loci.",ctype,".",datestamp,".",build,".chr.bed"),sep="\t",row.names=F,quote=F,col.names=F)
    #   write.table(cancer_somatic_hotspot_loci_tsv,file=paste0("data/cancer_somatic_hotspot_loci.",ctype,".",datestamp,".",build,".tsv"),sep="\t",row.names=F,quote=F,col.names=T)
    #   system(paste0("cp data/cancer_somatic_hotspot_loci.",ctype,".",datestamp,".",build,".tsv data/cancer_somatic_hotspot_loci.",ctype,".",build,".tsv"))
    #   system(paste0("cp data/cancer_somatic_hotspot_loci.",ctype,".",datestamp,".",build,".bed data/cancer_somatic_hotspot_loci.",ctype,".",build,".bed"))
    #   system(paste0("cp data/cancer_somatic_hotspot_loci.",ctype,".",datestamp,".",build,".chr.bed data/cancer_somatic_hotspot_loci.",ctype,".",build,".chr.bed"))
    #   cat(ctype,sep='\n')
    # }
  
}

for(build in c('grch37','grch38')){
  global_bed <- sort_bed_regions(rbind(all_tracks[['bed']][[build]][['clinvar_path']], all_tracks[['bed']][[build]][['civic']], all_tracks[['bed']][[build]][['hotspot']]))
  write.table(global_bed,file=paste0("data/cacao.",datestamp,".",build,".bed"),sep="\t",row.names=F,quote=F,col.names=F)
  global_bed$chrom <- paste0('chr',global_bed$chrom)
  write.table(global_bed,file=paste0("data/cacao.",datestamp,".",build,".chr.bed"),sep="\t",row.names=F,quote=F,col.names=F)
  system(paste0("cp data/cacao.",datestamp,".",build,".chr.bed data/cacao.",build,".chr.bed"))
  system(paste0("cp data/cacao.",datestamp,".",build,".bed data/cacao.",build,".bed"))
  
  for(ctype in c('clinvar_path','civic','hotspot')){
    write.table(all_tracks[['tsv']][[build]][[ctype]],file=paste0("data/cacao.",ctype,".",datestamp,".",build,".tsv"),sep="\t",row.names=F,quote=F,col.names=T)
    system(paste0("cp data/cacao.",ctype,".",datestamp,".",build,".tsv data/cacao.",ctype,".",build,".tsv"))
  }
  
}

 