#!/usr/bin/env python

import csv
import re
import argparse
import os
import subprocess
import logging
import sys

cacao_version = '0.2.0'


def __main__():
   
   parser = argparse.ArgumentParser(description='cacao - assessment of sequencing coverage at pathogenic and actionable loci in cancer',formatter_class=argparse.ArgumentDefaultsHelpFormatter, usage="%(prog)s [options] <BAM-or-CRAM> <BED_TARGET> <ALN_FNAME_HOST> <TARGET_FNAME_HOST> <BED_TRACK_DIRECTORY> OUTPUT_DIR> <GENOME_ASSEMBLY> <CANCER_MODE> <SAMPLE_ID> <MAPQ> <THREADS> <CALLABILITY_LEVELS_GERMLINE> <CALLABILITY_LEVELS_SOMATIC>")
   parser.add_argument('aln',help='Alignment file (BAM/CRAM)')
   parser.add_argument('target',help='BED file')
   parser.add_argument('host_name_aln',help='Host alignment filename (BAM/CRAM)')
   parser.add_argument('host_name_target',help='Host target filename (BED)')
   parser.add_argument('bed_track_directory', help='Directory with BED tracks of pathogenic/actionable cancer loci for grch37/grch38')
   parser.add_argument('output_directory', help='Output directory')
   parser.add_argument('genome_assembly',choices = ['grch37','grch38'], help='Human genome assembly build: grch37 or grch38')
   parser.add_argument('mode',choices=['hereditary','somatic','any'],help="Choice of clinical cancer context (hereditary/somatic/any)")
   parser.add_argument('sample_id', help="Sample identifier - prefix for output files")
   parser.add_argument('mapq', default = 0, type=int, help="mapping quality threshold")
   parser.add_argument('threads',default = 0, type=int, help='Number of mosdepth BAM decompression threads. (use 4 or fewer)')
   parser.add_argument('callability_levels_germline', default="0:10:100", help="Intervals of sequencing depth denoting e.g. NO_COVERAGE (0), LOW_COVERAGE (1-9), CALLABLE (10-99), HIGH_COVERAGE (>= 100)")
   parser.add_argument('callability_levels_somatic', default="0:30:200", help="Intervals of sequencing depth denoting e.g. NO_COVERAGE (0), LOW_COVERAGE (1-29), CALLABLE (30-199), HIGH_COVERAGE (>= 200)")
   parser.add_argument('--prbase', action = "store_true", help='By default, the cacao workflow will only assess coverage at clinically relevant variant loci in cancer. This option enables a full coverage calculation of the alignment (will increase runtime substantially)')


   args = parser.parse_args()

   logger = getlogger('cacao-check-input')
   track_info = check_input_files(args.aln, args.target, args.sample_id, args.output_directory, args.genome_assembly, args.mode, logger)
   logger = getlogger('cacao-coverage-assessment')
   coverage_bed_track = "NA"
   coverage_tsv_tracks = {}
   sample_postfix = "_" + str(args.genome_assembly) + "_coverage_cacao"
   sample_postfix_prbase = "_prbase_" + str(args.genome_assembly) + "_coverage_cacao"


   for m in ['hereditary','somatic_actionable','somatic_hotspot']:       
      coverage_tsv_tracks[m] = os.path.join(str(args.bed_track_directory),track_info['tsv'][m])
   coverage_description = ""
   if args.mode == "hereditary":
      coverage_description = "pathogenic loci in hereditary cancer"
   if args.mode == "somatic":
      coverage_description = "actionable somatic loci in cancer and cancer mutation hotspots"
   if args.mode == "any":
      coverage_description = "actionable somatic loci in cancer, cancer mutation hotspots, and pathogenic loci in hereditary cancer"
   
   if track_info['cacao_loci_bed'] != "NA":
      logger.info("Determination of sequencing coverage with https://github.com/brentp/mosdepth:")
      logger.info("Loci subject to coverage assessment: " + str(coverage_description))
      mosdepth_cmd = 'mosdepth --no-per-base --by ' + os.path.join(str(args.bed_track_directory),str(track_info['cacao_loci_bed'])) + ' --mapq ' + str(args.mapq) + ' --threads ' + str(args.threads) + ' ' + str(args.sample_id) + str(sample_postfix) + ' ' + str(args.aln)
      mosdepth_prbase_cmd = 'mosdepth --mapq ' + str(args.mapq) + ' --threads ' + str(args.threads) + ' ' + str(args.sample_id) + str(sample_postfix_prbase) + ' ' + str(args.aln)
      logger.info('command: ' + str(mosdepth_cmd))
      check_subprocess(mosdepth_cmd)
      if args.prbase is True:
         check_subprocess(mosdepth_prbase_cmd)
   bed_coverage_track_gz = str(args.sample_id) + str(sample_postfix) + '.regions.bed.gz'
   coverage_bed_track = os.path.join(str(args.output_directory),str(args.sample_id) + str(sample_postfix) + '.regions.bed')
   if os.path.exists(bed_coverage_track_gz) and os.path.getsize(bed_coverage_track_gz) > 0:
      logger.info('Decompressing BED file (' + str(bed_coverage_track_gz) + ') with coverage pr. loci')
      check_subprocess('bgzip -dc ' + str(bed_coverage_track_gz) + ' > ' + coverage_bed_track)
   else:
      error_message("Resulting mosdepth coverage file " + str(bed_coverage_track_gz) + "is non-existent or empty", logger)

   coverage_bed_track_target = os.path.join(str(args.output_directory),str(args.sample_id) + str(sample_postfix) + '.regions_target.bed')
   if track_info['query_target_bed'] != "NA":
      logger.info('Limiting coverage assessment to query target regions')
      #intersect_target_bed = 'bedtools intersect -wa -u -a ' + str(coverage_bed_track) + ' -b ' + str(track_info['query_target_bed'] ) + ' > ' + str(coverage_bed_track_target)
      annotate_target_bed = 'bedtools annotate -i ' + str(coverage_bed_track) + ' -files ' + str(track_info['query_target_bed'] ) + ' > ' + str(coverage_bed_track_target)
      check_subprocess(annotate_target_bed)
   else:
      cp_cmd = 'cp ' + str(coverage_bed_track) + ' ' + str(coverage_bed_track_target)
      check_subprocess(cp_cmd)

    
   cacao_report_parameters = []
   cacao_report_parameters.append(coverage_bed_track_target)
   cacao_report_parameters.append(args.host_name_aln)
   cacao_report_parameters.append(args.host_name_target)
   cacao_report_parameters.append(coverage_tsv_tracks['hereditary'])
   cacao_report_parameters.append(coverage_tsv_tracks['somatic_actionable'])
   cacao_report_parameters.append(coverage_tsv_tracks['somatic_hotspot'])
   cacao_report_parameters.append(str(args.sample_id))
   cacao_report_parameters.append(str(args.mode))
   cacao_report_parameters.append(str(args.callability_levels_germline))
   cacao_report_parameters.append(str(args.callability_levels_somatic))
   cacao_report_parameters.append(str(args.mapq))
   cacao_report_parameters.append(str(args.genome_assembly))
   cacao_report_parameters.append(str(cacao_version))
   cacao_report_parameters.append(str(args.output_directory))

   report_R_command = "/cacao.R " + " ".join(cacao_report_parameters)
   #print(str(report_R_command))
   check_subprocess(report_R_command)

   logger.info('Finished')


def check_input_files(aln_fname, target_bed_fname, sample_id, output_directory, genome_assembly, mode, logger):

   logger.info('Determination of BED region file - considering genome assembly, cancer mode, and chromosome naming convention')
   idxstats_fname = os.path.join(output_directory, str(sample_id) + '_' + str(genome_assembly) + '.idxstats.tsv')
   chrnaming_cmd = 'samtools idxstats ' + str(aln_fname) + ' > ' + idxstats_fname
   check_subprocess(chrnaming_cmd)

   chr_naming = False
   if os.path.exists(idxstats_fname) and os.path.getsize(idxstats_fname) > 0:
      f = open(idxstats_fname,'r')
      for line in f:
         chrom_stats = line.rstrip().split('\t')
         if len(chrom_stats) > 0:
            if chrom_stats[0].startswith('chr'):
                chr_naming = True
      
      f.close()
   
   chr_naming_bed = 'NA'
   target_bed = False
   if os.path.exists(target_bed_fname) and os.path.getsize(target_bed_fname) > 0:
      target_bed = True
      chr_naming_bed = False
      f = open(target_bed_fname,'r')
      for line in f:
         target_bed_data = line.rstrip().split('\t')
         if len(target_bed_data) >= 2:
            if not (target_bed_data[1].isdigit() and target_bed_data[2].isdigit()):
               #valid_bed = False
               msg = "Non-integer coordinates in target BED file: " + str(target_bed_data[0]) + " " + str(target_bed_data[1]) + " " + str(target_bed_data[2])
               error_message(msg, logger)
            if target_bed_data[0].startswith('chr'):
               chr_naming_bed = True

      f.close()

   if (chr_naming is False and chr_naming_bed is True) or (chr_naming is True and chr_naming_bed is False):
      error_message('Chromosome naming in target BED file and alignment file does not correspond', logger)

   #if valid_bed is True:
      #error_message('Chromosome naming in target BED file and alignment file does not correspond', logger)
      #aln_target_fname = 'balle.target.bam'
      #intersect_target_bam = 'bedtools intersect -wa -u -a ' + str(aln_fname) + ' -b ' + str(target_bed_fname) + ' > ' + str(aln_target_fname)
      #os.system(intersect_target_bam)
      ## intersect BAM file with target BED

   

   track_info = {}
   track_info['cacao_loci_bed'] = "NA"
   track_info['query_target_bed'] = "NA"
   if target_bed is True:
      track_info['query_target_bed'] = target_bed_fname
   track_info['tsv'] = {}
   for m in ['hereditary','somatic_hotspot','somatic_actionable']:
      track_info['tsv'][m] = "NA"

   track_info['cacao_loci_bed'] = 'cacao.' + str(genome_assembly) + '.bed'
   if chr_naming is True:
      track_info['cacao_loci_bed'] = 'cacao.' + str(genome_assembly) + '.chr.bed'

   if mode == 'hereditary' or mode == 'any':
      tsv_fname = 'cacao.clinvar_path.' + str(genome_assembly) + '.tsv'
      track_info['tsv']['hereditary'] = tsv_fname
   if mode == 'somatic' or mode == 'any':
      tsv_fname = 'cacao.civic.' + str(genome_assembly) + '.tsv'
      track_info['tsv']['somatic_actionable'] = tsv_fname
      tsv_fname = 'cacao.hotspot.' + str(genome_assembly) + '.tsv'
      track_info['tsv']['somatic_hotspot'] = tsv_fname

   return track_info


def check_subprocess(command):
   #print(command)
   try:
      output = subprocess.check_output(str(command), stderr=subprocess.STDOUT, shell=True)
      if len(output) > 0:
         print (str(output.decode()).rstrip())
   except subprocess.CalledProcessError as e:
      print (e.output.decode())
      exit(0)

def error_message(message, logger):
   logger.error('')
   logger.error(message)
   logger.error('')
   exit(0)

def warn_message(message, logger):
   logger.warning(message)


def getlogger(logger_name):
   logger = logging.getLogger(logger_name)
   logger.setLevel(logging.DEBUG)

   # create console handler and set level to debug
   ch = logging.StreamHandler(sys.stdout)
   ch.setLevel(logging.DEBUG)

   # add ch to logger
   logger.addHandler(ch)

   # create formatter
   formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s", "20%y-%m-%d %H:%M:%S")

   #add formatter to ch
   ch.setFormatter(formatter)

   return logger

if __name__=="__main__": __main__()

