#!/usr/bin/env python

import csv
import re
import argparse
import os
import subprocess
import logging
import sys

cacao_version = '1.0.0'


def __main__():
   
   parser = argparse.ArgumentParser(description='cacao - assessment of sequencing coverage at pathogenic and actionable loci in cancer',formatter_class=argparse.ArgumentDefaultsHelpFormatter, usage="%(prog)s [options] <BAM-or-CRAM> <BED_TRACK_DIRECTORY> OUTPUT_DIR> <GENOME_ASSEMBLY> <CANCER_MODE> <SAMPLE_ID> <MAPQ> <THREADS> <CALLABILITY_LEVELS>")
   parser.add_argument('aln',help='Alignment file (BAM/CRAM)')
   parser.add_argument('bed_track_directory', help='Directory with BED tracks of pathogenic/actionable cancer loci for grch37/grch38')
   parser.add_argument('output_directory', help='Output directory')
   parser.add_argument('genome_assembly',choices = ['grch37','grch38'], help='Human genome assembly build: grch37 or grch38')
   parser.add_argument('mode',choices=['hereditary','somatic','any'],help="Choice of clinical cancer context (hereditary/somatic/any)")
   parser.add_argument('sample_id', help="Sample identifier - prefix for output files")
   parser.add_argument('mapq', default = 0, type=int, help="mapping quality threshold")
   parser.add_argument('threads',default = 0, type=int, help='Number of mosdepth BAM decompression threads. (use 4 or fewer)')
   parser.add_argument('callability_levels_germline', default="0:10:150", help="Intervals of sequencing depth denoting e.g. NO_COVERAGE (0), LOW_COVERAGE (1-9), CALLABLE (10-149), HIGH_COVERAGE (>= 150)")
   parser.add_argument('callability_levels_somatic', default="0:30:200", help="Intervals of sequencing depth denoting e.g. NO_COVERAGE (0), LOW_COVERAGE (1-29), CALLABLE (30-199), HIGH_COVERAGE (>= 200)")


   args = parser.parse_args()

   logger = getlogger('cacao-get-track')
   track_info = get_track_file(args.aln, args.sample_id, args.output_directory, args.genome_assembly, args.mode, logger)
   logger = getlogger('cacao-coverage-assessment')
   coverage_bed_tracks = {}
   coverage_tsv_tracks = {}
   for m in ['hereditary','somatic_actionable','somatic_hotspot']:
       
       coverage_bed_tracks[m] = "NA"
       coverage_tsv_tracks[m] = os.path.join(str(args.bed_track_directory),track_info['tsv'][m])
       sample_postfix = "_" + str(m)

       coverage_description = "pathogenic loci in hereditary cancer"
       if m == "somatic_actionable":
           coverage_description = "actionable somatic mutations - implications for prognosis/diagnosis/drug response"
       if m == "somatic_hotspot":
           coverage_description = "somatic mutation hotspots - driver events in cancer"
       if track_info['bed'][m] != "NA":
          logger.info("Determination of coverage in target regions with https://github.com/brentp/mosdepth: " + str(coverage_description))
          mosdepth_cmd = 'mosdepth --no-per-base --by ' + os.path.join(str(args.bed_track_directory),str(track_info['bed'][m])) + ' --mapq ' + str(args.mapq) + ' --threads ' + str(args.threads) + ' ' + str(args.sample_id) + str(sample_postfix) + ' ' + str(args.aln)
          logger.info('command: ' + str(mosdepth_cmd))
          check_subprocess(mosdepth_cmd)
       bed_coverage_track_gz = str(args.sample_id) + str(sample_postfix) + '.regions.bed.gz'
       coverage_bed_tracks[m] = os.path.join(str(args.output_directory),str(args.sample_id) + str(sample_postfix) + '.regions.bed')

       if os.path.exists(bed_coverage_track_gz) and os.path.getsize(bed_coverage_track_gz) > 0:
           logger.info('Decompressing BED file (' + str(bed_coverage_track_gz) + ') with coverage pr. loci')
           check_subprocess('bgzip -dc ' + str(bed_coverage_track_gz) + ' > ' + coverage_bed_tracks[m])
       else:
           coverage_bed_tracks[m] = "NA"
    
   #print(str(coverage_bed_tracks)) 
   report_R_command = '/cacao.R ' + str(coverage_bed_tracks['hereditary']) + ' ' + str(coverage_tsv_tracks['hereditary']) + ' ' + str(coverage_bed_tracks['somatic_actionable']) + ' ' + str(coverage_tsv_tracks['somatic_actionable']) + ' ' + str(coverage_bed_tracks['somatic_hotspot']) + ' ' + str(coverage_tsv_tracks['somatic_hotspot']) + ' ' + str(args.sample_id) + ' ' + str(args.mode) + ' ' + str(args.callability_levels_germline) + ' ' + str(args.callability_levels_somatic) + ' ' + str(args.mapq) + ' ' + str(args.genome_assembly) + ' ' + str(cacao_version) + ' ' + str(args.output_directory)
   check_subprocess(report_R_command)

   logger.info('Finished')


def get_track_file(aln_fname, sample_id, output_directory, genome_assembly, mode, logger):

   logger.info('Determination of BED region file - considering genome assembly, cancer mode, and chromosome naming convention')
   idxstats_fname = os.path.join(output_directory, str(sample_id) + '.idxstats.tsv')
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
   
   track_info = {}
   track_info['bed'] = {}
   track_info['tsv'] = {}
   for m in ['hereditary','somatic_hotspot','somatic_actionable']:
      track_info['bed'][m] = "NA"
      track_info['tsv'][m] = "NA"

   if mode == 'hereditary' or mode == 'any':
      bed_track = 'cancer_hereditary_pathogenic_loci.' + str(genome_assembly) + '.bed'
      if chr_naming is True:
          bed_track = 'cancer_hereditary_pathogenic_loci.' + str(genome_assembly) + '.chr.bed'
      tsv_fname = 'cancer_hereditary_pathogenic_loci.' + str(genome_assembly) + '.tsv'
      track_info['bed']['hereditary'] = bed_track
      track_info['tsv']['hereditary'] = tsv_fname
   if mode == 'somatic' or mode == 'any':
      bed_track = 'cancer_somatic_actionable_loci.' + str(genome_assembly) + '.bed'
      if chr_naming is True:
          bed_track = 'cancer_somatic_actionable_loci.' + str(genome_assembly) + '.chr.bed'
      tsv_fname = 'cancer_somatic_actionable_loci.' + str(genome_assembly) + '.tsv'
      track_info['bed']['somatic_actionable'] = bed_track
      track_info['tsv']['somatic_actionable'] = tsv_fname

      bed_track = 'cancer_somatic_hotspot_loci.' + str(genome_assembly) + '.bed'
      if chr_naming is True:
          bed_track = 'cancer_somatic_hotspot_loci.' + str(genome_assembly) + '.chr.bed'
      tsv_fname = 'cancer_somatic_hotspot_loci.' + str(genome_assembly) + '.tsv'
      track_info['bed']['somatic_hotspot'] = bed_track
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

