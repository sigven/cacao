#!/usr/bin/env python

import csv
import re
import argparse
import os
import subprocess
import logging
import sys
import toml

cacao_version = '1.0.0'


def __main__():
   
   parser = argparse.ArgumentParser(description='cacao - assessment of sequencing coverage at pathogenic and actionable loci in cancer',formatter_class=argparse.ArgumentDefaultsHelpFormatter, usage="%(prog)s [options] <BAM-or-CRAM> <TRACK_DIRECTORY> OUTPUT_DIR> <GENOME_ASSEMBLY> <MODE> <SAMPLE_ID>")
   parser.add_argument('--mapq', dest = "mapq", default = 0, type=int, help='mapping quality threshold')
   parser.add_argument('--threads', dest = "threads", default = 0, type=int, help='Number of mosdepth BAM decompression threads. (use 4 or fewer)')
   parser.add_argument('--callability_levels_germline', dest="callability_levels_germline",default="0:10:100", help="Simple colon-separated string that defines four levels of variant callability: NO_COVERAGE (0), LOW_COVERAGE (0-9), CALLABLE (10-99), HIGH_COVERAGE ( > 100). Initial value must be 0.")
   parser.add_argument('--callability_levels_somatic', dest="callability_levels_somatic",default="0:30:200", help="Simple colon-separated string that defines four levels of variant callability: NO_COVERAGE (0), LOW_COVERAGE (0-29), CALLABLE (30-199), HIGH_COVERAGE (> 200). Initial value must be 0.")
   parser.add_argument('--target', dest = "target", help="BED file with target regions subject to sequencing")
   parser.add_argument('--force_overwrite', action = "store_true", help='By default, the script will fail with an error if any output file already exists. You can force the overwrite of existing result files by using this flag')
   parser.add_argument('--no-docker', action='store_true', dest='no_docker', default=False, help='Run the cacao workflow in a non-Docker mode (see install_no_docker/ folder for instructions')
   parser.add_argument('--version', action='version', version='%(prog)s ' + str(cacao_version))
   parser.add_argument('query_alignment_fname',help='Alignment file (BAM/CRAM)')
   parser.add_argument('track_directory', help='Directory with BED tracks of pathogenic/actionable cancer loci for grch37/grch38')
   parser.add_argument('output_directory', help='Output directory')
   parser.add_argument('genome_assembly',choices = ['grch37','grch38'], help='Human genome assembly build: grch37 or grch38')
   parser.add_argument('mode',choices=['hereditary','somatic','any'],help="Choice of loci and clinical cancer context (cancer predisposition/tumor sequencing)")
   parser.add_argument('sample_id', help="Sample identifier - prefix for output files")

   docker_image_version = 'sigven/cacao:' + str(cacao_version)
   args = parser.parse_args()
   logger = getlogger('cacao-run')
   logger.info('Start')
   host_directories = verify_arguments(args.query_alignment_fname, args.track_directory, args.output_directory, args.callability_levels_germline, args.callability_levels_somatic, args.force_overwrite, logger)

   input_aln_host = "NA"
   input_aln_index_host = "NA"
   if host_directories['input_alndir'] != 'NA':
      input_aln_host = os.path.join(host_directories['input_alndir'], host_directories['input_aln_basename'])
   if host_directories['input_aln_index_basename'] != 'NA':
      input_aln_index_host = os.path.join(host_directories['input_alndir'], host_directories['input_aln_index_basename']) 

   logger.info('Running cacao workflow - assessment of coverage at actionable and pathogenic loci')
   docker_command = 'docker run --rm -ti -w=/workdir/output -v=' + str(host_directories['track_dir']) + ':/workdir/tracks -v=' + str(host_directories['output_dir']) + ':/workdir/output -v=' + str(input_aln_host) + ':/workdir/query.bam -v=' + str(input_aln_index_host) + ':/workdir/query.bam.bai ' + str(docker_image_version)
   cacao_command = '/cacao.py /workdir/query.bam /workdir/tracks /workdir/output ' + str(args.genome_assembly) + ' ' + str(args.mode) + ' ' + str(args.sample_id) + ' ' + str(args.mapq) + ' ' + str(args.threads) + ' ' + str(args.callability_levels_germline) + ' ' + str(args.callability_levels_somatic)
   run_cacao_cmd = str(docker_command) + ' sh -c "' + str(cacao_command) + '"'
   check_subprocess(run_cacao_cmd)

   logger.info('Finished')

def error_message(message, logger):
   logger.error('')
   logger.error(message)
   logger.error('')
   exit(0)

def warn_message(message, logger):
   logger.warning(message)

def verify_arguments(query_alignment_fname, track_directory, output_dir, callability_levels_germline, callability_levels_somatic, overwrite, logger):
   """
   Function that checks the input files and directories provided by the user and checks for their existence
   """
 
   logger.info('Validating input files and command-line parameters')
   if not re.match(r'^0:[0-9]{1,}:[0-9]{1,}$',callability_levels_germline):
      err_msg = "Wrong fomatting of argument 'callability_levels_germline': It should be a string of three integers separated by ':', starting with 0 (NO_COVERAGE) (see help)"
      error_message(err_msg,logger)
   if not re.match(r'^0:[0-9]{1,}:[0-9]{1,}$',callability_levels_somatic):
      err_msg = "Wrong fomatting of argument 'callability_levels_somatic': It should be a string of three integers separated by ':', starting with 0 (NO_COVERAGE) (see help)"
      error_message(err_msg,logger)

   input_aln_dir = "NA"
   input_aln_basename = "NA"
   input_aln_index_basename = "NA"
   output_dir_full = "NA"
   track_dir_full = "NA"

   ## check the existence of given output folder
   output_dir_full = os.path.abspath(output_dir)
   if not os.path.isdir(output_dir_full):
      err_msg = "Output directory (" + str(output_dir_full) + ") does not exist"
      error_message(err_msg,logger)

   # check the existence of track folder
   track_dir_full = os.path.abspath(track_directory)
   if not os.path.isdir(track_dir_full):
      err_msg = "Track directory (" + str(track_dir_full) + ") does not exist"
      error_message(err_msg,logger)
   
   bam = 0   
   ## check if input vcf exist
   if not query_alignment_fname is None:

      if not query_alignment_fname.endswith('.bam') and not query_alignment_fname.endswith('.cram'):
         err_msg = "Query alignment file (" + str(query_alignment_fname) + ") should have a .bam or .cram suffix (BAM/CRAM)"
         error_message(err_msg,logger)
      
      if not os.path.exists(os.path.abspath(query_alignment_fname)):
         err_msg = "Query alignment file (" + str(query_alignment_fname) + ") does not exist"
         error_message(err_msg,logger)

      if query_alignment_fname.endswith('.bam'):
         bam = 1
         if not os.path.exists(query_alignment_fname + str('.bai')):
            err_msg = "BAM index file (" + str(query_alignment_fname) + ".bai) does not exist"
            error_message(err_msg,logger)
      
      if query_alignment_fname.endswith('.cram'):
         if not os.path.exists(query_alignment_fname + str('.crai')):
            err_msg = "CRAM index file (" + str(query_alignment_fname) + ".crai) does not exist"
            error_message(err_msg,logger)

      input_aln_basename = os.path.basename(str(query_alignment_fname))
      input_aln_dir = os.path.dirname(os.path.abspath(query_alignment_fname))
      input_aln_index_basename = input_aln_basename + '.bai'
      if bam == 0:
         input_aln_index_basename = input_aln_basename + '.crai'

      ##MISSING OVERWRITE CHECK FOR EXISTING OUTPUT FILES
       
   host_directories = {}
   host_directories['input_alndir'] = input_aln_dir
   host_directories['track_dir'] = track_dir_full
   host_directories['output_dir'] = output_dir_full
   host_directories['input_aln_basename'] = input_aln_basename
   host_directories['input_aln_index_basename'] = input_aln_index_basename
   
   return host_directories
   

def check_subprocess(command):
   #print(command)
   try:
      output = subprocess.check_output(str(command), stderr=subprocess.STDOUT, shell=True)
      if len(output) > 0:
         print (str(output.decode()).rstrip())
   except subprocess.CalledProcessError as e:
      print (e.output.decode())
      exit(0)

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

