#!/usr/bin/env python

import csv
import re
import argparse
import os
import subprocess
import logging
import sys
from argparse import RawTextHelpFormatter


CACAO_VERSION = '0.3.1'

def __main__():
   

   program_description = "cacao - assessment of sequencing coverage at pathogenic and actionable loci in cancer"
   program_options = "   --query_aln <QUERY_ALIGNMENT>\n   --track_directory <TRACK_DIR>\n   --output_dir <OUTPUT_DIR>\n   --genome_assembly " + \
      "<grch37|grch38>\n   --sample_id <SAMPLE_ID>\n   --mode <hereditary|somatic|any>"

   parser = argparse.ArgumentParser(description = program_description,
                                    formatter_class=RawTextHelpFormatter, usage="\n  %(prog)s -h [options]\n" + str(program_options))
   parser._action_groups.pop()
   required = parser.add_argument_group('Required arguments')
   optional = parser.add_argument_group('Optional arguments')

   optional.add_argument('--mapq', dest = "mapq", default = 0, type=int, help='mapping quality threshold')
   optional.add_argument('--threads', dest = "threads", default = 0, type=int, help='Number of mosdepth BAM decompression threads. (use 4 or fewer)')
   optional.add_argument('--callability_levels_germline', dest="callability_levels_germline",default="0:10:100", help="Simple colon-separated string that defines four levels of variant callability: NO_COVERAGE (0), LOW_COVERAGE (1-9), CALLABLE (10-99), HIGH_COVERAGE (>= 100). Initial value must be 0.")
   optional.add_argument('--callability_levels_somatic', dest="callability_levels_somatic",default="0:30:200", help="Simple colon-separated string that defines four levels of variant callability: NO_COVERAGE (0), LOW_COVERAGE (1-29), CALLABLE (30-199), HIGH_COVERAGE (>= 200). Initial value must be 0.")
   optional.add_argument('--query_target', dest = "query_target", help="BED file with genome target regions subject to sequencing")
   optional.add_argument('--is_rna',action = "store_true",help='By default, CACAO assumes a BAM file coming from DNA sequencing. Indicate whether the alignment is from RNA-seq with this flag')
   optional.add_argument('--force_overwrite', action = "store_true", help='By default, the script will fail with an error if any output file already exists. You can force the overwrite of existing result files by using this flag')
   #parser.add_argument('--no-docker', action='store_true', dest='no_docker', default=False, help='Run the cacao workflow in a non-Docker mode (see install_no_docker/ folder for instructions')
   optional.add_argument('--version', action='version', version='%(prog)s ' + str(CACAO_VERSION))
   required.add_argument('--query_aln',help='Query alignment file (BAM/CRAM)')
   required.add_argument('--track_dir', help='Directory with BED tracks of pathogenic/actionable cancer loci for grch37/grch38')
   required.add_argument('--output_dir', help='Output directory')
   required.add_argument('--genome_assembly',choices = ['grch37','grch38'], help='Human genome assembly build: grch37 or grch38')
   required.add_argument('--mode',choices=['hereditary','somatic','any'],help="Choice of loci and clinical cancer context (cancer predisposition/tumor sequencing)")
   required.add_argument('--sample_id', help="Sample identifier - prefix for output files")
   
   args = parser.parse_args()
   arg_dict = vars(args)
   
   docker_image_version = 'sigven/cacao:' + str(CACAO_VERSION)
   #args = parser.parse_args()
   logger = getlogger('cacao-run')
   logger.info('Start')
   overwrite = 0
   if args.force_overwrite is True:
      overwrite = 1


   ## Required arguments
   ## Check the existence of required arguments
   if arg_dict['track_dir'] is None or not os.path.exists(arg_dict['track_dir']):
      err_msg = "Required argument --track_dir has no/undefined value (" + str(arg_dict['track_dir']) + "). Type cacao_wflow.py --help to view all options and required arguments"
      error_message(err_msg,logger)
   
   if arg_dict['output_dir'] is None or not os.path.exists(arg_dict['output_dir']):
      err_msg = "Required argument --output_dir has no/undefined value (" + str(arg_dict['output_dir']) + "). Type cacao_wflow.py --help to view all options and required arguments"
      error_message(err_msg,logger)
   
   if arg_dict['genome_assembly'] is None:
      err_msg = "Required argument --genome_assembly has no/undefined value (" + str(arg_dict['genome_assembly']) + "). Type cacao_wflow.py --help to view all options and required arguments"
      error_message(err_msg,logger)
   
   if arg_dict['query_aln'] is None:
      err_msg = "Required argument --query_aln has no/undefined value (" + str(arg_dict['query_aln']) + "). Type cacao_wflow.py --help to view all options and required arguments"
      error_message(err_msg,logger)
   
   if arg_dict['sample_id'] is None:
      err_msg = "Required argument --sample_id has no/undefined value (" + str(arg_dict['sample_id']) + "). Type cacao_wflow.py --help to view all options and required arguments"
      error_message(err_msg,logger)
   
   if arg_dict['mode'] is None:
      err_msg = "Required argument --mode has no/undefined value (" + str(arg_dict['mode']) + "). Type cacao_wflow.py --help to view all options and required arguments"
      error_message(err_msg,logger)

   host_directories = verify_input_files(arg_dict, logger)

   input_aln_host = "NA"
   input_aln_index_host = "NA"
   input_target_host = "NA"
   if host_directories['input_alndir'] != 'NA':
      input_aln_host = os.path.join(host_directories['input_alndir'], host_directories['input_aln_basename'])
   if host_directories['input_aln_index_basename'] != 'NA':
      input_aln_index_host = os.path.join(host_directories['input_alndir'], host_directories['input_aln_index_basename']) 

   logger.info('Running cacao workflow - assessment of coverage at actionable and pathogenic loci')
   docker_command = 'docker run --rm -ti -w=/workdir/output -v=' + str(host_directories['track_dir']) + ':/workdir/tracks -v=' + \
      str(host_directories['output_dir']) + ':/workdir/output -v=' + str(input_aln_host) + ':/workdir/query.bam -v=' + \
      str(input_aln_index_host) + ':/workdir/query.bam.bai ' + str(docker_image_version)
   if host_directories['input_target_dir'] != "NA" and host_directories['input_target_basename'] != "NA":
      input_target_host = os.path.join(host_directories['input_target_dir'], host_directories['input_target_basename'])
      docker_command = 'docker run --rm -ti -w=/workdir/output -v=' + str(host_directories['track_dir']) + ':/workdir/tracks -v=' + \
         str(host_directories['output_dir']) + ':/workdir/output -v=' + str(input_aln_host) + ':/workdir/query.bam -v=' + \
         str(input_aln_index_host) + ':/workdir/query.bam.bai -v' + str(input_target_host) + ':/workdir/query_target.bed ' + str(docker_image_version)
   cacao_command = '/cacao.py /workdir/query.bam /workdir/query_target.bed ' + str(host_directories['input_aln_basename']) + ' ' + \
      str(host_directories['input_target_basename']) + ' /workdir/tracks /workdir/output ' + str(arg_dict['genome_assembly']) + ' ' + \
      str(arg_dict['mode']) + ' ' + str(arg_dict['sample_id']) + ' ' + str(arg_dict['mapq']) + ' ' + str(arg_dict['threads']) + ' ' + \
      str(arg_dict['callability_levels_germline']) + ' ' + str(arg_dict['callability_levels_somatic'])
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

def verify_input_files(arg_dict, logger):
   """
   Function that checks the input files and directories provided by the user and checks for their existence
   """
 
   logger.info('Validating input files and command-line parameters')
   if not re.match(r'^0:[0-9]{1,}:[0-9]{1,}$',arg_dict['callability_levels_germline']):
      err_msg = "Wrong fomatting of argument 'callability_levels_germline': It should be a string of three integers separated by ':', starting with 0 (NO_COVERAGE) (see help)"
      error_message(err_msg,logger)
   if not re.match(r'^0:[0-9]{1,}:[0-9]{1,}$',arg_dict['callability_levels_somatic']):
      err_msg = "Wrong fomatting of argument 'callability_levels_somatic': It should be a string of three integers separated by ':', starting with 0 (NO_COVERAGE) (see help)"
      error_message(err_msg,logger)

   input_aln_dir = "NA"
   input_aln_basename = "NA"
   input_target_dir = "NA"
   input_target_basename = "NA"
   input_aln_index_basename = "NA"
   output_dir_full = "NA"
   track_dir_full = "NA"

   overwrite = 0
   if arg_dict['force_overwrite'] is True:
      overwrite = 1

   ## check the existence of given output folder
   output_dir_full = os.path.abspath(arg_dict['output_dir'])
   if not os.path.isdir(output_dir_full):
      err_msg = "Output directory (" + str(output_dir_full) + ") does not exist"
      error_message(err_msg,logger)

   # check the existence of track folder
   track_dir_full = os.path.abspath(arg_dict['track_dir'])
   if not os.path.isdir(track_dir_full):
      err_msg = "Track directory with actionable/pathogenic loci (" + str(track_dir_full) + ") does not exist"
      error_message(err_msg,logger)
   
   # check if target file (BED) exists
   if not arg_dict['query_target'] is None:

      if not arg_dict['query_target'].endswith('.bed'):
         err_msg = "Query target file (" + str(arg_dict['query_target']) + ") should have a .bed suffix (BED format)"
         error_message(err_msg,logger)

      if not os.path.exists(os.path.abspath(arg_dict['query_target'])):
         err_msg = "Query target file (" + str(arg_dict['query_target']) + ") does not exist"
         error_message(err_msg,logger)
      
      input_target_basename = os.path.basename(str(arg_dict['query_target']))
      input_target_dir = os.path.dirname(os.path.abspath(arg_dict['query_target']))

   ## if output html exist and overwrite not set
   output_html = os.path.join(str(arg_dict['output_dir']), str(arg_dict['sample_id']) + '_' + str(arg_dict['genome_assembly']) + '_coverage_cacao.html')
   if os.path.exists(output_html) and overwrite == 0:
      err_msg = "Output files (e.g. " + str(output_html) + ") already exist - please specify different sample_id or add option --force_overwrite"
      error_message(err_msg,logger)

   bam = 0 
   ## check if input BAM/CRAM exists
   if not arg_dict['query_aln'] is None:

      if not arg_dict['query_aln'].endswith('.bam') and not arg_dict['query_aln'].endswith('.cram'):
         err_msg = "Query alignment file (" + str(arg_dict['query_aln']) + ") should have a .bam or .cram suffix (BAM/CRAM)"
         error_message(err_msg,logger)
      
      if not os.path.exists(os.path.abspath(arg_dict['query_aln'])):
         err_msg = "Query alignment file (" + str(arg_dict['query_aln']) + ") does not exist"
         error_message(err_msg,logger)

      if arg_dict['query_aln'].endswith('.bam'):
         bam = 1
         if not os.path.exists(arg_dict['query_aln'] + str('.bai')):
            err_msg = "BAM index file (" + str(arg_dict['query_aln']) + ".bai) does not exist"
            error_message(err_msg,logger)
      
      if arg_dict['query_aln'].endswith('.cram'):
         if not os.path.exists(arg_dict['query_aln'] + str('.crai')):
            err_msg = "CRAM index file (" + str(arg_dict['query_aln']) + ".crai) does not exist"
            error_message(err_msg,logger)

      input_aln_basename = os.path.basename(str(arg_dict['query_aln']))
      input_aln_dir = os.path.dirname(os.path.abspath(arg_dict['query_aln']))
      input_aln_index_basename = input_aln_basename + '.bai'
      if bam == 0:
         input_aln_index_basename = input_aln_basename + '.crai'

      ##MISSING OVERWRITE CHECK FOR EXISTING OUTPUT FILES
       
   host_directories = {}
   host_directories['track_dir'] = track_dir_full
   host_directories['output_dir'] = output_dir_full
   host_directories['input_alndir'] = input_aln_dir
   host_directories['input_aln_basename'] = input_aln_basename
   host_directories['input_aln_index_basename'] = input_aln_index_basename
   host_directories['input_target_dir'] = input_target_dir
   host_directories['input_target_basename'] = input_target_basename
   
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

