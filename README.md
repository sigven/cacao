# *CA*llable *CA*ncer L*O*ci (CACAO)

The *CACAO* framework will provide software and data to assess *sequencing depth for clinically actionable/pathogenic loci in cancer* for a given sequence alignment (BAM/CRAM). Most importantly, the software will pinpoint genomic loci of clinical relevance in cancer that has sufficient **sequencing coverage** for reliable variant calling. In combination with the actual variants that have been identified, it may thus serve to confirm **negative findings**, a matter of significant clinical value that is underappreciated in current cancer sequencing analysis. The specific requirements to denote loci as *callable* (i.e. depth & alignment quality) can be configured by the user, and should thus reflect how the input are used for variant calling (RNA/DNA, germline/somatic calling)

Technically, *CACAO* combines the speed of [mosdepth](https://github.com/brentp/mosdepth) with the powerful [R markdown framework](https://rmarkdown.rstudio.com/) for interactive data reporting. It employs the [Docker](https://www.docker.com) technology for software encapsulation to ease the installation process.

Three clinical genomic tracks in BED format have been created:

* Loci with pathogenic and likely pathogenic variants in protein-coding genes related to cancer predisposition and inherited cancer syndromes (BRCA1, BRCA2, ATM etc.)
	* Variants have been retrieved from [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar) (October 2018 release)
* Loci associated with actionable somatic variants (related to prognosis, diagnosis, or drug sensitivity, e.g. BRAF V600E)
	* Variants have been retrieved from [CIViC](https://civicdb.org) (data obtained October 17th 2018)
* Loci identified as somatic mutational hotspots (i.e. likely driver alterations) in cancer
	* Variants have been retrieved from [cancerhotspots.org](https://www.cancerhotspots.org) (v2)

At each variant identified from the three sources above, we have used a surrounding sequence window of 5bp for which the mean depth is calculated and representing the *loci coverage*.

All three tracks (*hereditary*, *somatic_actionable*, and *somatic_hotspot*) are available for GRCh37 and GRCh38, and there is also tab-separated files that link each locus to its associated
   * variants and phenotypes (ClinVar),
   * clinical evidence items (therapeutic context, evidence level, from CIViC)
   * tumor types (cancerhotspots.org)

## Example report
An [example report](https://folk.uio.no/sigven/test.cacao.grch37.html) from the CACAO workflow showing callable cancer loci in an RNA sequence alignment.

## Getting started

### Installation

* Prerequisites:
	* Make sure that [Docker](https://www.docker.com/) is installed and running
	* The CACAO workflow script requires that Python3 is installed
* Clone the repository `git clone https://github.com/sigven/cacao.git`
* Pull the docker image `docker pull sigven/cacao:1.0.0`

### Usage

Run the CACAO workflow with the `cacao_wflow.py` Python script, which takes the following parameters and options:

cacao_wflow.py [options] <BAM-or-CRAM> <TRACK_DIRECTORY> <OUTPUT_DIR> <GENOME_ASSEMBLY> <CANCER_MODE> <SAMPLE_ID>

	cacao - assessment of sequencing coverage at pathogenic and actionable loci in
	cancer

	positional arguments:
	  query_alignment_fname
	                        Alignment file (BAM/CRAM)
	  track_directory       Directory with BED tracks of pathogenic/actionable
	                        cancer loci for grch37/grch38 (i.e. data/ from repo)
	  output_directory      Output directory
	  {grch37,grch38}       Human genome assembly build: grch37 or grch38
	  {hereditary,somatic,any}
	                        Choice of loci and clinical cancer context (cancer
	                        predisposition/tumor sequencing)
	  sample_id             Sample identifier - prefix for output files

	optional arguments:
	  -h, --help            show this help message and exit
	  --mapq MAPQ           mapping quality threshold (default: 0)
	  --threads THREADS     Number of mosdepth BAM decompression threads. (use 4
	                        or fewer) (default: 0)
	  --callability_levels_germline CALLABILITY_LEVELS_GERMLINE
	                        Simple colon-separated string that defines four levels
	                        of variant callability: NO_COVERAGE (0), LOW_COVERAGE
	                        (0-9), CALLABLE (10-99), HIGH_COVERAGE ( > 100).
	                        Initial value must be 0. (default: 0:10:100)
	  --callability_levels_somatic CALLABILITY_LEVELS_SOMATIC
	                        Simple colon-separated string that defines four levels
	                        of variant callability: NO_COVERAGE (0), LOW_COVERAGE
	                        (0-29), CALLABLE (30-199), HIGH_COVERAGE (> 200).
	                        Initial value must be 0. (default: 0:30:200)
	  --target TARGET       BED file with target regions subject to sequencing
	                        (default: None)
	  --force_overwrite     By default, the script will fail with an error if any
	                        output file already exists. You can force the
	                        overwrite of existing result files by using this flag
	                        (default: False)
	  --version             show program's version number and exit

## Documentation

COMING SOON

## Contact

sigven AT ifi.uio.no
