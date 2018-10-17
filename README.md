# *CA*llable *CA*ncer L*O*ci (CACAO)

The *CACAO* framework will provide software and data to assess *sequencing depth for clinically actionable/pathogenic loci in cancer* for a given sequence alignment (BAM/CRAM). Most importantly, the software will pinpoint genomic loci of clinical relevance in cancer that has sufficient **sequencing coverage** for reliable variant calling. In combination with the actual variants that have been identified, it may thus serve to confirm **negative findings**. The specific requirements to denote loci as *callable* (i.e. depth & alignment quality) can be configured by the user, and should thus reflect how the input are used for variant calling (RNA/DNA, germline/somatic calling)

Technically, *CACAO* combines the speed of [mosdepth](https://github.com) with the powerful [R markdown framework](https://rmarkdown.rstudio.com/) for interactive data reporting. It employs the [Docker](https://www.docker.com) technology for software encapsulation to ease the installation

Three clinical genomic tracks in BED format have been created:

* Loci with pathogenic and likely pathogenic variants in protein-coding genes related to cancer predisposition (BRCA1, BRCA2, ATM etc.)
	* Variants have been retrieved from [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar) (October 2018 release)
* Loci associated with actionable somatic variants (related to prognosis, diagnosis, or drug sensitivity, e.g. BRAF V600E)
	* Variants have been retrieved from [CIViC](https://civicdb.org) (data obtained September 22nd 2018)
* Loci identified as somatic mutational hotspots in cancer
	* Variants have been retrieved from [cancerhotspots.org](https://cancerhotspots.org) (v2)

At each variant identified from the three sources above, we have used a surrounding sequence window of 5bp to denote the *genomic locus*.

All three tracks (*hereditary*, *somatic_actionable*, and *somatic_hotspot*) are available for GRCh37 and GRCh38, and there is also tab-separated files that link each locus to its associated
   * variants and phenotypes (ClinVar),
   * clinical evidence items (therapeutic context, evidence level, from CIViC)
   * tumor types (cancerhotspots.org)
