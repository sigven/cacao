#cp -L /Users/sigven/research/software/clinseqcov/code/clinseqcov.py .
#cp -L /Users/sigven/research/software/clinseqcov/output/cancer_hereditary_pathogenic_loci.g* clinical_tracks/
#cp -L /Users/sigven/research/software/clinseqcov/output/cancer_somatic_actionable_loci.g* clinical_tracks/
#cp -L /Users/sigven/research/software/clinseqcov/code/clinseqcov.R .
#cp -L /Users/sigven/research/software/clinseqcov/code/clinseqcov.Rmd .
#tar cvzf clinical_tracks.tgz clinical_tracks/
echo "Build the Docker Image"
TAG=`date "+%Y%m%d"`
docker build -t sigven/cacao:$TAG --rm=true .
