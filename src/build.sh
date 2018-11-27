echo "Build the Docker Image"
TAG=`date "+%Y%m%d"`
docker build -t sigven/cacao:$TAG --rm=true .
