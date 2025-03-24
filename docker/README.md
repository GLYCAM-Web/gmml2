A pre-built image exists on dockerhub here:
https://hub.docker.com/repository/docker/olyball/gmml2/general

To build the image yourself using the files in the already cloned gmml2 repository: 
cd gmml2/
docker build -t gmml2 .

To build by cloning the gmml2 github repo during the build:
cd gmml2/docker/dockerhub/
docker build -t gmml2 .

Usage instructions: gmml2/docker/README_USAGE.md
