FROM ubuntu:22.04
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update
RUN apt-get install -y python python3 python3-pip gfortran build-essential libhdf5-openmpi-dev openmpi-bin pkg-config libopenmpi-dev openmpi-bin libblas-dev liblapack-dev libpnetcdf-dev git julia
RUN pip3 install numpy matplotlib h5py scipy
ENV USER=jenkins
ENV LOGNAME=jenkins
