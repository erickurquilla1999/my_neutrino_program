FROM ubuntu:22.04
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update
RUN apt-get install -y python3 python3-pip git julia
RUN pip3 install numpy matplotlib h5py scipy
ENV USER=jenkins
ENV LOGNAME=jenkins
