FROM digitalproteomes/gosu

LABEL maintainer="Patrick Pedrioli"
LABEL Description="A container for the TP data vizualization tools" Version="1.0"

## Let apt-get know we are running in noninteractive mode
ENV DEBIAN_FRONTEND noninteractive

## Install python and R
RUN apt-get update
RUN apt-get install -y python3 r-base

## Copy scripts and data into the container
COPY src /usr/local/src
COPY data /usr/local/data

RUN ln -s /usr/local/src/plot_proteo_fmi_table/plot_proteo_fmi_table.py /usr/local/bin
