FROM ubuntu:18.04

LABEL maintainer="Patrick Pedrioli"
LABEL Description="A container for the TP data vizualization tools" Version="1.0"

## Let apt-get know we are running in noninteractive mode
ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update
RUN apt-get install -y python3 r-base

COPY src /local/src
COPY data /local/data
