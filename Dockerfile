FROM rocker/shiny:latest
MAINTAINER Markus List <markus.list@mpi-inf.mpg.de>

#install system packages
RUN apt-get update && apt-get install -y \
libxml2-dev \ 
libcurl4-gnutls-dev \
libssl-dev 

#install R packages
ADD install.R install.R
RUN R -e "source('install.R')"

#copy shiny app to work-dir
WORKDIR /srv/
RUN mkdir pathclass
ADD . pathclass

#update shiny server conf and configure it to run pathclass in single app mode
RUN sed -i 's/site_dir \/srv\/shiny-server;/app_dir \/srv\/pathclass;/g' /etc/shiny-server/shiny-server.conf

