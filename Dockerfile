FROM rocker/shiny:3.6.3
MAINTAINER Markus List <markus.list@wzw.tum.de>

#install system packages
RUN apt-get update && apt-get install -y \
libxml2-dev \ 
libcurl4-gnutls-dev \
libssl-dev 

#copy shiny app to work-dir
RUN mkdir /srv/pathclass
WORKDIR /srv/pathclass
ADD . .

#install R packages
ENV RENV_VERSION 0.10.0-2
RUN R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org')); \ 
          remotes::install_github('rstudio/renv@${RENV_VERSION}'); \
          renv::restore()"

#update shiny server conf and configure it to run pathclass in single app mode
RUN sed -i 's/site_dir \/srv\/shiny-server;/app_dir \/srv\/pathclass;/g' /etc/shiny-server/shiny-server.conf

