FROM rocker/binder:3.6.3

ARG NB_USER
ARG NB_UID

USER root

COPY . ${HOME}
RUN chown -R ${NB_USER} ${HOME}

RUN pwd
RUN ls -la
RUN ls -la .binder 

#RUN R -e 'source("~/.binder/installPackages.R"); installPackages("~/.binder/r-requirements.dcf")'
RUN R -e 'devtools::install_deps(".", dependencies=TRUE)'
RUN R -e 'devtools::install(".")'

USER ${NB_USER}

