FROM rocker/shiny:4.2.1

LABEL maintainer="SoniaGR <s.ruiz@ucl.ac.uk>"

RUN sudo apt-get update 

####### System libraries #######

RUN sudo apt-get install -y --no-install-recommends \
    sudo \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libssh2-1-dev \
    git-core \
    curl \
    libsodium-dev \
    libxml2-dev \
    libicu-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN apt-get update && apt-get install -y --no-install-recommends \
	libmariadb-dev

RUN sudo apt-get update && apt-get install -y --no-install-recommends \
	libfreetype6-dev \
	libpng-dev \
	libtiff5-dev \
	libjpeg-dev

RUN apt-get update && apt-get install -y --no-install-recommends \
  libudunits2-dev \
  libharfbuzz-dev \
  libfribidi-dev 
  
RUN apt-get update && apt-get install -y --no-install-recommends \
	cmake
	


####### R dependencies #######

RUN R -e 'install.packages(c("AER", "concaveman", "nloptr","units", "pbkrtest", "sf", "car", "terra", "lme4"), repos="http://cran.rstudio.com/", dependencies = T)'

RUN R -e 'install.packages(c("shiny", "shinydashboard", "ggstance", "bslib", "tidyverse", "BiocManager"),repos="http://cran.rstudio.com/", dependencies=T)'

RUN R -e 'BiocManager::install("GenomicRanges")'

RUN R -e 'install.packages(c("DBI", "shinyjs", "shinycssloaders", "shinyBS", "shinydashboard"), repos="http://cran.rstudio.com/", dependencies=T)'

RUN R -e 'install.packages(c("stringr", "data.table", "ggforce", "gridExtra", "sandwich"), repos="http://cran.rstudio.com/", dependencies = T)'

RUN R -e 'install.packages(c("DT"), repos="http://cran.rstudio.com/", dependencies = T)'

RUN R -e 'install.packages(c("shinylogs"), repos="http://cran.rstudio.com/", dependencies = T)'

######## ggtranscript ##########

RUN R -e 'install.packages(c("devtools"), repos="http://cran.rstudio.com/", dependencies = T)'
RUN R -e 'devtools::install_github("dzhang32/ggtranscript")'

####### other R libraries ######

RUN R -e 'install.packages(c("pbkrtest"), source="https://cran.rstudio.com/", dependencies = T)'
RUN R -e 'install.packages(c("ggpubr"), repos="http://cran.rstudio.com/", dependencies = T)'
RUN R -e 'install.packages(c("ggrepel"), repos="http://cran.rstudio.com/", dependencies = T)'

####### COPY shinyapp #########

COPY Rprofile.site /usr/lib/R/etc/
RUN mkdir /root/introverse
COPY ./introverse /root/introverse


####### EXPOSE #########

EXPOSE 3838

CMD ["R", "-e", "shiny::runApp('/root/introverse', host='0.0.0.0', port=3838)"]
