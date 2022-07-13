FROM rocker/shiny:latest

LABEL maintainer="SoniaGR <s.ruiz@ucl.ac.uk>"

RUN apt-get update && apt-get install -y --no-install-recommends \
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


####### R dependencies #######
RUN R -e 'install.packages(c("AER", "concaveman", "nloptr","units", "pbkrtest", "sf", "car", "terra", "lme4"), repos="http://cran.rstudio.com/", dependencies = T)'

RUN R -e 'install.packages(c("shiny", "shinydashboard", "ggstance", "bslib", "tidyverse", "BiocManager"),repos="http://cran.rstudio.com/", dependencies=T)'

RUN R -e 'BiocManager::install("GenomicRanges")'

RUN apt-get update && apt-get install -y --no-install-recommends \
	libmariadb-dev

RUN R -e 'install.packages(c("DBI", "shinyjs", "shinycssloaders", "shinyBS", "shinydashboard"), repos="http://cran.rstudio.com/", dependencies=T)'

RUN apt-get update && apt-get install -y --no-install-recommends \
	cmake

RUN R -e 'install.packages(c("stringr", "data.table", "ggforce", "gridExtra", "sandwich", "aws.s3"), repos="http://cran.rstudio.com/", dependencies = T)'

RUN R -e 'install.packages(c("DT"), repos="http://cran.rstudio.com/", dependencies = T)'

######## ggtranscript ##########
RUN apt-get update && apt-get install -y --no-install-recommends \
	libharfbuzz-dev \
	libfribidi-dev \
	libfreetype6-dev \
	libpng-dev \
	libtiff5-dev \
	libjpeg-dev


RUN R -e 'install.packages(c("devtools"), repos="http://cran.rstudio.com/", dependencies = T)'
RUN R -e 'devtools::install_github("dzhang32/ggtranscript")'


####### COPY shinyapp #########
COPY --chown=shiny:shiny ./introverse /srv/shiny-server/introverse/


RUN sed -i -e 's/\blisten 3838\b/listen 3848/g' /etc/shiny-server/shiny-server.conf

EXPOSE 3848

CMD ["/usr/bin/shiny-server"]


#CMD ["R", "-e", "shiny::runApp('/srv/shiny-server/introverse')"]
