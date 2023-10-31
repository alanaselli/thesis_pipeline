FROM rocker/tidyverse
# Builds on Ubuntu LTS

ENV DEBIAN_FRONTEND noninteractive

WORKDIR /home

RUN mkdir /home/01_genotypes/ /home/02_LD/ /home/03_ROH/ \
/home/04_PCA/ /home/scripts/

RUN apt-get -y update
RUN apt-get -y install plink1.9
RUN apt-get -y install python3
RUN apt-get -y install python3-pip
RUN apt-get -y install less
RUN pip install pandas==2.0.3
RUN apt-get install nano -y
RUN mv /usr/bin/plink1.9 /usr/bin/plink

ADD scripts/ /home/scripts/

RUN /usr/local/bin/R --no-restore --file=R_packages.R > R_packages.txt
RUN chmod 755 run_pipeline.sh

CMD /bin/bash