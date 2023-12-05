FROM rocker/tidyverse
# Builds on Ubuntu LTS

ENV DEBIAN_FRONTEND noninteractive

WORKDIR /home

RUN mkdir /home/scripts/ /home/01_genotypes/ /home/02_LD/ /home/03_ROH/ \
/home/04_PCA/ /home/05_BLUPF90/ /home/scenario_1/ /home/scenario_2/ \
/home/scenario_3/ /home/scenario_4/

RUN apt-get -y update
RUN apt-get -y install plink1.9
RUN apt-get -y install python3
RUN apt-get -y install python3-pip
RUN apt-get -y install less
RUN apt-get -y install bsdmainutils
RUN pip install pandas==2.0.3
RUN apt-get install nano -y
RUN mv /usr/bin/plink1.9 /usr/bin/plink

ADD scripts/ /home/scripts/
RUN mv /home/scripts/extract_ped.sh /home/05_BLUPF90/extract_ped.sh
RUN mv /home/scripts/prepare_snp_file.sh /home/05_BLUPF90/prepare_snp_file.sh
RUN mv /home/scripts/renumf90 /home/05_BLUPF90/renumf90
RUN mv /home/scripts/blupf90 /home/05_BLUPF90/blupf90
RUN mv /home/scripts/renum.txt /home/05_BLUPF90/renum.txt

RUN /usr/local/bin/R --no-restore --file=scripts/R_packages.R > scripts/R_packages.txt
RUN chmod 755 scripts/run_pipeline.sh
RUN chmod 755 scripts/run_simulation.sh 
RUN chmod 755 05_BLUPF90/renumf90
RUN chmod 755 05_BLUPF90/blupf90

CMD /bin/bash