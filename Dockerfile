############################################################
# Dockerfile to build bisukit images
# Based on Ubuntu
############################################################


# Download ubuntu from the docker hub
FROM ubuntu:latest

# LABEL maintainer="jawon@tacc.utexas.edu"

# Set up environments
RUN apt-get update;\
    apt-get -y install build-essential;\
    apt-get -y install gfortran;\
    apt-get -y install fort77;\
    apt-get -y install xorg-dev;\ 
    apt-get -y install liblzma-dev;\
    apt-get -y install libblas-dev;\
    apt-get -y install gcc-multilib;\ 
    apt-get -y install gobjc++;\
    apt-get -y install xorg openbox xvfb xvfb-run;\
    apt-get -y install libreadline-dev;\
    apt-get -y install python2.7 python2.7-dev;\
    apt-get -y install python-setuptools python-dev build-essential;\
    apt-get -y install wget git curl apache2 python-pip 

# Make directory
RUN mkdir -p /home/bisukit
WORKDIR /home/bisukit
#COPY . /home/bisukit

# Executable
#ADD https://github.com/al2na/methylKit/archive/v0.9.5.tar.gz /home/bisukit
ADD https://github.com/wonaya/BisuKit/raw/master/R.tar.gz /home/bisukit
ADD https://github.com/wonaya/BisuKit/raw/master/lib.tar.gz /home/bisukit
RUN cd /home/bisukit/; curl https://cran.r-project.org/src/base/R-2/R-2.15.1.tar.gz | tar xz;\
  cd R-2.15.1; chmod +x configure; ./configure --enable-R-shlib --enable-R-static-lib --enable-BLAS-shlib ; make ; make install 

ENV PYTHONPATH $PYTHONPATH:/usr/lib/python2.7

#Example files
ADD https://raw.githubusercontent.com/wonaya/BisuKit/master/CpG_OB_sample1.txt /home/bisukit
ADD https://raw.githubusercontent.com/wonaya/BisuKit/master/CpG_OB_sample2.txt /home/bisukit
ADD https://raw.githubusercontent.com/wonaya/BisuKit/master/CpG_OT_sample1.txt /home/bisukit
ADD https://raw.githubusercontent.com/wonaya/BisuKit/master/CpG_OT_sample2.txt /home/bisukit
ADD https://raw.githubusercontent.com/wonaya/BisuKit/master/TAIR_4.fa /home/bisukit

#Get R packages
RUN pwd
#COPY lib.tar.gz /home/bisukit/.
#COPY R /home/bisukit/.
RUN tar -zxf /home/bisukit/lib.tar.gz
RUN tar -zxf /home/bisukit/R.tar.gz
#RUN tar -zxvf /home/bisukit/v0.9.5.tar.gz
RUN ls /home/bisukit
#RUN Rscript -e 'source("http://bioconductor.org/biocLite.R")' -e 'biocLite("GenomicRanges")'
#RUN cp -r /home/bisukit/methylKit-0.9.5 /usr/local/lib/R/library/.
RUN cp -r /home/bisukit/lib/tempRlibs/* /usr/local/lib/R/library/.
ENV R_LIBS /home/bisukit/lib/tempRlibs:/usr/local/lib/R/library/:$R_LIBS
ENV PATH /usr/bin:$PATH
ENV PYTHON_VERSION 2.7.12

#RUN tar -zxvf methylKit-0.9.5.tar.gz
#RUN cp -r methylKit-0.9.5 /usr/local/lib/R/library/.
#RUN wget https://cran.r-project.org/src/contrib/chron_2.3-52.tar.gz
#RUN R CMD INSTALL chron_2.3-52.tar.gz
#RUN wget https://cran.r-project.org/src/contrib/iterators_1.0.9.tar.gz
#RUN R CMD INSTALL iterators_1.0.9.tar.gz
#RUN wget https://cran.r-project.org/src/contrib/foreach_1.4.4.tar.gz
#RUN R CMD INSTALL foreach_1.4.4.tar.gz
#RUN wget https://cran.r-project.org/src/contrib/doMC_1.3.5.tar.gz
#RUN R CMD INSTALL doMC_1.3.5.tar.gz
#RUN wget https://cran.r-project.org/src/contrib/Archive/data.table/data.table_1.9.6.tar.gz
#RUN R CMD INSTALL data.table_1.9.6.tar.gz
#RUN wget https://cran.r-project.org/src/contrib/bigmemory.sri_0.1.3.tar.gz
#RUN R CMD INSTALL bigmemory.sri_0.1.3.tar.gz
#RUN wget https://cran.r-project.org/src/contrib/Archive/bigmemory/bigmemory_4.3.0.tar.gz
#RUN R CMD INSTALL bigmemory_4.3.0.tar.gz 
#RUN pip install --upgrade pip
RUN pip install -U "matplotlib==2.0.2" "psutil" "rpy2==2.2.1"
#RUN pip install -U "matplotlib==2.0.2" "psutil" 
ADD https://raw.githubusercontent.com/wonaya/BisuKit/master/bisukit.py /home/bisukit/

ENV LD_LIBRARY_PATH /usr/local/lib/R/lib:$LD_LIBRARY_PATH
CMD export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib/R/lib
CMD export R_LIBS=/home/bisukit/lib/tempRlibs:$R_LIBS
CMD export PYTHONPATH=/home/bisukit/lib/python:$PYTHONPATH
#RUN cd /home/bisukit && /usr/bin/python2.7 bisukit.py --name1 sample1 --name2 sample2 --ot1 CpG_OT_sample1.txt --ot2 CpG_OT_sample2.txt --ob1 CpG_OB_sample1.txt --ob2 CpG_OB_sample2.txt --context CpG --genome TAIR_4.fa --specie TAIR --cores 4

ENTRYPOINT ["/usr/bin/python2.7", "/home/bisukit/bisukit.py"]
#
