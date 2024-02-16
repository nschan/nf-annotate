from continuumio/miniconda3:latest

RUN apt-get --allow-releaseinfo-change update && apt-get install -y libgomp1

RUN git clone https://github.com/thecgs/genblastG_extension.git

RUN mv genblastG_extension /opt/

RUN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.14.0+-x64-linux.tar.gz && tar -xf ncbi-blast-2.14.0+-x64-linux.tar.gz && rm https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.14.0+-x64-linux.tar.gz

RUN mkdir -p /opt/bin/ncbi-blast/

RUN mv ncbi-blast-2.14.0+/bin/* /opt/bin/ncbi-blast/

RUN export PATH="/opt/genblastG_extensions/:/opt/bin/ncbi-blast/:$PATH"