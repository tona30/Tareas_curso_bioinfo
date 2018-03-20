#!/bin/bash

##Este script contiene un jemplo de busqueda en BLAST de una secuencia de proteina.

# ubicar la carpeta de trabajo
/home/tona/Escritorio/Blastej

# bajar secuencia de prion human
wget http://www.uniprot.org/uniprot/P04156.fasta

# bajar la base de datos de zebrafish
curl -O ftp://ftp.ncbi.nih.gov/refseq/D_rerio/mRNA_Prot/zebrafish.1.protein.faa.gz

# descomprimir la carpeta
gunzip zebrafish.1.protein.faa.gz

# prepara la base de datos para la busqueda, montar un volumen al directorio deseado
sudo docker run -v /home/tona/Escritorio/Blastej:/data/ biocontainers/blast makeblastdb -in zebrafish.1.protein.faa -dbtype prot 

# prepara alineamientos finales y monta un volumen al directorio deseado
sudo docker run -v /home/tona/Escritorio/Blastej:/data/ biocontainers/blast blastp -query P04156.fasta -db zebrafish.1.protein.faa -out results.txt

