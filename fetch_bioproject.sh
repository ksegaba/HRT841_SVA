#!/bin/bash
# Description: Retrieve Bioproject IDs from NCBI's SRA database
# Author: Kenia Segura Abá
# Code source: https://www.biostars.org/p/444581/ 

while read i;
do
	VAR="$(esearch -db sra -query ${i} < /dev/null | efetch  -format runinfo -mode xml | xtract -pattern SraRunInfo -element BioProject)"
	echo "${i}, ${VAR}" > bioproject_ids.txt
done <sra.txt
