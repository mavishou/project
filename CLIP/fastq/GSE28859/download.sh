#!/bin/bash

while read s
do
	wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP006/SRP006474/${s}/${s}.sra
done
