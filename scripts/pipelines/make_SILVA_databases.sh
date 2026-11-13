#!/bin/bash

# Script to set up a SILVA SSU and LSU NR99 BLASTand kraken2 DB

# Usage: make_SILVA_databases.sh path_to_new_kraken2_DB_directory path_to_new_BLAST_DB_directory
# Dependencies:
#   Python
#   BLAST
#   Kraken2
module load blast
module load kraken

# Set variable to FULL directory path for generated SILVA kraken2 and BLAST DB
# Note: kraken2 won't work if the DB is moved after building
kraken2DB=/mnt/cluster_storage/data/ref/metagenomics_vs_totalRNA/kraken_db #$1
blastDB=/mnt/cluster_storage/data/ref/metagenomics_vs_totalRNA/blast_db #$2

# ended up trying slightly newer version
VER="138.2"

cwd=$(pwd)
mkdir -p make_SILVA_databases/
cd make_SILVA_databases/

# Get SILVA fasta and taxonomy files
#wget      https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.1_SSURef_NR99_tax_silva_trunc.fasta.gz
wget https://www.arb-silva.de/fileadmin/silva_databases/release_$VER/Exports/SILVA_${VER}_SSURef_NR99_tax_silva_trunc.fasta.gz

#wget      https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.1_LSURef_NR99_tax_silva_trunc.fasta.gz
wget https://www.arb-silva.de/fileadmin/silva_databases/release_$VER/Exports/SILVA_${VER}_LSURef_NR99_tax_silva_trunc.fasta.gz

#wget      https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/taxonomy/taxmap_slv_lsu_ref_nr_138.1.txt.gz
wget https://www.arb-silva.de/fileadmin/silva_databases/release_$VER/Exports/taxonomy/taxmap_slv_lsu_ref_nr_$VER.txt.gz

#wget      https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/taxonomy/taxmap_slv_ssu_ref_nr_138.1.txt.gz
wget https://www.arb-silva.de/fileadmin/silva_databases/release_$VER/Exports/taxonomy/taxmap_slv_ssu_ref_nr_$VER.txt.gz

#wget      https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/taxonomy/tax_slv_ssu_138.1.txt.gz
wget https://www.arb-silva.de/fileadmin/silva_databases/release_$VER/Exports/taxonomy/tax_slv_ssu_$VER.txt.gz

gunzip *.gz

# Get NCBI taxonomy files
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip
unzip taxdmp.zip
mv names.dmp names_ncbi.dmp

# Edit SILVA taxonomy file:
# SILVA taxonomy file needs to have two columns, first one with accession number
# and second one with taxonomy path
sed -e 's/\t/./' -e 's/\t/./' taxmap_slv_ssu_ref_nr_$VER.txt > taxmap_slv_ssu_lsu_ref_nr_$VER.txt
tail -n +2 taxmap_slv_ssu_lsu_ref_nr_$VER.txt | cut -f1 | sed s'/\./\\\./g' > grep_list.txt
sed -e 's/\t/./' -e 's/\t/./' taxmap_slv_lsu_ref_nr_$VER.txt | tail -n +2 \
| grep -v -f grep_list.txt >> taxmap_slv_ssu_lsu_ref_nr_$VER.txt
sed -i "s/ <[a-zA-Z -,.&:'0-9]*>//" taxmap_slv_ssu_lsu_ref_nr_$VER.txt
# 1. concatenate ssu and lsu taxmap files (with access id and start and stop merged)
#    into one ssu_lsu file with one header
#	  	 note: as of 04 Sep 2020, the available SILVA LSU and SSU taxmap files contain
#   	 duplicate accession IDs with different taxids in the SSU and LSU files. We're
#	  	 going to use all accession IDs from the SSU file, check which additional ones
# 	  are in the LSU file (about 23,000 are not in the SSU file) and just take these
# 	  extra ones from the LSU taxmap file to not overwrite SSU taxids with LSU taxids:
# 2. removes <genus>, <family> etc. from taxonomic paths (otherwise the SILVA
#    taxonomy won't match the NCBI taxonomy)


# As of 04 Sep 2020, the available SILVA LSU and SSU version contain duplicate
# sequences. I contacted the SILVA support, who said that should not have
# happened and that we can remove these duplicates. So we're going to remove
# them from the SSU file with a chunk of code that I got from the support
# (filter-fasta.awk is a separate subscript):
RE="^>(KY7649(2[1278]|58)|LNRQ01000003)\\\\."
cat SILVA_${VER}_SSURef_NR99_tax_silva_trunc.fasta \
| filter-fasta.awk -v expression=$RE -v printMatch=0 \
> SILVA_${VER}_SSURef_NR99_tax_silva_trunc_filtered.fasta

# Prepare the SSU and LSU SILVA DB
cat SILVA_${VER}_SSURef_NR99_tax_silva_trunc_filtered.fasta \
SILVA_${VER}_LSURef_NR99_tax_silva_trunc.fasta \
> SILVA_${VER}_SSU_LSURef_NR99_tax_silva_trunc.fasta

# Turn fasta into tab delimited file and edit so that accession numbers can be merged later
${cwd}/subscripts/fasta_to_tab SILVA_${VER}_SSU_LSURef_NR99_tax_silva_trunc.fasta \
| sed 's/ /\t/' > SILVA_${VER}_SSU_LSURef_NR99_tax_silva_trunc_edited.tab

# Edit names_ncbi.dmp file into a file only containing scientific names
sed "s/ <[a-zA-Z -,.&:'0-9]*>//g" names_ncbi.dmp | grep 'scientific name' \
| cut -f 1,3 | awk -F $'\t' ' { t = $1; $1 = $2; $2 = t; print; } ' OFS=$'\t' \
| grep -v 'environmental\|uncultured\|unidentified\|metagenome' \
> NCBI_staxids_scientific.txt
  # 1. Remove <genus>, <family>, and other strings in <> brackets from taxonomy
  #    (otherwise the NCBI taxonomy won't match the SILVA taxonomy)
  # 2. Extract only scientific names and staxids
  # 3. Cut out columns we need
  # 4. Invert columns for script to work
  # 5. Removes lines containing "environmental", "uncultured", "unidentified",
  #    and "metagenome"
    # Needed because
      # SILVA taxonomy can have the same taxonomic ranks (e.g., "environmental
      # sample") for different higher ranks (e.g., "nematode; environmental
      # sample" and "bacteria;environmental sample"), which would, however, be
      # assigned to the same staxid because the lower rank "environmental sample"
      # is similar
      # NCBI taxonomy can have different staxids for the same taxonomic name,
      # which will cause issues when matching
# We match SILVA taxonomy against this file (against scientific names) first

# Edit names_ncbi.dmp file into a file only containing non-scientific names
sed "s/ <[a-zA-Z -,.&:'0-9]*>//g" names_ncbi.dmp | grep -v 'scientific name' \
| cut -f 1,3 | awk -F $'\t' ' { t = $1; $1 = $2; $2 = t; print; } ' OFS=$'\t' \
| grep -v 'environmental\|uncultured\|unidentified\|metagenome' \
> NCBI_staxids_non_scientific.txt
	# 1. Remove <genus>, <family>, and other strings in <> brackets from taxonomy
  #    (otherwise the NCBI taxonomy won't match the SILVA taxonomy)
	# 2. Extract only non-scientific names and staxids
	# 3. Cut out columns we need
	# 4. Invert columns for script to work
	# 5. Removes lines containing "environmental", "uncultured", "unidentified",
  #    and "metagenome"
		# Needed because
			# SILVA taxonomy can have the same taxonomic ranks (e.g., "environmental
      # sample") for different higher ranks (e.g., "nematode; environmental
      # sample" and "bacteria;environmental sample"), wich would, however, be
      # assigned to the same staxid because the lower rank "environmental sample"
      # is similar
			# NCBI taxonomy can different staxids for the same taxonomic name, which
      # will cause issue when matching
# If SILVA taxonomy is not in exact scientific names, then we match against these
# to check for synonyms etc.


# To make a kraken2 DB, we need files with taxonomical info in the fashion of
# NCBI (like names.dmp, nodes.dmp, and seqid2taxid from NCBI). We're going to do
# that in python:

python3 $cwd/make_SILVA_databases.py

# Switched back to bash


# Merge accession numbers of fasta file and taxid file and add kraken2-specific
# info to each sequence name (when adding fasta sequences to a kraken2 DB, they
# need to contain the taxid in the name):
$cwd/subscripts/mergeFilesOnColumn.pl SILVA_${VER}_SSU_LSURef_NR99_tax_silva_trunc_edited.tab \
SILVA_accession_numbers_and_NCBI_taxids.txt 1 1 \
| sed 's/^/>/g' | awk -F "\t" '{ print $1 "|kraken:taxid|" $5 " " $2 "\n" $3}' \
> SILVA_${VER}_SSU_LSURef_NR99_tax_silva_trunc_kraken2.fasta


# Make the kraken2 DB
mkdir -p $kraken2DB
mkdir -p $kraken2DB/taxonomy
cp ${cwd}/make_SILVA_databases/nodes.dmp ${cwd}/make_SILVA_databases/names.dmp $kraken2DB/taxonomy # Prepare the taxonomy for the kraken2 DB
kraken2-build --add-to-library ${cwd}/make_SILVA_databases/SILVA_${VER}_SSU_LSURef_NR99_tax_silva_trunc_kraken2.fasta --db $kraken2DB
kraken2-build --build --db $kraken2DB

# Make the blast DB
mkdir -p $blastDB
cd $blastDB
makeblastdb -dbtype 'nucl' -in ${cwd}/make_SILVA_databases/SILVA_${VER}_SSU_LSURef_NR99_tax_silva_trunc.fasta \
-parse_seqids -taxid_map ${cwd}/make_SILVA_databases/SILVA_accession_numbers_and_NCBI_taxids.txt \
-out $(basename $blastDB)
cd $cwd
rm -r make_SILVA_databases/

echo "Database setup done."
