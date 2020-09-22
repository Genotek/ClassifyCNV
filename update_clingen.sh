#!/usr/bin/env bash

# download the latest files from ClinGen
wget -N ftp://ftp.clinicalgenome.org/ClinGen_haploinsufficiency_gene_GRCh37.bed 2>/dev/null || \
curl -o ClinGen_haploinsufficiency_gene_GRCh37.bed \
ftp://ftp.clinicalgenome.org/ClinGen_haploinsufficiency_gene_GRCh37.bed || exit 1

wget -N ftp://ftp.clinicalgenome.org/ClinGen_haploinsufficiency_gene_GRCh38.bed 2>/dev/null || \
curl -o ClinGen_haploinsufficiency_gene_GRCh38.bed \
ftp://ftp.clinicalgenome.org/ClinGen_haploinsufficiency_gene_GRCh38.bed || exit 1

wget -N ftp://ftp.clinicalgenome.org/ClinGen_triplosensitivity_gene_GRCh37.bed 2>/dev/null || \
curl -o ClinGen_triplosensitivity_gene_GRCh37.bed \
ftp://ftp.clinicalgenome.org/ClinGen_triplosensitivity_gene_GRCh37.bed || exit 1

wget -N ftp://ftp.clinicalgenome.org/ClinGen_triplosensitivity_gene_GRCh38.bed 2>/dev/null || \
curl -o ClinGen_triplosensitivity_gene_GRCh38.bed \
ftp://ftp.clinicalgenome.org/ClinGen_triplosensitivity_gene_GRCh38.bed || exit 1

wget -N ftp://ftp.clinicalgenome.org/ClinGen_region_curation_list_GRCh37.tsv 2>/dev/null || \
curl -o ClinGen_region_curation_list_GRCh37.tsv \
ftp://ftp.clinicalgenome.org/ClinGen_region_curation_list_GRCh37.tsv || exit 1

wget -N ftp://ftp.clinicalgenome.org/ClinGen_region_curation_list_GRCh38.tsv 2>/dev/null || \
curl -o ClinGen_region_curation_list_GRCh38.tsv \
ftp://ftp.clinicalgenome.org/ClinGen_region_curation_list_GRCh38.tsv || exit 1

# parse the tsv files into bed
python3 parse_clingen_tsv.py --infile ClinGen_region_curation_list_GRCh37.tsv || exit 1
python3 parse_clingen_tsv.py --infile ClinGen_region_curation_list_GRCh38.tsv || exit 1

# move the resulting files to Resources, delete the tsv files
mv ClinGen_haploinsufficiency_gene_GRCh37.bed Resources/hg19/ClinGen_haploinsufficiency_gene.bed || exit 1
mv ClinGen_triplosensitivity_gene_GRCh37.bed Resources/hg19/ClinGen_triplosensitivity_gene.bed || exit 1
mv ClinGen_haploinsufficiency_gene_GRCh38.bed Resources/hg38/ClinGen_haploinsufficiency_gene.bed || exit 1
mv ClinGen_triplosensitivity_gene_GRCh38.bed Resources/hg38/ClinGen_triplosensitivity_gene.bed || exit 1

mv ClinGen_region_curation_list_GRCh37.HI.bed Resources/hg19/ClinGen_region_curation_list.HI.bed || exit 1
mv ClinGen_region_curation_list_GRCh37.TS.bed Resources/hg19/ClinGen_region_curation_list.TS.bed || exit 1
mv ClinGen_region_curation_list_GRCh38.HI.bed Resources/hg38/ClinGen_region_curation_list.HI.bed || exit 1
mv ClinGen_region_curation_list_GRCh38.TS.bed Resources/hg38/ClinGen_region_curation_list.TS.bed || exit 1

rm ClinGen_region_curation_list_GRCh37.tsv
rm ClinGen_region_curation_list_GRCh38.tsv

cd Resources/hg19/
bedtools intersect -a ClinGen_region_curation_list.TS.bed -b refGenes.parsed.SelectTranscript.bed -wo |
awk '($4 == 40)' | grep 'NM_' | cut -f1-3,9 > Benign_TS_region_genelist.regions.bed
cat ClinGen_triplosensitivity_gene.bed | awk '($5 == 40)' | cut -f1-4 > Benign_TS_region_genelist.genes.bed
cat Benign_TS_region_genelist.regions.bed Benign_TS_region_genelist.genes.bed > Benign_TS_region_genelist.bed

cd ../hg38/
bedtools intersect -a ClinGen_region_curation_list.TS.bed -b refGenes.parsed.SelectTranscript.bed -wo |
awk '($4 == 40)' | grep 'NM_' | cut -f1-3,9 > Benign_TS_region_genelist.regions.bed
cat ClinGen_triplosensitivity_gene.bed | awk '($5 == 40)' | cut -f1-4 > Benign_TS_region_genelist.genes.bed
cat Benign_TS_region_genelist.regions.bed Benign_TS_region_genelist.genes.bed > Benign_TS_region_genelist.bed

cd ../../

# check that all files are not empty
if [ -s Resources/hg19/ClinGen_region_curation_list.HI.bed ] && [ -s Resources/hg19/ClinGen_region_curation_list.TS.bed ] \
&& [ -s Resources/hg38/ClinGen_region_curation_list.HI.bed ] && [ -s Resources/hg38/ClinGen_region_curation_list.TS.bed ] \
&& [ -s Resources/hg19/ClinGen_haploinsufficiency_gene.bed ] && [ -s Resources/hg19/ClinGen_triplosensitivity_gene.bed ] \
&& [ -s Resources/hg38/ClinGen_haploinsufficiency_gene.bed ] && [ -s Resources/hg38/ClinGen_triplosensitivity_gene.bed ] \
&& [ -s Resources/hg19/Benign_TS_region_genelist.bed ] && [ -s Resources/hg38/Benign_TS_region_genelist.bed ]
then
    echo "The ClinGen files are updated."
else
    echo "There was an error when updating the ClinGen files."
fi