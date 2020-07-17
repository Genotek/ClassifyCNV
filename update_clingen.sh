#!/usr/bin/env bash

# download the latest files from ClinGen
wget ftp://ftp.clinicalgenome.org/ClinGen_haploinsufficiency_gene_GRCh37.bed
wget ftp://ftp.clinicalgenome.org/ClinGen_haploinsufficiency_gene_GRCh38.bed
wget ftp://ftp.clinicalgenome.org/ClinGen_triplosensitivity_gene_GRCh37.bed
wget ftp://ftp.clinicalgenome.org/ClinGen_triplosensitivity_gene_GRCh38.bed
wget ftp://ftp.clinicalgenome.org/ClinGen_region_curation_list_GRCh37.tsv
wget ftp://ftp.clinicalgenome.org/ClinGen_region_curation_list_GRCh38.tsv

# parse the tsv files into bed
python3 parse_clingen_tsv.py --infile ClinGen_region_curation_list_GRCh37.tsv
python3 parse_clingen_tsv.py --infile ClinGen_region_curation_list_GRCh38.tsv

# move the resulting files to Resources, delete the tsv files
mv ClinGen_haploinsufficiency_gene_GRCh37.bed Resources/hg19/ClinGen_haploinsufficiency_gene.bed
mv ClinGen_triplosensitivity_gene_GRCh37.bed Resources/hg19/ClinGen_triplosensitivity_gene.bed
mv ClinGen_haploinsufficiency_gene_GRCh38.bed Resources/hg38/ClinGen_haploinsufficiency_gene.bed
mv ClinGen_triplosensitivity_gene_GRCh38.bed Resources/hg38/ClinGen_triplosensitivity_gene.bed

mv ClinGen_region_curation_list_GRCh37.HI.bed Resources/hg19/ClinGen_region_curation_list.HI.bed
mv ClinGen_region_curation_list_GRCh37.TS.bed Resources/hg19/ClinGen_region_curation_list.TS.bed
mv ClinGen_region_curation_list_GRCh38.HI.bed Resources/hg38/ClinGen_region_curation_list.HI.bed
mv ClinGen_region_curation_list_GRCh38.TS.bed Resources/hg38/ClinGen_region_curation_list.TS.bed

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

echo "The ClinGen files are updated"