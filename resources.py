#!/usr/bin/env python3
# coding=utf-8

"""
This module contains a list of resources used by the pipeline
"""

import os
import sys
from datetime import datetime
from collections import OrderedDict
import random
import string

# Default results path
home_dir = os.path.dirname(os.path.realpath(sys.argv[0]))  # directory where the script is actually located
main_results_folder = 'ClassifyCNV_results'
run_folder_prefix = 'Result_'
random_string = ''.join(random.choices(string.ascii_lowercase + string.digits, k=8))
run_results_folder = run_folder_prefix + datetime.now().strftime("%d-%b-%Y-%H-%M-%S-") + random_string
default_results_folder = os.path.join(main_results_folder, run_results_folder)
intermediate_folder = 'Intermediate_files'

# Filename for the cleaned input
cleaned_bed = 'infile.cleaned.bed'
cleaned_bed_path = os.path.join(intermediate_folder, cleaned_bed)

# Resources and databases
main_resources_folder = 'Resources'
common_resources_folder = 'common'
refgenes_db = 'refGenes.parsed.SelectTranscript.bed'
promoters_db = 'promoters.500bp.bed'
enhancers_db = 'Enhancers.3sources.merged.bed'
clingen_hi_db = 'ClinGen_haploinsufficiency_gene.bed'
clingen_ts_db = 'ClinGen_triplosensitivity_gene.bed'
clingen_regions_hi_db = 'ClinGen_region_curation_list.HI.bed'
clingen_regions_ts_db = 'ClinGen_region_curation_list.TS.bed'
gene_features_db = 'gene_features.bed'
penultimate_exon_50bp_db = '50bp_penultimate_exon.bed'
pop_freqs_db = 'population_freqs.bed'
decipher_HI_db = 'DECIPHER_HI_Predictions_Version3.bed'
pLI_db = 'ExAC_pLI.txt'
loeuf_db = 'gnomad.v2.1.1.lof_metrics.by_gene.txt'
decipher_HI_path = os.path.join(home_dir, main_resources_folder, common_resources_folder, decipher_HI_db)
pLI_path = os.path.join(home_dir, main_resources_folder, common_resources_folder, pLI_db)
loeuf_path = os.path.join(home_dir, main_resources_folder, common_resources_folder, loeuf_db)
benign_region_genes_db = 'Benign_TS_region_genelist.bed'

# Output files created by bedtools intersect
refgenes_intersect = 'refgenes_intersect.bed'
promoters_intersect = 'promoters_intersect.bed'
enhancers_intersect = 'enhancers_intersect.bed'
clingen_hi_intersect = 'clingen_hi_intersect.bed'
clingen_ts_intersect = 'clingen_ts_intersect.bed'
clingen_regions_hi_intersect = 'clingen_regions_hi_intersect.bed'
clingen_regions_ts_intersect = 'clingen_regions_ts_intersect.bed'
gene_features_intersect = 'gene_features_intersect.bed'
pop_freqs_intersect = 'population_freqs_intersect.bed'

refgenes_intersect_path = os.path.join(intermediate_folder, refgenes_intersect)
promoters_intersect_path = os.path.join(intermediate_folder, promoters_intersect)
enhancers_intersect_path = os.path.join(intermediate_folder, enhancers_intersect)
clingen_hi_intersect_path = os.path.join(intermediate_folder, clingen_hi_intersect)
clingen_ts_intersect_path = os.path.join(intermediate_folder, clingen_ts_intersect)
clingen_regions_hi_intersect_path = os.path.join(intermediate_folder, clingen_regions_hi_intersect)
clingen_regions_ts_intersect_path = os.path.join(intermediate_folder, clingen_regions_ts_intersect)
gene_features_intersect_path = os.path.join(intermediate_folder, gene_features_intersect)
pop_freqs_intersect_path = os.path.join(intermediate_folder, pop_freqs_intersect)

databases = {
    'genes': {'source': refgenes_db, 'result_path': refgenes_intersect_path},
    'promoters': {'source': promoters_db, 'result_path': promoters_intersect_path},
    'enhancers': {'source': enhancers_db, 'result_path': enhancers_intersect_path},
    'ClinGen_HI': {'source': clingen_hi_db, 'result_path': clingen_hi_intersect_path},
    'ClinGen_TS': {'source': clingen_ts_db, 'result_path': clingen_ts_intersect_path},
    'ClinGen_regions_HI': {'source': clingen_regions_hi_db, 'result_path': clingen_regions_hi_intersect_path},
    'ClinGen_regions_TS': {'source': clingen_regions_ts_db, 'result_path': clingen_regions_ts_intersect_path},
    'gene_features': {'source': gene_features_db, 'result_path': gene_features_intersect_path},
    'pop_freqs': {'source': pop_freqs_db, 'result_path': pop_freqs_intersect_path}
}

rubric = OrderedDict([
    ('1A-B', 0.0), ('2A', 0.0), ('2B', 0.0), ('2C', 0.0), ('2D', 0.0), ('2E', 0.0), ('2F', 0.0), ('2G', 0.0), ('2H', 0.0),
    ('2I', 0.0), ('2J', 0.0), ('2K', 0.0), ('2L', 0.0), ('3', 0.0), ('4A', 0.0), ('4B', 0.0), ('4C', 0.0), ('4D', 0.0),
    ('4E', 0.0), ('4F-H', 0.0), ('4I', 0.0), ('4J', 0.0), ('4K', 0.0), ('4L', 0.0), ('4M', 0.0), ('4N', 0.0),
    ('4O', 0.0), ('5A', 0.0), ('5B', 0.0), ('5C', 0.0), ('5D', 0.0), ('5E', 0.0), ('5F', 0.0), ('5G', 0.0), ('5H', 0.0)
])

# Printed results
scoresheet_filename = 'Scoresheet.txt'
scoresheet_header = '\t'.join(['VariantID', 'Chromosome', 'Start', 'End', 'Type', 'Classification', 'Total score']) + '\t'
scoresheet_header += '\t'.join(rubric.keys()) + '\t' + 'Known or predicted dosage-sensitive genes' + \
                     '\t' + 'All protein coding genes'

# Cutoffs for pathogenicity
pathogenic = 0.99
likely_pathogenic = 0.9
likely_benign = -0.9
benign = -0.99
