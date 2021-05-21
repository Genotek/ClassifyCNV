#!/usr/bin/env python3
# coding=utf-8

import argparse
from resources import *
from subprocess import Popen
import copy
from multiprocessing import Pool
import time


__version__ = '1.1.1'

parser = argparse.ArgumentParser()

parser.add_argument('--infile', required=True,
                    help='Input file in BED format; the first four columns should be chromosome, start position, '
                         'end position, CNV type (DEL or DUP).')
parser.add_argument('--GenomeBuild', required=True, choices=['hg19', 'hg38'],
                    help='Human assembly version (hg19 or hg38).')
parser.add_argument('--cores', type=int, default=1, help='Maximum number of threads to use. Default: 1')
parser.add_argument('--precise', action='store_true',
                    help='Specify this flag if the CNV breakpoints are precise. WARNING: if the breakpoints are not '
                         'precise, specifying the flag could lead to incorrect results. Default = False')
parser.add_argument('--outdir', default=default_results_folder,
                    help='Specify path to the run output directory. If no output directory is provided, results will '
                         'be saved to ClassifyCNV_results/Result_dd_Mon_yyyy-hh-mm-ss-{random}')
args = parser.parse_args()


def make_results_folder():
    """Creates a directory where all results and technical files will be saved to if it doesn't already exist.
    """
    # check if the main results folder already exists and if it is empty
    # if not, create it before proceeding
    if os.path.isdir(args.outdir):
        assert not os.listdir(args.outdir), "Results directory is not empty"
    else:
        if args.outdir == default_results_folder:
            os.makedirs(args.outdir)
        else:
            os.mkdir(args.outdir)
    # change to the results directory
    os.chdir(args.outdir)
    # make the intermediate folder where the technical files will be saved to
    os.mkdir(intermediate_folder)


def run_in_parallel(function, params_list, cores):
    """Runs a function in parallel.

    Args:
        function: A function to run.
        params_list: A list of arguments for the function.
        cores: Number of threads, taken from the command line arguments.

    Returns:
        returncodes: A list of return codes.

    """
    tasks = min([cores, len(params_list)])
    pool = Pool(tasks)
    results = [None] * len(params_list)
    for i, params in enumerate(params_list):
        results[i] = pool.apply_async(function, params)
    returncodes = [result.get() for result in results]
    pool.close()
    pool.join()
    return returncodes


def parse_infile(infile):
    """Parses the BED infile.
    Adds "chr" to the chromosome number if needed, removes duplicate entries.
    Creates a new "clean" BED file for further manipulations.

    Args:
        infile: Original CNV file submitted by the user.

    Returns:
        parsed_list: A list of CNVs in the chr_start_end_type format, for example, chr1_1000_2500_DEL.

    """
    parsed_list = set()
    bed_infile = open(infile, 'r')
    parsed_outfile = open(cleaned_bed_path, 'w')
    for line in bed_infile:
        fields = line.strip().split()
        if len(fields) < 4:
            sys.exit("ERROR: the input file must have at least 4 columns")
        # check that the line starts with the chromosome number, skip the line if it does not
        if fields[0][:3] not in ['chr', 'X', 'Y', 'M', '1', '2', '3', '4', '5', '6', '7', '8', '9',
                                 '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20',
                                 '21', '22']:
            continue
        # don't keep alternative contigs
        if fields[0].endswith('alt'):
            continue
        # stop execution if the fourth column does not contain the CNV type
        if fields[3] not in ['DEL', 'DUP']:
            sys.exit("ERROR: the 4th column of the input file does not contain the CNV type (DEL/DUP).")
        # each chromosome number should start with 'chr'
        if not fields[0].startswith('chr'):
            fields[0] = 'chr' + fields[0]
        # add the CNV id to the set
        cnv_to_add = "_".join([fields[0], fields[1], fields[2], fields[3]])
        if cnv_to_add not in parsed_list:
            parsed_list.add(cnv_to_add)
            parsed_outfile.write('\t'.join([fields[0], fields[1], fields[2], fields[3]]) + '\n')
    bed_infile.close()
    parsed_outfile.close()
    return parsed_list


def initialize_cnv_genes(cnv_list):
    """Creates a dictionary with an empty list of genes for each CNV.

    Args:
        cnv_list: A list of all CNVs to be analyzed with each CNV in the chr_start_end_type format.
    Returns:
        cnv_genes: A dictionary with an empty list for each CNV (to be populated later).
    """
    cnv_genes = dict()
    for cnv in cnv_list:
        cnv_genes.setdefault(cnv, [])
    return cnv_genes


def run_bedtools_intersect(file_b_type):
    """Runs the BEDTools intersect command.
    Intersects the cleaned BED file that contains CNVs with a database specified in the file_b_type variable.

    Args:
        file_b_type: Contents of the database.

    """
    if file_b_type == 'genes':
        extra_params = '| cut -f1-4,6-9 | sort -u'  # save the CNV and hit coordinates, type, transcript, gene name
    else:
        extra_params = ''

    db_path = os.path.join(home_dir, main_resources_folder, args.GenomeBuild, databases[file_b_type]['source'])
    intersect_command = ' '.join(['bedtools intersect -a', cleaned_bed_path, '-b', db_path, '-wo', extra_params, '>',
                                  databases[file_b_type]['result_path']])
    intersect_proc = Popen(['/bin/bash', '-c', intersect_command])
    intersect_proc.wait()
    if not intersect_proc.returncode == 0:
        print("ERROR when running BEDTools on ", file_b_type)
        sys.exit(1)


def genes_promoters_enhancers_intersect():
    """Uses BEDTools intersect to intersect the original CNV list with 3 databases: genes, promoters, enhancers.
    Assigns points for section 1 of the rubric (based on whether or not the CNV intersects anything at all).
    Assigns points for section 3 of the rubric (points depend on how many genes the CNV intersects).
    The points are saved to the global detailed_results dictionary.

    """
    element_counter = dict()  # key = CNV id, value = sum of the number of genes, promoters, enhancers the CNV overlaps
    for db_type in ['genes', 'promoters', 'enhancers']:
        run_bedtools_intersect(db_type)
        # parse the results of BEDTools intersect
        with open(databases[db_type]['result_path'], 'r') as res:
            for line in res:
                fields = line.strip().split()
                cnv_id = "_".join([fields[0], fields[1], fields[2], fields[3]])
                if cnv_id not in element_counter:
                    element_counter[cnv_id] = 1
                else:
                    element_counter[cnv_id] +=1
                # store overlapping protein-coding gene names
                if db_type == 'genes' and fields[6].startswith("NM_"):
                    if cnv_id not in cnv_genes:
                        cnv_genes[cnv_id] = [fields[7]]  # store the gene name
                    else:
                        cnv_genes[cnv_id].append(fields[7])  # store the gene name
                    # if the precise flag is on and the variant is a duplication, check if both breakpoints are
                    # within the same protein-coding gene
                    if args.precise:
                        # for deletions record intragenic variants
                        if fields[3] == 'DEL' and int(fields[4]) <= int(fields[1]) <= int(fields[2]) <= int(fields[5]):
                            breakpoints[cnv_id] = 1  # record the CNV as disrupting a protein coding gene
                            # make a record if a full gene is deleted exactly
                            if int(fields[4]) == int(fields[1]) and int(fields[2]) == int(fields[5]):
                                breakpoints[cnv_id] = 2
                        # for duplications record variants with one or two breakpoints inside a gene
                        elif fields[3] == 'DUP':
                            if int(fields[4]) <= int(fields[1]) <= int(fields[5]) or \
                                    int(fields[4]) <= int(fields[2]) <= int(fields[5]):
                                breakpoints[cnv_id] = 1  # record the CNV as disrupting a protein coding gene

        res.close()
        # assign points for section 3 of the rubric
        if db_type == 'genes':
            for cnv in cnv_genes:
                # assign points for section 3 (number of genes within the CNV) and save to final detailed results dict
                detailed_results[cnv]['3'] = assign_section3_points(len(cnv_genes[cnv]), cnv[-3:])

    # after we assembled the results from all 3 databases, assign points for section 1 (did anything overlap each CNV?)
    for cnv in cnv_list:
        if cnv not in element_counter:
            detailed_results[cnv]['1A-B'] = -0.6


def assign_section3_points(gene_number, cnv_type):
    """Assign points for section 3 of the ACMG rubric based on the number of genes included in the CNV.

    Args:
        gene_number: Number of genes included in the CNV.
        cnv_type: DEL or DUP.

    Returns:
        number_points: Number of points to assign.

    """
    number_points = 0.0
    if cnv_type.upper() == 'DEL':
        if gene_number < 25:
            number_points = 0.0
        elif gene_number < 35:
            number_points = 0.45
        elif gene_number > 34:
            number_points = 0.9
    elif cnv_type.upper() == 'DUP':
        if gene_number < 35:
            number_points = 0.0
        elif gene_number < 50:
            number_points = 0.45
        elif gene_number > 49:
            number_points = 0.9
    return number_points


def load_sensitive_genes(filename):
    """Saves the names of dosage sensitive genes.

    Args:
        filename: The name of the ClinGen genes file.

    Returns:
        gene_set: A set of gene names.

    """
    gene_set = set()
    file_path = os.path.join(home_dir, main_resources_folder, args.GenomeBuild, filename)
    clingen_in = open(file_path, 'r')
    for line in clingen_in:
        if line.startswith('chr'):
            fields = line.strip().split()
            if fields[4] == '3':
                gene_set.add(fields[3])
    clingen_in.close()
    return gene_set


def load_dosage_predictors():
    """Makes a list of gene names for which both pLI and DECIPHER scores suggest haploinsufficiency.
    This list is only used for evaluating deletions.

    Returns:
        predictors: A set of genes that are predicted to be dosage sensitive by both DECIPHER and ExAC.

    """
    predictors = set()  # only contains genes that are HI according to both DECIPHER and pLI
    decipher_set = set()  # a set of genes that are HI according to DECIPER
    pli_set = set()  # a set of genes that are HI according to pLI
    loeuf_set = set()  # a set of genes for which the upper bound of the o/e conf. interval indicates intolerance to LoF

    # load DECIPHER
    decipher_in = open(decipher_HI_path, 'r')
    for line in decipher_in:
        if line.startswith('chr'):
            fields = line.strip().split()
            subfields = fields[3].strip().split('|')  # parse out score info
            subfields[2] = float(subfields[2].replace('%', ''))
            if subfields[2] <= 10.0:  # check if the score is less than 10% (the gene is HI)
                decipher_set.add(subfields[0])  # add the gene name to the DECIPHER results
    decipher_in.close()

    # load pLI scores
    pli_in = open(pLI_path, 'r')
    for line in pli_in:
       if not line.startswith('transcript'):  # skip header
            fields = line.strip().split()
            if float(fields[19]) >= 0.9:  # check that the gene is intolerant
                pli_set.add(fields[1])  # add the gene name to the set
    pli_in.close()

    # load LOEUF
    loeuf_in = open(loeuf_path, 'r')
    for line in loeuf_in:
        if not line.startswith('gene'):  # skip header
            fields = line.strip().split()
            if not fields[29] == 'NA':  # don't evaluate empty scores
                if float(fields[29]) < 0.35:  # check that the gene is intolerant
                    loeuf_set.add(fields[0])  # add the gene name to the set
    loeuf_in.close()

    # only save genes that are present in all three datasets
    for cand_gene in decipher_set:
        if cand_gene in pli_set and cand_gene in loeuf_set:
            predictors.add(cand_gene)
    return predictors


def parse_established_regions(results_dict, cnv_type, file_regions, effect_column):
    """Parses the results of BEDTools intersect between the CNV list and the ClinGen databases.
    Checks for complete overlap with established dosage sensitive or benign regions/genes.

    Args:
        results_dict: A nested dictionary to save the results to.
        cnv_type: DEL or DUP.
        file_regions: The BEDTools intersect output to parse.
        effect_column: The number of the column that contains the dosage effect score.

    Returns:
        Modifies the results_dict dictionary.

    """
    regions_in = open(file_regions, 'r')
    for line in regions_in:
        fields = line.strip().split()
        cnv_id = '_'.join([fields[0], fields[1], fields[2], fields[3]])
        # we only want to look at one type of CNV at a time and only if it's an established dosage sensitive region
        if fields[3] == cnv_type and fields[effect_column] == '3':
            # check that the region is completely included inside the CNV and if so, save the result
            if int(fields[1]) <= int(fields[5]) and int(fields[2]) >= int(fields[6]):
                if cnv_id in results_dict:
                    results_dict[cnv_id]['full'] = 1
                else:
                    results_dict[cnv_id] = {}
                    results_dict[cnv_id]['full'] = 1
                # if it is a specific gene that overlapped, save its name to sensitive genes for this CNV
                if effect_column == 8:
                    if cnv_id in sensitive_genes:
                        sensitive_genes[cnv_id].append(fields[7])
                    else:
                        sensitive_genes[cnv_id] = [fields[7]]

        # check if the overlap is with an established benign region
        elif fields[3] == cnv_type and fields[effect_column] == '40':
            if cnv_type == 'DEL':
                # check if the CNV is completely contained inside the benign region and if so, save the result
                if int(fields[1]) >= int(fields[5]) and int(fields[2]) <= int(fields[6]):
                    if cnv_id in results_dict:
                        results_dict[cnv_id]['benign'] = 1
                    else:
                        results_dict[cnv_id] = {}
                        results_dict[cnv_id]['benign'] = 1
            elif cnv_type == 'DUP':
                # for duplications we need to keep track of whether the CNV is smaller or larger than the benign region
                # check if the CNV is identical to the benign region
                if int(fields[1]) == int(fields[5]) and int(fields[2]) == int(fields[6]):
                    if cnv_id in results_dict:
                        results_dict[cnv_id]['benign'] = 1
                    else:
                        results_dict[cnv_id] = {}
                        results_dict[cnv_id]['benign'] = 1
                # check if it is smaller than the benign region
                region_coords = '_'.join([fields[4], fields[5], fields[6]])  # save the coordinates of the benign region
                if int(fields[1]) >= int(fields[5]) and int(fields[2]) <= int(fields[6]):
                    if cnv_id in results_dict:
                        if 'benign_smaller' in results_dict[cnv_id]:
                            results_dict[cnv_id]['benign_smaller'].append(region_coords)
                        else:
                            results_dict[cnv_id]['benign_smaller'] = [region_coords]
                    else:
                        results_dict[cnv_id] = {}
                        results_dict[cnv_id]['benign_smaller'] = [region_coords]
                # check if it is larger than the benign region
                elif int(fields[1]) <= int(fields[5]) and int(fields[2]) >= int(fields[6]):
                    if cnv_id in results_dict:
                        if 'benign_larger' in results_dict[cnv_id]:
                            results_dict[cnv_id]['benign_larger'].append(region_coords)
                        else:
                            results_dict[cnv_id]['benign_larger'] = [region_coords]
                    else:
                        results_dict[cnv_id] = {}
                        results_dict[cnv_id]['benign_larger'] = [region_coords]

    regions_in.close()


def parse_gene_features(del_dict):
    """Analyzes which gene features are included in each CNV and saves them to the global detailed_results dictionary
    Ultimately, each gene that overlaps the CNV will have a list of features (5'UTR, 3'UTR, exons, cds) that are
    included in the CNV.

    Args:
        del_dict: Dosage results for deletions.

    Returns:
        Modifies the del_dict dictionary; adds gene names as inner keys and a list of features as values.

    """
    # load sensitive genes
    hi_genes = load_sensitive_genes(clingen_hi_db)
    file_in = open(gene_features_intersect_path, 'r')
    for line in file_in:
        fields = line.strip().split()
        if fields[8] in hi_genes:
            cnv_id = '_'.join([fields[0], fields[1], fields[2], fields[3]])
            # for deletions we need to know which features are affected in each gene
            if fields[3] == 'DEL':
                if cnv_id in del_dict:
                    if fields[8] in del_dict[cnv_id]:
                        del_dict[cnv_id][fields[8]].append(fields[9])  # inner key = gene name, value = feature
                    else:
                        del_dict[cnv_id][fields[8]] = [fields[9]]
                else:
                    del_dict[cnv_id] = {}
                    del_dict[cnv_id][fields[8]] = [fields[9]]
    file_in.close()


def load_benign_regions():
    """Loads lists of genes for each benign TS region.

    Returns:
        region_dictionary: A dictionary where key = benign region coordinates, value = a list of genes that belong
                           in the region.

    """
    region_dictionary = dict()
    file_path = os.path.join(home_dir, main_resources_folder, args.GenomeBuild, benign_region_genes_db)
    file_in = open(file_path, 'r')
    for line in file_in:
        fields = line.strip().split()
        # get the coordinates of the benign region
        id = '_'.join([fields[0], fields[1], fields[2]])
        # add the gene name to the list
        if id not in region_dictionary:
            region_dictionary[id] = [fields[3]]
        else:
            region_dictionary[id].append(fields[3])
    file_in.close()
    return region_dictionary


def assign_dup_points_s2(results):
    """For each duplication that has dosage information, assigns points for section 2 of the rubric.

    Args:
        results: A dictionary where key = CNV id, value = parsed dosage information.

    Returns:
        Modifies the global detailed_results dictionary.

    """
    # for each benign TS region get a list of genes that are expected within the benign region
    # anything that is not included here is outside of the benign region
    benign_region_genes = load_benign_regions()

    # iterate through each CNV
    for cnv in results:
        # check if the CNV fully overlapped with a dosage sensitive region
        if 'full' in results[cnv]:
            detailed_results[cnv]['2A'] = 1.0
        # check if the CNV overlapped completely with a known benign region
        elif 'benign' in results[cnv]:
            detailed_results[cnv]['2C'] = -1.0
        # check if the CNV is smaller than a known benign region and doesn't break genes (only if --precise flag is on)
        elif args.precise and 'benign_smaller' in results[cnv]:
            if cnv not in breakpoints:
                detailed_results[cnv]['2D'] = -1.0
        # check if the CNV overlapped with a known benign region but is larger than the region
        elif 'benign_larger' in results[cnv]:
            # get all gene names that belong inside the benign regions in this CNV
            benign_list = list()
            for matched_region in results[cnv]['benign_larger']:
                if matched_region in benign_region_genes:
                    benign_list.extend(benign_region_genes[matched_region])
            # go through the protein-coding genes that overlap the CNV and check if there are any that are not part
            # of the benign region
            tracker = 0
            for gene in cnv_genes[cnv]:
                if gene not in benign_list:
                    tracker = 1
            if tracker == 0:
                detailed_results[cnv]['2F'] = -1.0


def assign_del_points_s2(results):
    """For each deletion that has any dosage information, assigns points for section 2 of the rubric.

    Args:
        results: A dictionary where key = CNV id, value = parsed dosage information.

    Returns:
        Modifies the global detailed_results dictionary.

    """
    # load the names of HI genes
    for cnv in results:
        # check if a haploinsufficient region overlapped completely with the CNV
        if 'full' in results[cnv]:
            detailed_results[cnv]['2A'] = 1.0
        # check if the CNV is completely in a benign region
        elif 'benign' in results[cnv]:
            detailed_results[cnv]['2F'] = -1.0
        # if we already know this CNV is an intragenic variant and assigned points for it, skip the rest of evaluation
        elif detailed_results[cnv]['2E'] > 0:
            for gene_key in results[cnv].keys():
                if cnv in sensitive_genes:
                    sensitive_genes[cnv].append(gene_key)
                else:
                    sensitive_genes[cnv] = [gene_key]
            continue
        # if none of these conditions are true, the result can only contain individual gene information
        else:
            gene_keys = results[cnv].keys()
            for gene_key in gene_keys:
                results[cnv][gene_key] = set(results[cnv][gene_key])
                # check if the 5' part of the gene is deleted (5'UTR or first exon)
                if {'utr5', 'first_exon'} & results[cnv][gene_key]:
                    # check if coding sequence is involved
                    if 'cds' in results[cnv][gene_key]:
                        detailed_results[cnv]['2C'] = 0.9

                # check if the 3' part of the gene is deleted (3'UTR or the last exon)
                elif {'utr3', 'last_exon'} & results[cnv][gene_key]:
                    # check if other exons are involved
                    if {'exon', 'first_exon'} & results[cnv][gene_key]:
                        detailed_results[cnv]['2D'] = 0.9
                    # check if cds in the last exon is involved
                    elif 'cds' in results[cnv][gene_key]:
                        detailed_results[cnv]['2D'] = 0.3

                # save the gene to sensitive genes
                if cnv in sensitive_genes:
                    sensitive_genes[cnv].append(gene_key)
                else:
                    sensitive_genes[cnv] = [gene_key]


def analyze_intragenic_deletions(dosage_sensitive_cnv):
    """For intragenic deletions checks if the reading frame is disrupted and NMD is expected.
    The function is only executed if the --precise flag is on.

    Args:
        dosage_sensitive_cnv: A dictionary that contains dosage-sensitive deletions.

    Returns:
        intragenic_results: A set that contains variants that undergo NMD and where deletion causes a frameshift.

    """
    # load coordinates of the 3'-most 50 bp for the penultimate exon for each transcript
    # a gene is not expected to undergo NMD if only the last exon and/or this region are deleted
    nmd_coords_check = dict()  # contains the start and end of the 3' 50bp of second to last exon for each transcript
    nmd_path = os.path.join(home_dir, main_resources_folder, args.GenomeBuild, penultimate_exon_50bp_db)
    nmd_in = open(nmd_path, 'r')
    for region in nmd_in:
        parts = region.strip().split()
        nmd_coords_check[parts[3]] = {}
        nmd_coords_check[parts[3]]['start'] = int(parts[1])
        nmd_coords_check[parts[3]]['end'] = int(parts[2])
    nmd_in.close()

    intragenic_lengths = dict()  # technical dictionary to store deletion lengths
    nmd_results = dict()
    file_in = open(gene_features_intersect_path, 'r')  # bedtools intersect of CNVs with gene features
    for line in file_in:
        fields = line.strip().split()
        cnv_id = '_'.join([fields[0], fields[1], fields[2], fields[3]])
        # only look at deletions that we know have both breakpoints within the same gene,
        # are known to be dosage-sensitive and check only exon hits
        if cnv_id in breakpoints and fields[3] == 'DEL' and cnv_id in dosage_sensitive_cnv and 'exon' in fields[9]:
            # check for NMD; store 1 if CNV start or end is inside the last exon (evidence for NO NMD)
            if fields[9] == 'last_exon':
                # check if the start coordinate is inside the last exon
                if int(fields[5]) <= int(fields[1]) <= int(fields[6]):
                    if cnv_id in nmd_results:
                        nmd_results[cnv_id]['start'] = 1
                    else:
                        nmd_results[cnv_id] = {}
                        nmd_results[cnv_id]['start'] = 1
                # check if the end coordinate is inside the last exon
                if int(fields[5]) <= int(fields[2]) <= int(fields[6]):
                    if cnv_id in nmd_results:
                        nmd_results[cnv_id]['end'] = 1
                    else:
                        nmd_results[cnv_id] = {}
                        nmd_results[cnv_id]['end'] = 1
            else:
                # check if CNV start or end is inside the last 50 bp of the penultimate exon
                if fields[7] in nmd_coords_check:
                    if nmd_coords_check[fields[7]]['start'] <= int(fields[1]) <= nmd_coords_check[fields[7]]['end']:
                        if cnv_id in nmd_results:
                            nmd_results[cnv_id]['start'] = 1
                        else:
                            nmd_results[cnv_id] = {}
                            nmd_results[cnv_id]['start'] = 1
                    if nmd_coords_check[fields[7]]['start'] <= int(fields[2]) <= nmd_coords_check[fields[7]]['end']:
                        if cnv_id in nmd_results:
                            nmd_results[cnv_id]['end'] = 1
                        else:
                            nmd_results[cnv_id] = {}
                            nmd_results[cnv_id]['end'] = 1
            # count the number of deleted bases
            if cnv_id in intragenic_lengths:
                intragenic_lengths[cnv_id] += int(fields[10])
            else:
                intragenic_lengths[cnv_id] = int(fields[10])

    file_in.close()

    intragenic_results = set()
    # identify PVS1 variants - look for frameshift and variants that don't have two 1's in nmd_results
    for variant in intragenic_lengths:
        # check if the variant undergoes NMD, if not - skip
        if variant in nmd_results:
            if {'start', 'end'} <= nmd_results[variant].keys():
                continue
        # check if frameshift occurs, if not - skip
        if intragenic_lengths[variant] % 3 == 0:
            continue
        # now the variant is definitely causing a frameshift and NMD, add it to the results list
        intragenic_results.add(variant)
    return intragenic_results


def assign_points_intragenic_del_2e(pvs1_list):
    """Assigns points for section 2E of gene loss rubric.
    Args:
        pvs1_list: A list of deletions that are known to cause frameshift and NMD.

    Returns:
        Modifies the detailed_results dictionary by adding points for section 2E.

    """
    for variant in pvs1_list:
        detailed_results[variant]['2E'] = 0.9

    # check if there are any dosage-sensitive genes that were deleted entirely and assign points
    for variant in breakpoints:
        if breakpoints[variant] == 2:
            detailed_results[variant]['2E'] = 0.9


def dosage_sensitivity():
    """Runs BEDTools intersect for CNVs against dosage sensitivity databases, parses results, assigns points for
    section 2 of the rubric.

    """
    # run BEDTools intersect in parallel
    list_to_run = list()
    for db_type in ['ClinGen_HI', 'ClinGen_TS', 'ClinGen_regions_HI', 'ClinGen_regions_TS', 'gene_features']:
        list_to_run.append(db_type)
    dosage_params = [[a] for a in list_to_run]
    run_in_parallel(run_bedtools_intersect, dosage_params, args.cores)

    # initiate dosage sensitivity result dictionaries (step 2 of the rubric)
    dosage_res_dup = dict()
    dosage_res_del = dict()

    # parse complete overlaps with established HI and TS regions and genes
    parse_established_regions(dosage_res_del, 'DEL', clingen_regions_hi_intersect_path, 7)
    parse_established_regions(dosage_res_dup, 'DUP', clingen_regions_ts_intersect_path, 7)
    parse_established_regions(dosage_res_del, 'DEL', clingen_hi_intersect_path, 8)
    parse_established_regions(dosage_res_dup, 'DUP', clingen_ts_intersect_path, 8)

    # for protein coding genes, check which gene features overlap with each CNV
    parse_gene_features(dosage_res_del)

    # if the --precise flag is used, check the intragenic deletions and see if reading frames are disrupted
    if args.precise:
        intragenic_deletions = analyze_intragenic_deletions(dosage_res_del)  # creates a set of PVS1 variants
        assign_points_intragenic_del_2e(intragenic_deletions)

    # assign points for section 2 for duplications
    assign_dup_points_s2(dosage_res_dup)
    # assign points for section 2 for deletions
    assign_del_points_s2(dosage_res_del)

    # assign points for section 2H (predicted HI)
    assign_HI_predictor_points()


def assign_HI_predictor_points():
    """Assigns points for section 2H of the deletion rubric.
    For each CNV we check if any of the genes are predicted to be haploinsufficient.
    If DECIPHER and gnomAD both consider a gene HI, points are assigned.

    Returns:
        Modifies the global detailed_results dictionary.

    """
    # load a set of genes that are predicted to be haploinsufficient by both gnomAD and DECIPHER
    predicted_hi_genes = load_dosage_predictors()

    # iterate through CNVs and check which ones have any of the predicted HI genes
    for test_cnv in cnv_genes:
        if test_cnv.endswith('DEL'):  # we are only looking at deletions
            for test_gene in cnv_genes[test_cnv]:
                # don't evaluate if this is a known dosage sensitive gene
                if test_cnv in sensitive_genes:
                    if test_gene in sensitive_genes[test_cnv]:
                        continue
                # don't evaluate if we already assigned points for an internal or an established benign variant
                if detailed_results[test_cnv]['2E'] > 0 or detailed_results[test_cnv]['2F'] < 0:
                    continue
                if test_gene in predicted_hi_genes:  # check the gene is in the predicted HI list
                    detailed_results[test_cnv]['2H'] = 0.15
                    # save this gene to the list of dosage sensitive genes for this CNV
                    if test_cnv in sensitive_genes:
                        sensitive_genes[test_cnv].append(test_gene)
                    else:
                        sensitive_genes[test_cnv] = [test_gene]


def analyze_pop_freqs():
    """Runs BEDTools intersect against the population frequency database to identify common CNVs.
    Parses results, assigns points for section 4O. Only looks at overlaps >= 80% of the CNV.
    Doesn't analyze CNVs that contain known or predicted dosage sensitive genes.

    Returns:
        Modifies the global detailed_results dictionary by adding points for section 4O.

    """
    # run BEDTools intersect
    run_bedtools_intersect('pop_freqs')

    # parse BEDTools results
    pop_freqs_res = dict()  # save results for each CNV
    pop_freqs_in = open(pop_freqs_intersect_path, 'r')
    for line in pop_freqs_in:
        fields = line.strip().split()
        # check if the overlap is with the same type of CNV
        if fields[3] == fields[7]:
            cnv_id = '_'.join([fields[0], fields[1], fields[2], fields[3]])

            # don't analyze these data for CNVs that overlap known or predicted dosage sensitive regions
            if max(detailed_results[cnv_id]['2A'], detailed_results[cnv_id]['2C'], detailed_results[cnv_id]['2D'],
                   detailed_results[cnv_id]['2E'], detailed_results[cnv_id]['2H']) > 0:
                continue

            overlap_perc = int(fields[9]) * 100 / (int(fields[2]) - int(fields[1]))
            # don't keep very short overlaps
            if overlap_perc < 80:
                continue

            # add the variant frequency to the results for the CNV
            if cnv_id in pop_freqs_res:
                pop_freqs_res[cnv_id].append(float(fields[8]))
            else:
                pop_freqs_res[cnv_id] = [float(fields[8])]
    pop_freqs_in.close()

    # iterate through results and assign points
    for cnv_id in pop_freqs_res:
        # if there was only one pop frequency variant match, use its frequency
        if len(pop_freqs_res[cnv_id]) == 1:
            if pop_freqs_res[cnv_id][0] > 1.0:
                detailed_results[cnv_id]['4O'] = -1.0
        # if there are more than 1 pop frequency variant overlaps, take the average of frequencies of all overlaps
        else:
            average = sum(pop_freqs_res[cnv_id])/len(pop_freqs_res[cnv_id])
            if average > 1.0:
                detailed_results[cnv_id]['4O'] = -1.0


def generate_results():
    """Uses the detailed_results dictionary to calculate the final score, converts it to pathogenicity status
    and prints full results to file.

    """
    results_out = open(scoresheet_filename, 'w')
    results_out.write(scoresheet_header + '\n')
    for cnv in sorted(cnv_list):
        # add up individual scores for each element in the rubric to get the final score
        final_score = sum(detailed_results[cnv].values())
        # convert final score to pathogenicity
        if final_score >= pathogenic:
            status = 'Pathogenic'
        elif final_score >= likely_pathogenic:
            status = 'Likely pathogenic'
        elif final_score > likely_benign:
            status = 'Uncertain significance'
        elif final_score > benign:
            status = 'Likely benign'
        else:
            status = 'Benign'

        # make a list of dosage sensitive genes contained within the CNV
        if cnv in sensitive_genes:
            genes_to_print = set(sensitive_genes[cnv])
            genes_to_print_str = ', '.join(genes_to_print)
        else:
            genes_to_print_str = ''
        # make a list of all protein coding genes contained within the CNV
        if cnv in cnv_genes:
            all_genes = ', '.join(cnv_genes[cnv])
        else:
            all_genes = ''
        # convert each value in the results dictionary to a string
        for k in detailed_results[cnv]:
            detailed_results[cnv][k] = str(detailed_results[cnv][k])

        # print results to file
        cnv_info = cnv.strip().split('_')
        formatted_result = '\t'.join([cnv, cnv_info[0], cnv_info[1], cnv_info[2], cnv_info[3], status,
                                      str(format(final_score, '.2f'))]) + '\t'
        formatted_result += '\t'.join(detailed_results[cnv].values()) + '\t'
        formatted_result += genes_to_print_str + '\t' + all_genes
        results_out.write(formatted_result + '\n')
    results_out.close()


if __name__ == "__main__":
    t_start = time.perf_counter()  # time the run
    print(__file__, 'Version', __version__)

    # initialize results dictionaries
    detailed_results = dict()  # contains a breakdown of the final pathogenicity score
    sensitive_genes = dict()  # contains a list of dosage sensitive genes for each CNV
    if args.precise:
        breakpoints = dict()  # stores intragenic CNVs

    infile_path = os.path.abspath(args.infile)
    make_results_folder()  # create a folder where the results will be stored if it doesn't already exist
    cnv_list = parse_infile(infile_path)  # save each CNV as chr_start_end_type and print a new file for BEDTools

    # make empty result dictionaries
    for cnv in cnv_list:
        detailed_results[cnv] = copy.deepcopy(rubric)
    # initialize the cnv_genes dictionary which will contain a list of protein-coding genes that are included in
    # each CNV
    cnv_genes = initialize_cnv_genes(cnv_list)
    # intersect CNVs with genes, promoters, enhancers and assign points for steps 1,3
    genes_promoters_enhancers_intersect()
    # check if CNVs are in dosage sensitive regions and assign points for step 2
    dosage_sensitivity()
    # check if any of the CNVs are frequent and assign points for section 4O
    analyze_pop_freqs()
    # calculate the total score, determine pathogenicity, print results to file
    generate_results()

    t_stop = time.perf_counter()
    t_fact = t_stop - t_start
    print('Results saved to', args.outdir)
    print('Elapsed time:', '{0:.2f}'.format(t_fact), 'seconds')
