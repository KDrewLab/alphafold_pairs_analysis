'''
Code for extracting score metrics from AlphaFold3 runs
Written by Miles Woodcock-Girard for Drew Lab at UIC
'''

import json
import csv
import os
import sys
import glob
import subprocess
import argparse

import pdockq as pDockQ
import pdockq2 as pDockQ2
import uniprot_client
#import pae
import utils

import pandas as pd
import numpy as np

dc2_csv = "dc2_af3.csv"
dc2_df = pd.read_csv(dc2_csv)


parser = argparse.ArgumentParser(description="Obtain pairwise scoring metrics for AlphaFold output")
#parser.add_argument("--af_model", action="store", dest="model_file", required=True,
#                    help="Directory of AlphaFold output. Can be from server or local run.")

parser.add_argument("--variants", action="store_true", help="Whether to query UniProt for disease variants")

parser.add_argument("--screen_interactions", action="store_true", help="Whether to query IntAct for interaction support")

parser.add_argument("--scores_outfmt", type=str, default="gene_id1 gene_id2 iptm dc2 pdockq pdockq2 ipsae",
                    help="Space-separated list of output score fields to include in output\n  e.g., 'gene_id1 gene_id2 iptm dc2 pdockq pdockq2 ipsae lis'")

parser.add_argument("--contact_threshold", action="store", dest="contact_threshold", default=8,
                    help="Maximum interchain CA-CA distance (Å) to be considered 'in contact'")

parser.add_argument("--pae_threshold", action="store", dest="pae_threshold", default=12,
                    help="Maximum predicted aligned error (PAE) value to consider when predicting ipSAE, pDockQ2, LIS")

parser.add_argument("--output_filename", action="store", dest='outfile', required=True,
                    help="Name for output metrics file")

args = parser.parse_args()

def parse_args_for_outfields(args):
    # Initialize base 'out_fields' list for determining output format
    out_fields = ['interaction_id']

    # Parse user-defined '--outfmt' flag, adding desired information/scores
    if args.scores_outfmt:
        out_fields += args.scores_outfmt.strip().split()
    else:
        out_fields += ['iptm', 'dc2', 'pdockq', 'pdockq2', 'ipsae', 'lis']

    # Parse '--variants' flag for querying UniProt for disease variants with publications
    if args.variants:
        out_fields += ['has_pub1', 'has_pub2']

    # Parse '--screen_interaction' flag for querying IntAct for experimental support of interaction
    if args.screen_interactions:
        out_fields += ['exp_support']

    #print(out_fields)
    #exit()

    return out_fields

        

out_fields = parse_args_for_outfields(args)


# Function for retrieving DC2 score, given a pair of proteins
def get_dc2_score_from_pair(prot1, prot2):
    # Get DirectContacts2 score
    dc2_row = dc2_df[dc2_df["fset"] == f"frozenset({{\'{prot1.upper()}\', \'{prot2.upper()}\'}})"]
    if dc2_row.empty:
        dc2_row = dc2_df[dc2_df["fset"] == f"frozenset({{\'{prot2.upper()}\', \'{prot1.upper()}\'}})"]

    dc2_score_str = dc2_row.iloc[0]["score"]
    dc2_score = float(dc2_score_str)

    return dc2_score



# Function for calculating pDockQ, given a pair mmCIF file
def get_pdockq_from_mmcif(mmcif_path, contact_threshold=8):
    # Obtain chain coordinates, pLDDTs from mmCIF file
    chain_coords, chain_plddt = pDockQ.read_model_file(mmcif_path)

    # Calculate pDockQ score for protein pair
    pdockq, ppv = pDockQ.calc_pdockq(chain_coords, chain_plddt, contact_threshold)

    return pdockq



# Function for calculating pDockQ2, given a pair mmCIF file, JSON, threshold
def get_pdockq2_from_mmcif(mmcif_path, json_path, contact_threshold=8):
    pdockq2_df = pDockQ2.main(mmcif_path, json_path, contact_threshold)

    pmidockqs = pdockq2_df[['pmidockq']].values

    pmidockqs_min = np.min(pmidockqs)

    return pmidockqs_min



# Function for reading ipSAE .txt, returning contained scores (ipSAE, pDockQ, pDockQ2)
def get_scores_from_ipsae_txt(ipsae_txt):
    try:
        ipsae_txt_file = open(ipsae_txt)
        ipsae_txt_data = ipsae_txt_file.readlines()
    except FileNotFoundError:
        print(f"Error: file '{ipsae_txt}' not found. Skipping ...")
        return

    ipsae_data = ipsae_txt_data[4].split()

    ipsae = ipsae_data[5]
    iptm = ipsae_data[8]
    pdockq = ipsae_data[10]
    pdockq2 = ipsae_data[11]
    lis = ipsae_data[12]

    return ipsae, iptm, pdockq, pdockq2, lis



# Function for calculating ipSAE, given a pair mmCIF file, JSON, thresholds
def get_scores_from_mmcif(mmcif_path, json_path, pae_threshold, contact_threshold):
    # Calculate ipSAE score for protein pair
    if len(str(pae_threshold)) == 1:
        pae_threshold = str(0) + str(pae_threshold)

    if len(str(contact_threshold)) == 1:
        contact_threshold = str(0) + str(contact_threshold)

    ipsae_output = mmcif_path[:-4:] + f"_{pae_threshold}_{contact_threshold}.txt"
    subprocess.run(["python3", "ipsae.py", json_path, mmcif_path, str(pae_threshold), str(contact_threshold)])
    #with open(ipsae_output, 'r') as ipsae_file:
    #    ipsae_data = ipsae_file.read().splitlines()[1:-1:][3].split()
    #    ipsae_max = float(ipsae_data[5])
    ipsae, iptm, pdockq, pdockq2, lis = get_scores_from_ipsae_txt(ipsae_output)

    return ipsae, iptm, pdockq, pdockq2, lis



# Function for checking if user-defined output filename exists, extracting already processed pairs into set
def get_already_processed_pairs(output_file):
    existing_pairs = set()

    csv_exists = os.path.exists(args.outfile)
    if csv_exists:
        with open(args.outfile, "r") as f:
            reader = csv.reader(f)
            next(reader)  # Skip header
            for row in reader:
                if row:   # Avoid blank lines
                    existing_pairs.add(row[0])

    return csv_exists, existing_pairs


# Function for getting user-specified information/scores from user-specified pair model file
def get_metrics_for_pair_file(af3_pair_output, pae_threshold, contact_threshold):

    if utils.usedAF3Server(af3_pair_output):

        if not glob.glob(f"./{af3_pair_output}/fold_{af3_pair_output}_*_summary_confidences_*.json"):
            print(f"No 'summary_confidences_#.json' found. Skipping ...")
            return 11 * [-1]

        best_model = utils.getBestModel(af3_pair_output)

        if not glob.glob(f"./{af3_pair_output}/fold_{af3_pair_output}_*_summary_confidences_{best_model}.json"):
            print(f"No 'summary_confidences_{best_model}.json' found. Skipping ...")
            return 11 * [-1]

        json_scores = glob.glob(f"./{af3_pair_output}/fold_{af3_pair_output}_*_full_data_{best_model}.json")[0]
        json_summary = glob.glob(f"./{af3_pair_output}/fold_{af3_pair_output}_*_summary_confidences_{best_model}.json")[0]
        mmcif_path = glob.glob(f"./{af3_pair_output}/fold_{af3_pair_output}_*_model_{best_model}.cif")[0]

    else:
        json_scores = f"./{af3_pair_output}/{af3_pair_output}_confidences.json"
        json_summary = f"./{af3_pair_output}/{af3_pair_output}_summary_confidences.json"
        mmcif_path = f"./{af3_pair_output}/{af3_pair_output}_model.cif"

    if not os.path.exists(json_summary):
        print(f"No 'summary_confidences' JSON found. Skipping ...")
        return 11 * [-1]

    if not os.path.exists(json_scores):
        print(f"No 'confidences' JSON found. Skipping ...")
        return 11 * [-1]

    if not os.path.exists(mmcif_path):
        print(f"No 'model' CIF found. Skipping ...")
        return 11 * [-1]

    # Get UniProt IDs of proteins in pairwise model. Assumes AlphaFold output directory name takes form: 'UniProtID1_UniProtID2_...'
    uniprot_id1 = af3_pair_output.split("_")[0]
    uniprot_id2 = af3_pair_output.split("_")[1]
    interaction_id = "_".join([uniprot_id1, uniprot_id2])

    # Get DirectContacts2 score for pair of uniprot IDs
    if 'dc2' in out_fields:
        dc2_score = get_dc2_score_from_pair(uniprot_id1, uniprot_id2)
    else:
        dc2_score = -1

    # Get pairwise scoring metrics using the ipSAE script
    ipsae, iptm, pdockq, pdockq2, lis = get_scores_from_mmcif(mmcif_path, json_scores, pae_threshold, contact_threshold)

    # Query UniProt for the corresponding gene names
    if 'gene_id1' in out_fields:
        gene_id1 = uniprot_client.get_gene_id(uniprot_id1)
    else:
        gene_id1 = ""

    if 'gene_id2' in out_fields:
        gene_id2 = uniprot_client.get_gene_id(uniprot_id2)
    else:
        gene_id2 = ""

    # Query UniProt for publications associated with pathogenic variants of genes
    if 'has_pub1' in out_fields and 'has_pub2' in out_fields:
        has_pub1 = uniprot_client.has_publication(uniprot_id1)
        has_pub2 = uniprot_client.has_publication(uniprot_id2)
    else:
        has_pub1 = -1
        has_pub2 = -1

    # Query IntAct for whether interaction has experimental support
    if 'exp_support' in out_fields:
        known_interact = uniprot_client.interaction_in_uniprot(uniprot_id1, uniprot_id2)
    else:
        known_interact = -1


    # Return all information for pair
    return interaction_id, gene_id1, gene_id2, iptm, dc2_score, pdockq, pdockq2, ipsae, has_pub1, has_pub2, known_interact






outfile_exists, processed_pairs = get_already_processed_pairs(args.outfile)
with open(args.outfile, 'a') as output_file:
    outfile_writer = csv.writer(output_file, delimiter = ',')

    if not outfile_exists:
        outfile_writer.writerow(out_fields)

    for i, af_pair_output in enumerate(os.listdir()):

        # Determine if pair has been processed
        if af_pair_output in processed_pairs:
            print(f"Skipping already processed pair: {af_pair_output}")
            continue

        if not os.path.isdir(af_pair_output):
            continue

        print(f"Now processing: {i} {af_pair_output}")

        interaction_id, gene_id1, gene_id2, iptm, dc2, pdockq, pdockq2, ipsae, has_pub1, has_pub2, known_interact = get_metrics_for_pair_file(af_pair_output,
                                                                                                                                              args.pae_threshold,
                                                                                                                                              args.contact_threshold)

        if interaction_id == -1:
            continue

        row_data = {
                "interaction_id": interaction_id,
                "gene_id1": gene_id1,
                "gene_id2": gene_id2,
                "iptm": iptm,
                "dc2": dc2,
                "pdockq": pdockq,
                "pdockq2": pdockq2,
                "ipsae": ipsae,
                "has_pub1": has_pub1,
                "has_pub2": has_pub2,
                "exp_support": known_interact
        }

        #print([row_data.get(field, "") for field in out_fields])

        outfile_writer.writerow([row_data.get(field, "") for field in out_fields])




'''
with open("metrics.csv", "a") as csv_file:
    csv_writer = csv.writer(csv_file, delimiter = ",")Å

    if not csv_exists:
        csv_writer.writerow(["Pair", "Gene1", "Gene2", "ipTM", "pTM", "Ranking Score", "DC2 Score", "pDockQ", "pDockQ2", "ipSAE", "Path. Var. Publication (A)", "Path. Var. Publication (B)", "Novel Interaction"])

    for i, pair in enumerate(os.listdir()):
        if os.path.isdir(pair) == False:
            continue

        # Determine if pair has been processed
        if pair in existing_pairs:
            print(f"Skipping already processed pair: {pair}")
            continue

        print(f"{i} {pair}")

        if utils.usedAF3Server(pair):

            if not glob.glob(f"./{pair}/fold_{pair}_*_summary_confidences_*.json"):
                print(f"Empty directory: {pair} Skipping ...")
                continue

            best_model = utils.getBestModel(pair)

            if not glob.glob(f"./{pair}/fold_{pair}_*_summary_confidences_{best_model}.json"):
                print(f"Empty directory: {pair} Skipping ...")
                continue

            curr_json_scores = glob.glob(f"./{pair}/fold_{pair}_*_full_data_{best_model}.json")[0]
            curr_json_summary = glob.glob(f"./{pair}/fold_{pair}_*_summary_confidences_{best_model}.json")[0]
            curr_mmcif_path = glob.glob(f"./{pair}/fold_{pair}_*_model_{best_model}.cif")[0]

        else:
            curr_json_scores = f"./{pair}/{pair}_confidences.json"
            curr_json_summary = f"./{pair}/{pair}_summary_confidences.json"
            curr_mmcif_path = f"./{pair}/{pair}_model.cif"

        if not os.path.exists(curr_json_summary) or not os.path.exists(curr_mmcif_path):
            continue


        prot1 = pair.split("_")[0]
        prot2 = pair.split("_")[1]

        with open(curr_json_summary) as curr_json_file:
            json_data = json.load(curr_json_file)
    
            # Obtain ipTM, pTM, Ranking Scores from confidences JSON
            iptm = json_data["iptm"]
            ptm = json_data["ptm"]
            ranking_score = json_data["ranking_score"]

            # Get DirectContacts2 score
            dc2_score = get_dc2_score_from_pair(prot1, prot2)
            
            # Calculate pDockQ score for protein pair
            #pdockq = get_pdockq_from_mmcif(curr_mmcif_path)

            # Calculate pDockQ2 score for protein pair
            #pdockq2 = get_pdockq2_from_mmcif(curr_mmcif_path, curr_json_scores)

            # Calculate ipSAE score for protein pair (also outputs pDockQ, pDockQ2, LIS)
            #ipsae = get_ipsae_from_mmcif(curr_mmcif_path, curr_json_scores, str(15), str(15))

            # Obtain scores from mmCIF file, data JSON
            ipsae, pdockq, pdockq2, lis = get_scores_from_mmcif(curr_mmcif_path, curr_json_scores,
                                                                pae_threshold, contact_threshold)

            # Get gene names from Uniprot
            geneID1 = uniprot_client.get_gene_id(prot1)
            geneID2 = uniprot_client.get_gene_id(prot2)

            # Discern whether disease variants exist for protein
            #numVars1 = uniprot_client.count_disease_var(prot1)
            #numVars2 = uniprot_client.count_disease_var(prot2)

            # Check if variants have associated publications
            has_pub1 = uniprot_client.has_publication(prot1)
            has_pub2 = uniprot_client.has_publication(prot2)

            # Determine if interaction can be found on UniProt
            known_interact = uniprot_client.interaction_in_uniprot(prot1, prot2)

            # Write pair data to row in CSV output
            csv_writer.writerow([pair, geneID1, geneID2, iptm, ptm, ranking_score, dc2_score, pdockq, pdockq2, ipsae, has_pub1, has_pub2, known_interact])

            # Update metrics in dc2_df
            fset_key_f = f"frozenset({{\'{prot1.upper()}', '{prot2.upper()}\'}})"
            fset_key_r = f"frozenset({{\'{prot2.upper()}', '{prot1.upper()}\'}})"

            dc2_index = dc2_df[dc2_df["fset"] == fset_key_f].index
            if dc2_index.empty:
                dc2_index = dc2_df[dc2_df["fset"] == fset_key_r].index

            if not dc2_index.empty:
                idx = dc2_index[0]
                dc2_df.at[idx, "AF3_ipTM"] = iptm
                dc2_df.at[idx, "AF3_pTM"] = ptm
                dc2_df.at[idx, "pDockQ"] = pdockq
            else:
                print(f"Warning: could not find fset for {pair}")
            
            curr_json_file.close()


    csv_file.close()
'''
# Dump final dataframe to CSV
dc2_df.to_csv(dc2_csv, index=False)
