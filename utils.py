# General utilities for processing AlphaFold3 output
# Written by Miles Woodcock-Girard for Drew Lab at UIC

import os
import glob
import json
import csv


# Function to convert mmCIF to PDB, using BeEM package
#def cif2pdb(mmcif_input, 


# Boolean function for determining whether AF3 model was generated using AF3 server
def usedAF3Server(afOutputDir):

    subdirs = [file for file in os.listdir(afOutputDir)
               if os.path.isdir(os.path.join(afOutputDir, file))]

    files = [file for file in os.listdir(afOutputDir)]

    # AlphaFold3 server output has no subdirectories
    if not subdirs:
        return True
    # AlphaFold3 local output has five subdirectories, one for each model (0 - 4)
    elif len(subdirs) == 5:
        return False
    # Sanity check -- this case should never occur
    else:
        print("ERROR: This should not happen")
        return False


# Function to get the best-scoring AF3 model number (0 - 4), per their "ranking_score"
# in their confidence summary JSON
def getBestModel(afOutputDir):

    usedServer = usedAF3Server(afOutputDir)

    # AF3 webserver outputs have simpler directory structure
    if usedServer:
        summaryFiles = glob.glob(os.path.join(afOutputDir, "*_summary_confidences_*.json"))

        if not summaryFiles:
            return -1

        ranking_scores = []
        for summaryFile in summaryFiles:
            with open(summaryFile, 'r') as summary_json:
                json_data = json.load(summary_json)
                ranking_scores.append(float(json_data["ranking_score"]))
        model_number = ranking_scores.index(max(ranking_scores))

        return model_number

    else:
        ranking_scores_csv = os.path.join(afOutputDir, "ranking_scores.csv")

        if not os.path.exists(ranking_scores_csv):
            return -1

        ranking_scores = []
        with open(ranking_scores_csv, 'r') as rankings_csv:
            reader = csv.DictReader(rankings_csv)
            for row in reader:
                print(row)
                ranking_scores.append(float(row["ranking_score"]))
        model_number = ranking_scores.index(max(ranking_scores))

        return model_number

