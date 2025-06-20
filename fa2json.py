# Code for converting FASTA file to AlphaFold3 JSON format for running pairwise models
# Written by Miles Woodcock-Girard for Drew Lab at UIC


import os
import json
import string


# Loop over all FASTAs in current directory
for file in os.listdir("."):

    filename = os.fsdecode(file)

    if filename.endswith(".fasta") == False:
        continue

    else:
        uc_alphabet = list(string.ascii_uppercase)
        with open(file, 'r') as curr_fasta:
            fasta_data = curr_fasta.readlines()
            head_line = fasta_data[0]
            seq_line = fasta_data[1]

            job_name = filename.split(".")[0]
            head = fasta_data[::2]
            seqs = fasta_data[1::2]
            fasta_tuples = list(zip(head, seqs))
            

            json_dict = {
                    "name": job_name,
                    "modelSeeds": [1],
                    "sequences": [{"protein": {
                                     "id": uc_alphabet[i],
                                     "sequence": entry[1].strip(),
                                     }
                                  } for i, entry in enumerate(fasta_tuples)],
                    "dialect": "alphafold3",
                    "version": 1
                    }


            with open(f"./{job_name}.json", "w") as output_json:
                json.dump(json_dict, output_json)
