'''
Code for calculating pDockQ for a given PDB or mmCIF model file
Written by Miles Woodcock-Girard for Drew Lab at UIC

Based on Patrick Bryant paper:

  Bryant, P., Pozzati, G. & Elofsson, A. Improved prediction of protein-protein interactions using AlphaFold2.
  Nat Commun 13, 1265 (2022). https://doi.org/10.1038/s41467-022-28865-w

Adapted from FoldDock repository:

  https://gitlab.com/ElofssonLab/FoldDock

'''


import argparse
import sys
import os
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser(description = 'Calculate pDockQ score for a predicted structure mmCIF file.')
parser.add_argument('--infile', nargs = 1, type = str, default = sys.stdin)


# Parse mmCIF model file
def parse_atom_record_cif(line):
    line = line.split()
    record = {}
    record['atom'] = line[2]
    record['atom_no'] = int(line[1])
    record['atom_name'] = line[3]
    record['res_name'] = line[5]
    record['chain'] = line[7]
    record['res_no'] = int(line[8])
    record['resid'] = line[15]
    record['x'] = float(line[10])
    record['y'] = float(line[11])
    record['z'] = float(line[12])
    record['occ'] = float(line[13])
    record['B'] = float(line[14])

    return record

# Parse PDB model file
def parse_atom_record_pdb(line):
    record = {}
    record['atom'] = line[0:6].strip()
    record['atom_no'] = int(line[6:11])
    record['atom_name'] = line[12:16].strip()
    record['atom_alt'] = line[17]
    record['res_name'] = line[17:20].strip()
    record['chain'] = line[21]
    record['res_no'] = int(line[22:26])
    record['insert'] = line[26].strip()
    record['resid'] = line[22:29]
    record['x'] = float(line[30:38])
    record['y'] = float(line[38:46])
    record['z'] = float(line[46:54])
    record['occ'] = float(line[54:60])
    record['B'] = float(line[60:66])

    return record
    
# Extract coordinate data from mmCIF or PDB file
def read_model_file(model_file):
    
    chain_coords, chain_plddt = {}, {}

    with open(model_file, 'r') as file:
        for line in file:
            if not line.startswith('ATOM'):
                continue

            # Parse differently depending on file type
            if model_file.endswith('.cif'):
                record = parse_atom_record_cif(line)
            elif model_file.endswith('.pdb'):
                record = parse_atom_record_pdb(line)

            # Extract coordinates, pLDDT for residues
            if record['atom_name'] == 'CB' or (record['atom_name'] == 'CA' and record['res_name'] == 'GLY'):
                if record['chain'] in [*chain_coords.keys()]:
                    chain_coords[record['chain']].append([record['x'], record['y'], record['z']])
                    chain_plddt[record['chain']].append(record['B'])
                else:
                    chain_coords[record['chain']] = [[record['x'], record['y'], record['z']]]
                    chain_plddt[record['chain']] = [record['B']]

    # Convert to arrays
    for chain in chain_coords:
        chain_coords[chain] = np.array(chain_coords[chain])
        chain_plddt[chain] = np.array(chain_plddt[chain])

    return chain_coords, chain_plddt


def calc_pdockq(chain_coords, chain_plddt, t):
    ch1, ch2 = [*chain_coords.keys()]
    coords1, coords2 = chain_coords[ch1], chain_coords[ch2]
    plddt1, plddt2 = chain_plddt[ch1], chain_plddt[ch2]

    # Calculate 2-norm
    mat = np.append(coords1, coords2, axis = 0)
    a_min_b = mat[:, np.newaxis, :] -mat[np.newaxis, :, :]
    dists = np.sqrt(np.sum(a_min_b.T ** 2, axis = 0)).T
    l1 = len(coords1)
    contact_dists = dists[:l1,l1:]
    contacts = np.argwhere(contact_dists<=t)

    if contacts.shape[0] < 1:
        pdockq = 0
        ppv = 0
    else:
        # Get the average interface pLDDT
        avg_if_plddt = np.average(np.concatenate([plddt1[np.unique(contacts[:,0])], plddt2[np.unique(contacts[:,1])]]))
        n_if_contact = contacts.shape[0]
        x = avg_if_plddt * np.log10(n_if_contact)

        # Define constants for pDockQ score
        L = 0.724
        x0 = 152.611
        k = 0.052
        b = 0.018

        # Evaluate pDockQ score for structure
        pdockq = L / (1 + np.exp(-1 * (k * (x - x0)))) + b

        # PPV
        PPV = np.array([0.98128027, 0.96322524, 0.95333044, 0.9400192 ,
            0.93172991, 0.92420274, 0.91629946, 0.90952562, 0.90043139,
            0.8919553 , 0.88570037, 0.87822061, 0.87116417, 0.86040801,
            0.85453785, 0.84294946, 0.83367787, 0.82238224, 0.81190228,
            0.80223507, 0.78549007, 0.77766077, 0.75941223, 0.74006263,
            0.73044282, 0.71391784, 0.70615739, 0.68635536, 0.66728511,
            0.63555449, 0.55890174])

        pdockq_thresholds = np.array([0.67333079, 0.65666073, 0.63254566, 0.62604391,
            0.60150931, 0.58313803, 0.5647381 , 0.54122438, 0.52314392,
            0.49659878, 0.4774676 , 0.44661346, 0.42628389, 0.39990988,
            0.38479715, 0.3649393 , 0.34526004, 0.3262589 , 0.31475668,
            0.29750023, 0.26673725, 0.24561247, 0.21882689, 0.19651314,
            0.17606258, 0.15398168, 0.13927677, 0.12024131, 0.09996019,
            0.06968505, 0.02946438])

        idx = np.argwhere(pdockq_thresholds >= pdockq)
        if len(idx) > 0:
            ppv = PPV[idx[-1]][0]
        else:
            ppv = PPV[0]

    return pdockq, ppv




# Main
if __name__ == "__main__":
    
    # Process command line arguments
    args = parser.parse_args()
    in_file = args.infile[0]

    # Obtain coordinates, plddt values from chains
    chain_coords, chain_plddt = read_model_file(in_file)

    # Ensure at least 2 chains exist for pDockQ calculation
    if len(chain_coords.keys()) < 2:
        print(f'Only one chain in structure file: {in_file}')
        sys.exit()

    # Calculate pDockQ score for structure
    t = 8 # Distance threshold, set at 8 Å
    pdockq, ppv = calc_pdockq(chain_coords, chain_plddt, t)

    print(f'pDockQ = {pdockq} for {in_file}')
    print(f'corresponds to PPV of at least {ppv}')
