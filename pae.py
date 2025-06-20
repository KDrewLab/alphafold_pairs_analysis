# Code for generating a predicted aligned error (PAE) plot for an AlphaFold3 confidence JSON
# Written by Miles Woodcock-Girard for Drew Lab at UIC

import json
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

import utils


def load_json(file_path):
    with open(file_path, 'r') as f:
        data = json.load(f)

        nres_c1 = data["token_chain_ids"].count("A")
        nres_c2 = data["token_chain_ids"].count("B")

    pae_data = np.array(data['pae'], dtype='float32')

    return pae_data, nres_c1, nres_c2

def plot_pae_heatmap(path, pae_data, nres_c1, nres_c2):
    plt.figure(figsize=(8, 6))
    #ax = sns.heatmap(pae_data, cmap='Greens_r', vmin=0, vmax=30, cbar=True)
    ax = sns.heatmap(pae_data, cmap='Greens_r', vmin=0, vmax=30, cbar=True)
    ax.axvline(x=nres_c1, linewidth=2, color="black")
    ax.axhline(y=nres_c1, linewidth=2, color="black")
    plt.title('Predicted Aligned Error (PAE) Heatmap')
    plt.xlabel('Residue Index')
    plt.ylabel('Residue Index')
    plt.savefig(path)
    plt.close()






for pair in os.listdir():
    print(pair)
    if not os.path.isdir(pair):
        continue

    if utils.usedAF3Server(pair):
        best_model = utils.getBestModel(pair)
        if best_model == -1:
            print("Output not familiar AF3 file structure. Skipping ...")
            continue
        conf_json = glob.glob(f"{pair}/fold_{pair}_*_full_data_{best_model}.json")[0]
    else:
        conf_json = f"{pair}/{pair}_confidences.json"
    
    pae_heatmap = f"pae/{pair}_pae.png"

    if not os.path.exists(conf_json):
        print(f"{pair} confidences JSON missing! Skipping.")
        continue

    if os.path.exists(pae_heatmap):
        print(f"{pair} PAE exists! Skipping.")
        continue

    pae_data, nres_c1, nres_c2 = load_json(conf_json)
    plot_pae_heatmap(pae_heatmap, pae_data, nres_c1, nres_c2)

'''
if __name__ == "__main__":
    json_path = "o95292_q9bz71_confidences.json"
    pair = "_".join(json_path.split("_")[0:1])
    pae_data, nres_c1, nres_c2 = load_json(json_path)
    plot_pae_heatmap(pair, pae_data, nres_c1, nres_c2)
'''
