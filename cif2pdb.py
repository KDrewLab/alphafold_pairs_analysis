import json
import os
import subprocess

for file in os.listdir():
    if not os.path.isdir(file):
        continue

    cif_file_list = [os.path.join(file, f"seed-1_sample-0/{file}_seed-1_sample-0_model.cif"),
                     os.path.join(file, f"seed-1_sample-1/{file}_seed-1_sample-1_model.cif"),
                     os.path.join(file, f"seed-1_sample-2/{file}_seed-1_sample-2_model.cif")]

    for i, cif in enumerate(cif_file_list):
        subprocess.run(["BeEM", cif, f"-p={file}/seed-1_sample-{i}/{file}"])
