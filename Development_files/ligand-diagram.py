# ligand_clustering_pipeline.py

import requests
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import rdFMCS
from rdkit.DataStructs.cDataStructs import ConvertToNumpyArray
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import os
import json

CACHE_FILE = "smiles_cache.json"

PREDEFINED_SUBSTRUCTURES = {
    # Basic tetracyclic steroid nucleus
    "Steroid": "C1CCC2C(C1)CCC3C2CCC4=CC=CCC34",
    # Thiazolidinedione ring (5-membered ring with two carbonyls and a nitrogen)
    "TZD": "C1C(=O)NC(=O)S1",
    # Generic fatty acid chain with terminal carboxylic acid
    "FattyAcid": "C(=O)OCCCC"
}

# === 1. Fetch CACTVS Canonical SMILES from RCSB with caching ===
def load_smiles_cache():
    if os.path.exists(CACHE_FILE):
        with open(CACHE_FILE, 'r') as f:
            return json.load(f)
    return {}

def save_smiles_cache(cache):
    with open(CACHE_FILE, 'w') as f:
        json.dump(cache, f, indent=2)

def get_cactvs_smiles(ccd_code: str, cache: dict) -> str:
    if ccd_code in cache:
        return cache[ccd_code]

    url = f'https://data.rcsb.org/rest/v1/core/chemcomp/{ccd_code.upper()}'
    response = requests.get(url)
    if response.status_code != 200:
        return None
    try:
        smiles = response.json()['rcsb_chem_comp_descriptor']['smiles']
        cache[ccd_code] = smiles
        return smiles
    except KeyError:
        return None

# === 2. Read CCD codes from CSV ===
def read_ccd_codes(csv_path: str):
    return pd.read_csv(csv_path, header=None)[0].dropna().unique().tolist()

# === 3. Generate RDKit molecules ===
def build_molecules(ccd_codes):
    smiles_cache = load_smiles_cache()
    ccd_smiles = {}
    mols = []
    labels = []
    for code in ccd_codes:
        smiles = get_cactvs_smiles(code, smiles_cache)
        if smiles:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                ccd_smiles[code] = smiles
                mols.append(mol)
                labels.append(code)
            else:
                print(f"Failed to parse SMILES for {code}")
        else:
            print(f"Failed to retrieve SMILES for {code}")
    save_smiles_cache(smiles_cache)
    return mols, labels, ccd_smiles

# === 4. Group by predefined substructures ===
def group_by_substructure(mols, labels):
    substructs = {name: Chem.MolFromSmarts(smarts) for name, smarts in PREDEFINED_SUBSTRUCTURES.items()}
    group_map = {name: [] for name in substructs}
    uncategorized = []

    for mol, label in zip(mols, labels):
        found = False
        for name, patt in substructs.items():
            if mol.HasSubstructMatch(patt):
                group_map[name].append((mol, label))
                found = True
                break
        if not found:
            uncategorized.append((mol, label))

    group_map.update({f"Uncategorized:{lbl}": [(mol, lbl)] for mol, lbl in uncategorized})
    return group_map

# === 5. Draw molecules per group ===
def draw_groups(group_map, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    group_items = sorted(group_map.items(), key=lambda x: -len(x[1]))
    for i, (group, mol_label_list) in enumerate(group_items):
        mols, legends = zip(*mol_label_list)
        img = Draw.MolsToGridImage(mols, molsPerRow=5, legends=legends, subImgSize=(200, 200))
        safe_group = group.replace(':', '_').replace('/', '_')[:40]
        fname = f"group_{i:03d}_{safe_group}_n{len(mols)}.png"
        img.save(os.path.join(out_dir, fname))

# === 6. Optional t-SNE plot ===
def compute_fingerprints(mols):
    fps = [AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048) for mol in mols]
    fp_array = []
    for fp in fps:
        arr = np.zeros((1,), dtype=int)
        ConvertToNumpyArray(fp, arr)
        fp_array.append(arr)
    return np.array(fp_array)

def draw_tsne(fp_array):
    tsne = TSNE(n_components=2, metric='cosine', perplexity=30)
    tsne_result = tsne.fit_transform(fp_array)
    plt.figure(figsize=(10,8))
    plt.scatter(tsne_result[:,0], tsne_result[:,1], s=15, alpha=0.6)
    plt.title("t-SNE Projection of Ligands")
    plt.savefig("tsne_scatter.png")
    plt.close()

# === 7. Main pipeline ===
def main(csv_path: str, out_dir: str):
    print("Reading CCD codes...")
    ccd_codes = read_ccd_codes(csv_path)

    print("Fetching SMILES and building molecules...")
    mols, labels, smiles_dict = build_molecules(ccd_codes)

    print("Grouping ligands using predefined substructures...")
    group_map = group_by_substructure(mols, labels)

    print("Drawing all groups (sorted by size)...")
    draw_groups(group_map, out_dir)

    print("Optional: computing t-SNE plot...")
    fp_array = compute_fingerprints(mols)
    draw_tsne(fp_array)

    print("Done. Group images and t-SNE saved in:", out_dir)

# === Run ===
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Automatically group PPARg ligands by SMARTS-defined categories.")
    parser.add_argument("csv_path", help="Path to CSV file with CCD codes (one per line)")
    parser.add_argument("out_dir", help="Output directory for images")
    args = parser.parse_args()
    main(args.csv_path, args.out_dir)
