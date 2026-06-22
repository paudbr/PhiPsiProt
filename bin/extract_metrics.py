#!/usr/bin/env python

import pickle
import os
import argparse
import json
#import torch moved to a conditional import since too bulky import if not used
import numpy as np
import csv
import string
from utils import plddt_from_struct_b_factor, get_chain_ids

# TODO: Issue #309, make into a proper separate process, it its own module so that dependencies can be managed better
# TODO: Need a sense of ranking, so that metrics can be traced back to correct model structure, even if they're not in sequential order. The enumerates() here are not sufficient.
#       Needs to be program-dependent, (see item below).
# TODO: look into have a --prog argument that could set filenames etc, logically seperate it?
# {name}_{prog}_{metric}.tsv might be easier for MultiQC to parse a complex workdir, than without the .prog
# TODO: read --prog from ${meta.model} in the NextFlow pipes. This also allows case switching in a proper EXTRACT_METRICS process.
# E.g. in main.nf of EXTACT_METRICS process, we could have:
# match ${meta.mode}:
#     case 'alphafold2':
#        ...
#     case 'rosettafold_all_atom':
#        ...
#...
# ^ overwrought with duplication, but can catch program specific weirdness, and lower barrier to adding new programs in the future.

# TODO: Chain-wise iPTM since the relevant interface might not always be the average of all.
# Would complete Issue #308
# Proposed format is pair-interfaces in rows, structure inference number in cols: https://github.com/nf-core/proteinfold/pull/312#issuecomment-2917709432
# KR - changed to have both sides of the matrix, because it's not symmetrical (see comment in Issue #306)

# Mapping of characters to integers for MSA parsing.
# 20 is for unknown characters, and 21 is for gaps.
AA_to_int = {
    "A": 0, "C": 1, "D": 2, "E": 3, "F": 4, "G": 5, "H": 6, "I": 7, "K": 8, "L": 9,
    "M": 10, "N": 11, "P": 12, "Q": 13, "R": 14, "S": 15, "T": 16, "V": 17, "W": 18, "Y": 19,
   ".": 20, "-": 21
}

def a3m_to_int(a3m_file):
    """
    Convert an A3M MSA representation into an integer representation (0-21).
    """
    with open(a3m_file, "r") as f:
        msa = f.read()

    # Convert each sequence in the MSA
    int_sequences = []
    for idx, line in enumerate(msa.splitlines()):
        if idx == 0 and not line.startswith(">"):  # If there's an additional header (non-FASTA) skip it. E.g ColabFold
            continue

        if not line.startswith(">"):  # Ignore header lines
            filtered_line = ''.join(char for char in line if not char.islower()) # Remove inserts (lower-case chars) in a3m
            int_sequence = [AA_to_int.get(char.upper(), 20) for char in filtered_line]
            int_sequences.append(int_sequence)

    int_sequences_array = np.array(int_sequences, dtype=object)
    return int_sequences_array

def format_msa_rows(msa_data):
    return [[str(x) for x in val] for val in msa_data]

def format_pae_rows(pae_data):
    return [[f"{num:.4f}" for num in row] for row in pae_data]

def format_iptm_rows(chain_pair_entries, chain_ids=None):
    """
    Format iPTM data into a list of rows for writing to a TSV file.
    Each row contains: the chain-pair in uppercase, e.g. "A:B", "B:A", A:C", etc. and then the iPTM value formatted to 4 decimal places.
    """
    def idx_to_letter(idx):
        """ Convert the index integer of the matrix to a letter representation that wraps to double representation, e.g. 0 -> A, 1 -> B, ..., 25 -> Z, 26 -> AA, 27 -> AB, etc.
            This is somewhat compatible with how protein structure chain names are numbered by biochemists.
            But we should move away from fixed-format PDB files -- we have nothing to lose but our chains."""
        result = ""
        while idx >= 0:
            result = string.ascii_uppercase[idx % 26] + result
            idx = idx // 26 - 1
            if idx < 0:
                break
        return result

    if chain_ids:
        #would be better with some model_id sorting
        iptm_rows = [[""]+[f"{chain_ids[idx[0]]}:{chain_ids[idx[1]]}" for idx, val in next(iter(chain_pair_entries.values()))]]
    else:
        iptm_rows = [[""]+[f"{idx_to_letter(idx[0])}:{idx_to_letter(idx[1])}" for idx, val in next(iter(chain_pair_entries.values()))]]

    for model_idx, chain_pair_entries_values in chain_pair_entries.items():
        iptm_rows.append([model_idx]+[f"{val:.4f}" for idx, val in chain_pair_entries_values])

    return [list(row) for row in zip(*iptm_rows)]


def chain_iptm_matrix_to_pairs(iptm_matrix):
    """
    Convert a chain-wise iPTM matrix to pair values by taking off-diagonal elements.
    """
    # From AlphaFold3 output docs:
    # 'chain_pair_iptm': An [num_chains, num_chains] array.
    # Off-diagonal element (i, j) of the array contains the ipTM restricted to tokens from chains i and j.
    # Diagonal element (i, i) contains the pTM restricted to chain i.
    return [(idx, val) for idx, val in np.ndenumerate(iptm_matrix) if idx[0] != idx[1]]

def chainwise_iptm_matrix_to_ptms(iptm_matrix):
    return [(idx, val) for idx, val in np.ndenumerate(iptm_matrix) if idx[0] == idx[1]]

def write_tsv(file_path, rows):
    with open(file_path, 'w') as out_f:
        writer = csv.writer(out_f, delimiter='\t')
        writer.writerows(rows)

def extract_structs_plddt_to_tsv(name, structures):
    """
    Write out a tsv file contain pLDDTs for reading by MultiQC in nf-core/proteinfold
    Uses utils function with BioPython PDB package to extract residue pLDDT values from the b-factor column.
    """
    plddt_cols = [plddt_from_struct_b_factor(structure) for structure in structures]
    res_counts = [len(plddt_col) for plddt_col in plddt_cols]

    if len(set(res_counts)) != 1:
        raise ValueError("Not all structures have the same number of residues!")

    rank_names = [f"rank_{i}" for i in range(len(structures))]
    # Create header as the first row
    plddt_rows =  [["Positions"] + rank_names]
    res_id_col = list(range(len(plddt_cols[0])))
    plddt_rows.extend(zip(res_id_col, *plddt_cols))  # Combine lists column-wise to make rows
    write_tsv(f"{name}_plddt.tsv", plddt_rows)

def read_pkl(name, pkl_files):
    """
    Adapted from the Galaxy AlphaFold tool (https://github.com/usegalaxy-au/tools-au/blob/de94df520c8dc7b8652aedb92e90f6ebb312f95f/tools/alphafold/scripts/outputs.py), originally authored by @neoformit and @graceahall and funded by Australian Biocommons and QCIF Australia.
    """
    ptm_data = {}
    iptm_data = {}
    for pkl_file in pkl_files:
        print(f"Processing {pkl_file}")
        data = pickle.load(open(pkl_file, "rb"))

        # Process MSA data
        if pkl_file.endswith("final_features.pkl"): # HelixFold3 - This one must be first
            write_tsv(f"{name}_msa.tsv", format_msa_rows(data["feat"]["msa"]))
        elif pkl_file.endswith("features.pkl"): # AlphaFold2.3
            try:
                N = data["num_alignments"][0] #monomer
            except:
                N = data["num_alignments"] #multimer
            write_tsv(f"{name}_msa.tsv", format_msa_rows(data["msa"][:N]))
        else:
            model_info = os.path.basename(pkl_file).replace("result_", "").replace(".pkl", "")
            #TODO: Make this explicit input
            with open(os.path.join(os.path.dirname(pkl_file),"ranking_debug.json")) as f:
                ranking_data = json.load(f)['order']
            model_id = ranking_data.index(model_info)
            if 'predicted_aligned_error' not in data.keys():
                print(f"No PAE output in {pkl_file}, it was likely a monomer calculation")
            else:
                write_tsv(f"{name}_{model_id}_pae.tsv", format_pae_rows(data["predicted_aligned_error"]))

            if 'ptm' not in data.keys():
                print(f"No pTM/iPTM output in {pkl_file}, it was likely a monomer calculation")
            else:
                ptm_data[model_id] = f"{np.round(data['ptm'],3)}\n"

            if 'iptm' in data:
                iptm_data[model_id] = f"{np.round(data['iptm'],3)}\n"
    if ptm_data:
        ptm_rows = sorted([[k, v.strip()] for k, v in ptm_data.items()], key=lambda x: x[0])
        write_tsv(f"{name}_ptm.tsv", ptm_rows)

    if iptm_data:
        iptm_rows = sorted([[k, v.strip()] for k, v in iptm_data.items()], key=lambda x: x[0])
        write_tsv(f"{name}_iptm.tsv", iptm_rows)

def read_paired_a3m(name, a3m_file):
    msa_rows = a3m_to_int(a3m_file)
    write_tsv(f"{name}_msa.tsv", format_msa_rows(msa_rows))

def read_a3m(name, a3m_files):
    # RosettaFold-All-Atom
    #TODO: DRY with unpaired below for Boltz
    msa_rows = {}
    for a3m_file in a3m_files: #Should already be alphabetical by chain
        msa_rows[a3m_file] = a3m_to_int(a3m_file)

    final_rows = []
    temp_row = []
    for a3m_file in a3m_files:
        temp_row.extend(msa_rows[a3m_file][0])
    final_rows.append(temp_row)

    # Un-paired TODO: get pairing code from RF-AA source
    # https://github.com/baker-laboratory/RoseTTAFold-All-Atom/blob/main/rf2aa/data/parsers.py#L405
    msa_widths = [len(msa_rows[chain][0]) for chain in a3m_files]
    msa_heights = [len(msa_rows[chain]) for chain in a3m_files]

    cum_total_rows = np.cumsum(msa_heights)
    for row_idx in range(cum_total_rows[-1]):
        temp_row = []

        for i, chain in enumerate(a3m_files):
            msa = msa_rows[chain]
            width = msa_widths[i]
            if i == 0:
                minrow = 0
            else:
                minrow = cum_total_rows[i-1]
            maxrow = cum_total_rows[i]

            if minrow <= row_idx < maxrow:
                msa_row_idx = row_idx - minrow
                temp_row.extend(msa[msa_row_idx])
            else:
                temp_row.extend(["21"] * width) #gap
        final_rows.append(temp_row)

    write_tsv(f"{name}_msa.tsv", format_msa_rows(final_rows))

def read_npz(name, npz_files):
   for idx, npz_file in enumerate(npz_files):
        data = np.load(npz_file)
       #Boltz PAE files if --write_full_pae is used
        if npz_file.split('/')[-1].startswith('pae') and npz_file.endswith('.npz'):
            model_id = os.path.basename(npz_file).split('_model_')[-1].split('.npz')[0]
            write_tsv(f"{name}_{model_id}_pae.tsv", format_pae_rows(data["pae"]))

# Boltz MSA processing
def read_csv(name, csv_files):
    if not csv_files or not os.path.isfile(csv_files[0]):
        return

    msa_rows = {}
    unpaired_msa_rows = {}

    for csv_file in sorted(csv_files, key=lambda x: int(x.split('_')[-1].split('.csv')[0])):
        msa_lines = []
        unpaired_msa_lines = []

        with open(csv_file) as f:
            f.readline()  # skip header

            for line in f:
                cols = line.strip("\n").split(",")

                if len(cols) < 2:
                    continue

                seq = ''.join(c for c in cols[1] if not c.islower())

                if cols[0] == "-1" and len(csv_files) > 1:
                    unpaired_msa_lines.append(seq)
                else:
                    msa_lines.append(seq)

        msa_id = csv_file.split("_")[-1].split(".csv")[0]

        msa_rows[msa_id] = [
            [str(AA_to_int.get(residue, 20)) for residue in line]
            for line in msa_lines
        ]

        unpaired_msa_rows[msa_id] = [
            [str(AA_to_int.get(residue, 20)) for residue in line]
            for line in unpaired_msa_lines
        ]

    with open(f'boltz_results_{name}/processed/manifest.json') as f:
        manifest = json.load(f)

    def get_msa_id(chain):
        msa_id = chain.get("msa_id", -1)

        if msa_id in (-1, None):
            return None

        msa_id = str(msa_id)

        if "_" in msa_id:
            msa_id = msa_id.split("_")[-1]

        return msa_id

    valid_chains = []

    for chain in manifest["records"][0]["chains"]:
        msa_id = get_msa_id(chain)

        if msa_id is None:
            continue

        if msa_id not in msa_rows:
            continue

        valid_chains.append((chain, msa_id))

    if len(valid_chains) == 0:
        print("WARNING: No valid MSA chains found. Skipping MSA TSV generation.")
        return

    final_rows = []

    # -----------------------
    # PAIRED MSA
    # -----------------------

    paired_height = min(
        len(msa_rows[msa_id])
        for _, msa_id in valid_chains
        if len(msa_rows[msa_id]) > 0
    )

    for row_idx in range(paired_height):
        temp_row = []

        for _, msa_id in valid_chains:
            temp_row.extend(msa_rows[msa_id][row_idx])

        final_rows.append(temp_row)

    # -----------------------
    # UNPAIRED MSA
    # -----------------------

    msa_widths = []

    for _, msa_id in valid_chains:
        if len(msa_rows[msa_id]) > 0:
            msa_widths.append(len(msa_rows[msa_id][0]))
        else:
            msa_widths.append(0)

    msa_heights = [
        len(unpaired_msa_rows.get(msa_id, []))
        for _, msa_id in valid_chains
    ]

    if sum(msa_heights) > 0:

        cum_total_rows = np.cumsum(msa_heights)

        for row_idx in range(cum_total_rows[-1]):

            temp_row = []

            for i, (_, msa_id) in enumerate(valid_chains):

                msa = unpaired_msa_rows.get(msa_id, [])
                width = msa_widths[i]

                minrow = 0 if i == 0 else cum_total_rows[i - 1]
                maxrow = cum_total_rows[i]

                if minrow <= row_idx < maxrow and len(msa) > 0:
                    msa_row_idx = row_idx - minrow
                    temp_row.extend(msa[msa_row_idx])
                else:
                    temp_row.extend(["21"] * width)

            final_rows.append(temp_row)

    write_tsv(f"{name}_msa.tsv", final_rows)

def read_json(name, json_files):
    ptm_data = {}
    iptm_data = {}
    chain_pair_entries = {}
    chainwise_ptms = {}
    chain_ids = []
    ranking_scores = {}

    for idx, json_file in enumerate(json_files):
        with open(json_file, 'r') as f:
            data = json.load(f)

        basename = os.path.basename(json_file)
        dirname  = os.path.dirname(json_file)

        # ── AF3: datos MSA ──────────────────────────────────────────
        if json_file.endswith("_data.json"):
            paired_msa_rows   = []
            unpaired_msa_rows = []
            for chain in data['sequences']:
                if 'protein' not in chain:
                    continue
                unpaired_MSA   = chain['protein'].get('unpairedMsa', '')
                unpaired_lines = [
                    ''.join(c for c in line if not c.islower())
                    for line in unpaired_MSA.split("\n")
                    if line.strip() and not line.startswith(">")
                ]
                unpaired_msa_rows.append([
                    [str(AA_to_int.get(r, 20)) for r in line]
                    for line in unpaired_lines
                ])
                paired_MSA   = chain['protein'].get('pairedMsa', '')
                paired_lines = [
                    ''.join(c for c in line if not c.islower())
                    for line in paired_MSA.split("\n")
                    if line.strip() and not line.startswith(">")
                ]
                paired_msa_rows.append([
                    [str(AA_to_int.get(r, 20)) for r in line]
                    for line in paired_lines
                ])

            chains     = len(paired_msa_rows)
            final_rows = []

            if chains > 0 and len(paired_msa_rows[0]) > 0:
                for i in range(len(paired_msa_rows[0])):
                    temp_row = []
                    for j in range(chains):
                        if j < len(paired_msa_rows) and i < len(paired_msa_rows[j]):
                            temp_row.extend(paired_msa_rows[j][i])
                        else:
                            width = len(paired_msa_rows[0][0]) if paired_msa_rows[0] else 0
                            temp_row.extend(["21"] * width)
                    final_rows.append(temp_row)

            msa_widths = []
            for chain_idx in range(chains):
                if paired_msa_rows[chain_idx]:
                    msa_widths.append(len(paired_msa_rows[chain_idx][0]))
                elif unpaired_msa_rows[chain_idx]:
                    msa_widths.append(len(unpaired_msa_rows[chain_idx][0]))
                else:
                    msa_widths.append(0)

            msa_heights = [len(unpaired_msa_rows[c]) for c in range(chains)]
            if sum(msa_heights) > 0:
                cum_total_rows = np.cumsum(msa_heights)
                for row_idx in range(cum_total_rows[-1]):
                    temp_row = []
                    for i in range(chains):
                        msa    = unpaired_msa_rows[i]
                        width  = msa_widths[i]
                        minrow = 0 if i == 0 else cum_total_rows[i - 1]
                        maxrow = cum_total_rows[i]
                        if minrow <= row_idx < maxrow and len(msa) > 0:
                            temp_row.extend(msa[row_idx - minrow])
                        else:
                            temp_row.extend(["21"] * width)
                    final_rows.append(temp_row)

            if final_rows:
                write_tsv(f"{name}_msa.tsv", final_rows)
            else:
                print(f"WARNING: Empty MSA for {name}, skipping MSA TSV.")

        # ── AF3: métricas de confianza ───────────────────────────────
        elif json_file.endswith("_summary_confidences.json"):
            model_id = 0
            if 'ranking_score' in data:
                ranking_scores[model_id] = data['ranking_score']
            if 'ptm' in data:
                ptm_data[model_id] = f"{np.round(data['ptm'], 3)}\n"
            if 'iptm' in data and data['iptm'] is not None:
                iptm_data[model_id] = f"{np.round(data['iptm'], 3)}\n"
            if 'chain_pair_iptm' in data:
                chain_iptm_matrix = np.array(data['chain_pair_iptm'])
                chain_pair_entries[model_id] = chain_iptm_matrix_to_pairs(chain_iptm_matrix)
                chainwise_ptms[model_id]     = chainwise_iptm_matrix_to_ptms(chain_iptm_matrix)
            if 'chain_pair_pae_min' in data:
                pae_min_matrix = data['chain_pair_pae_min']
                pae_min_rows   = [[""] + [f"chain_{j}" for j in range(len(pae_min_matrix[0]))]]
                for i, row in enumerate(pae_min_matrix):
                    pae_min_rows.append([f"chain_{i}"] + [f"{v:.4f}" for v in row])
                write_tsv(f"{name}_chain_pair_pae_min.tsv", pae_min_rows)
            extra = {k: data[k] for k in ('has_clash', 'fraction_disordered') if k in data}
            if extra:
                write_tsv(f"{name}_extra_metrics.tsv",
                          [list(extra.keys()), list(extra.values())])

        # ── AF3: PAE por modelo ──────────────────────────────────────
        elif json_file.endswith("_confidences.json"):
            model_id = 0
            if 'pae' in data:
                write_tsv(f"{name}_{model_id}_pae.tsv", format_pae_rows(data["pae"]))

        # ── Boltz1/Boltz2: confidence_*_model_*.json ────────────────
        elif 'predictions' in json_file or basename.startswith("confidence_"):
            model_id = basename.split('_model_')[-1].replace('.json', '')

            # ptm / iptm escalares — PRIMERO, antes de cualquier cosa que pueda fallar
            if 'ptm' in data:
                ptm_data[model_id] = f"{np.round(data['ptm'], 3)}\n"
            if 'iptm' in data and data['iptm']:
                iptm_data[model_id] = f"{np.round(data['iptm'], 3)}\n"

            # PAE
            if 'pae' in data:
                write_tsv(f"{name}_{model_id}_pae.tsv", format_pae_rows(data["pae"]))

            # chain_pair_iptm (Boltz2: pair_chains_iptm)
            if 'chain_pair_iptm' in data:
                chain_iptm_matrix = np.array(data['chain_pair_iptm'])
                chain_pair_entries[model_id] = chain_iptm_matrix_to_pairs(chain_iptm_matrix)
                chainwise_ptms[model_id]     = chainwise_iptm_matrix_to_ptms(chain_iptm_matrix)

            elif 'pair_chains_iptm' in data:
                cpd = data['pair_chains_iptm']
                chain_iptm_matrix = np.array([
                    [cpd[row][col] for col in sorted(cpd[row], key=int)]
                    for row in sorted(cpd, key=int)
                ])
                # Boltz2: confidence_NAME_model_0.json → NAME_model_0.pdb
                pdb_name = basename.replace("confidence_", "").replace(".json", ".pdb")
                pdb_path = os.path.join(dirname, pdb_name)
                if os.path.exists(pdb_path):
                    chain_ids = get_chain_ids(pdb_path)
                chain_pair_entries[model_id] = chain_iptm_matrix_to_pairs(chain_iptm_matrix)
                chainwise_ptms[model_id]     = chainwise_iptm_matrix_to_ptms(chain_iptm_matrix)

        else:
            # HF3 individual predictions
            if 'all_results' in json_file:
                model_id = int(os.path.dirname(json_file).split('-rank')[-1])
            else:
                model_id = idx

            if 'ptm' in data:
                ptm_data[model_id] = f"{np.round(data['ptm'], 3)}\n"
            if 'iptm' in data and data['iptm']:
                iptm_data[model_id] = f"{np.round(data['iptm'], 3)}\n"
            if 'pae' in data:
                write_tsv(f"{name}_{model_id}_pae.tsv", format_pae_rows(data["pae"]))

    # ── Escribir TSVs finales ────────────────────────────────────────
    if ranking_scores:
        write_tsv(f"{name}_ranking_scores.tsv",
                  [["model_id", "ranking_score"]] +
                  [[k, f"{v:.6f}"] for k, v in sorted(ranking_scores.items())])
    if chainwise_ptms:
        write_tsv(f"{name}_chainwise_ptm.tsv",
                  format_iptm_rows(chainwise_ptms, chain_ids=chain_ids))
    if chain_pair_entries:
        write_tsv(f"{name}_chainwise_iptm.tsv",
                  format_iptm_rows(chain_pair_entries, chain_ids=chain_ids))
    if ptm_data:
        ptm_rows = [[k, v.strip()] for k, v in sorted(ptm_data.items(), key=lambda x: str(x[0]))]
        write_tsv(f"{name}_ptm.tsv", ptm_rows)
    if iptm_data:
        iptm_rows = [[k, v.strip()] for k, v in sorted(iptm_data.items(), key=lambda x: str(x[0]))]
        write_tsv(f"{name}_iptm.tsv", iptm_rows)




def read_pt(name, pt_files):
    import torch # moved to a conditional import since too bulky import if not used
    #TODO: Handle this better when refactored - Is this just RFAA??
    for pt_file in pt_files:
        with open(pt_file, 'rb') as f:   # TODO: point to [protein]_aux.pt
            data = torch.load(f, map_location="cpu")
            if 'pae' in data:
                # The pt file contains a tensor that needs to be cast as an array
                # Squeeze leading dimension (batch?)
                write_tsv(f"{name}_0_pae.tsv", format_pae_rows(np.squeeze(data["pae"].numpy())))
        break

def read_colabfold_metrics(name, colabfold_metrics_fns):
    ptm_rows = []
    iptm_rows = []
    for fn in colabfold_metrics_fns:
        with open(fn) as f:
            data = json.load(f)
        rank_id = int(fn.split("rank_")[1].split("_")[0])-1
        model_id = int(fn.split("model_")[1].split("_")[0])
        seed_id = int(fn.split("seed_")[1].split(".")[0])
        if "pae" in data:
            write_tsv(f"{name}_{rank_id}_pae.tsv", format_pae_rows(data["pae"]))
        if "ptm" in data:
            ptm_rows.append((f"{rank_id}", data["ptm"]))
        if "iptm" in data:
            iptm_rows.append((f"{rank_id}", data["iptm"]))
    if len(ptm_rows)>0:
        write_tsv(f"{name}_ptm.tsv", sorted(ptm_rows, key = lambda x: x[0]))
    if len(iptm_rows)>0:
        write_tsv(f"{name}_iptm.tsv", sorted(iptm_rows, key = lambda x: x[0]))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pkls", dest="pkls", required=False, nargs="+") # For reading both HelixFold3 and AlphaFold2 MSA formats
    parser.add_argument("--npzs", dest="npzs", required=False, nargs="+") # For reading the Boltz-1 PAE formats. TODO: Boltz-1 MSA not implemented (go straight to .a3m file), implement
    parser.add_argument("--a3ms", dest="a3ms", required=False, nargs="+") # For reading the RosettaFold-All-Atom MSA formats
    parser.add_argument("--paired_a3m", dest="paired_a3m", required=False) # For reading the ColabFold MSA format
    parser.add_argument("--csvs", dest="csvs", required=False, nargs="+") # For reading boltz csvs
    parser.add_argument("--jsons", dest="jsons", required=False, nargs="+") # For reading the AF3 MSA & PAE, HF3 PAE
    parser.add_argument("--colabfold_metrics_fns", required=False, nargs="+")
    parser.add_argument("--pts", dest="pts", required=False, nargs="+") # For read RFAA pytorch model to get PAE data
    parser.add_argument("--structs", dest="structs", required=False, nargs="+")
    parser.add_argument("--name", default="untitled", dest="name") # might need a --name $meta.id
    args = parser.parse_args()

    if args.pkls:
        read_pkl(args.name, args.pkls)
    if args.a3ms:
        read_a3m(args.name, args.a3ms)
    if args.paired_a3m:
        read_paired_a3m(args.name, args.paired_a3m)
    if args.csvs:
        read_csv(args.name, args.csvs)
    if args.npzs:
        read_npz(args.name, args.npzs)
    if args.jsons:
        read_json(args.name, args.jsons)
    if args.pts:
        read_pt(args.name, args.pts)
    if args.structs:
        extract_structs_plddt_to_tsv(args.name, args.structs)
    if args.colabfold_metrics_fns:
        read_colabfold_metrics(args.name, args.colabfold_metrics_fns)

if __name__ == "__main__":
    main()
