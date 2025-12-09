#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 10 11:26:03 2025

@author: frazaodesouzam2
"""

import os
import numpy as np
import pandas as pd
from MDAnalysis import Universe
from MDAnalysis.lib.distances import distance_array
from sklearn.cluster import DBSCAN
from collections import defaultdict
import mdtraj as md
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import json

####### Coordination cutoffs
DIRECT_CUTOFF = 2.4   ### Cutoff for ions directly coordinating an RNA atom
WATER_CUTOFF = 3.4  ### Cutoff for identifying potential water molecules
HYDRATED_CUTOFF = 5.0  ### Cutoff for ions coordinated through a water bridge

Q_NUC_THRESHOLD = 0.5   ### nucleotide avarage Q-score  cutoff
Q_ION_THRESHOLD = 0.6   ### ion avarage Q-score  cutoff

CUTOFF = HYDRATED_CUTOFF  # general cutoff for selecting ions near RNA atom

ION_SELECTION = 'resname MG'
ATOM_SELECTION = (
    "nucleic and (name OP1 or name OP2 or "
    "name O2' or name O3' or name O4' or name O5' or "
    "name N1 or name N3 or name N7 or "
    "name O2 or name O4 or name O6)"
)

def load_cmm(path_CMM):
    site_info = []

    num = sorted(glob.glob(path_CMM))
    for i in num:
        f = f'{i}/CMM_results.json'
        with open(f, 'r') as file:
            cmm = json.load(file)
            pdb_id = os.path.basename(i)  # adjust if needed

            for site in cmm['sites']:
                chain, resid = site['ID'].split(':')
                site_info.append({
                    "pdb_id": pdb_id,
                    "chain": chain,
                    "resid": int(resid),
                    "geometry": site['Geometry'],
                    "ligands": ";".join(site['Ligands']),
                })

    site_df = pd.DataFrame(site_info)
    return site_df

def load_qscores(qscore_folder, pdb_filename, ion_chain):
    base_name = os.path.splitext(pdb_filename)[0]
    qscore_file = os.path.join(qscore_folder, f"{base_name}_qscore.csv")
    if not os.path.exists(qscore_file):
        print(f"[WARN] Missing Q-score file: {qscore_file}")
        return None, None

    try:
        df = pd.read_csv(qscore_file)
        df["Name"] = df["Name"].astype(str).str.strip()

        # RNA nucleotides chain A
        rna_df = df[df['Chain'] == 'A'].copy()
        rna_df = rna_df.rename(columns={"Number": "residue", "Qavg": "nucleotide_q"})
        rna_df.set_index("residue", inplace=True)

        # Ions in specified chain (B for apo, C for holo)
        ion_df = df[(df['Chain'] == ion_chain) & (df["Name"] == "MG")].copy()
        ion_df = ion_df.rename(columns={"Number": "residue", "Qavg": "ion_q"})
        ion_df.set_index("residue", inplace=True)

        return rna_df, ion_df
    except Exception as e:
        print(f"[ERROR] Failed to read Q-score file {qscore_file}: {e}")
        return None, None

def build_records(path, ion_chain, ion_selection, reference_selection, cutoff):
    pdb_files = sorted([f for f in os.listdir(path) if f.endswith('.pdb')])
    pdb_id = [os.path.splitext(f)[0] for f in pdb_files]
    ion_records = []
    for i, pdb_file in enumerate(pdb_files):
        pdb_path = os.path.join(path, pdb_file)
        u = Universe(pdb_path)
        mg_atoms = u.select_atoms(ion_selection)
        phosphate_atoms = u.select_atoms(reference_selection)

        if len(mg_atoms) == 0 or len(phosphate_atoms) == 0:
            print(f"[WARN] Skipping {pdb_file}: missing atoms.")
            continue

        rna_qdf, ion_qdf = load_qscores(path, pdb_file, ion_chain)

        for mg in mg_atoms:
            # --- Distance-based coordination ---
            min_dist = np.min(distance_array(phosphate_atoms.positions, mg.position[np.newaxis])[:, 0])
            classification = "Free"  # initial, will be refined later

            # --- Coordinating residues ---
            dists = distance_array(phosphate_atoms.positions, mg.position[np.newaxis])
            close_idxs = np.where(dists[:, 0] < cutoff)[0]
            if len(close_idxs) > 0:
                residues = phosphate_atoms[close_idxs].residues.resids
                unique_resids = sorted(set(residues))
            else:
                    unique_resids = []

            # --- Q-score for Mg ion ---
            ion_q = np.nan
            if ion_qdf is not None and mg.resid in ion_qdf.index:
                ion_q = ion_qdf.loc[mg.resid, 'ion_q']

            # --- Q-score for coordinating residues ---
            coord_qavg = np.nan
            if rna_qdf is not None and len(unique_resids) > 0:
                valid_res = [res for res in unique_resids if res in rna_qdf.index]
                if valid_res:
                    coord_qavg = rna_qdf.loc[valid_res, 'nucleotide_q'].mean()

            ion_records.append([
                mg.position[0], mg.position[1], mg.position[2], i, pdb_id[i],
                classification,
                ";".join(map(str, unique_resids)),
                mg.chainID,
                mg.resid,
                min_dist,
                ion_q,
                coord_qavg
                ])
    if not ion_records:
        raise ValueError("No Mg²⁺ ions found near RNA atoms.")
    return ion_records, pdb_files

def refined_coordination(row):
    geom = row['geometry']
    d = row['min_dist']

    if pd.isna(geom) or pd.isna(d):
        return 'unknown'

    if geom != 'Free' and geom != 'Poorly Coordinated':
        if d < DIRECT_CUTOFF:
            return 'direct'
    elif d < WATER_CUTOFF:
        return 'water'
    else:
        return 'hydrated'
    # For 'Free' or 'Poorly Coordinated' or unknown geometry
    if d < DIRECT_CUTOFF:
        return 'direct'
    elif d < WATER_CUTOFF:
        return 'water'
    else:
        return 'hydrated'

def build_df(ion_records, cmm_df):
    ion_array = np.round(np.array([rec[:4] for rec in ion_records], dtype=float), 3)
    pdb_names = [rec[4] for rec in ion_records]
    coord_types = [rec[5] for rec in ion_records]
    residues_str = [rec[6] for rec in ion_records]
    chains = [rec[7] for rec in ion_records]
    resids = [rec[8] for rec in ion_records]
    min_dists = [rec[9] for rec in ion_records]
    
    df =  pd.DataFrame({
        "x": ion_array[:, 0],
        "y": ion_array[:, 1],
        "z": ion_array[:, 2],
        "frame": ion_array[:, 3].astype(int),
        "pdb_id": pdb_names,
        "coordination": coord_types,
        "coordinating_residues": residues_str,
        "chain": chains,
        "resid": resids,
        "min_dist": min_dists,
        "ion_q": [rec[10] for rec in ion_records],
        "coord_nuc_qavg": [rec[11] for rec in ion_records],
        })
    df['pdb_id'] = df['pdb_id'].str.split('_').str[0]
    
    df = pd.merge(df, cmm_df, how='left', on=['pdb_id',  'resid'])
    df['coordination'] = df.apply(refined_coordination, axis=1)
    
    priority = {'direct': 0, 'hydrated': 1, 'water': 2}
    resid_coord_map = (
        df.groupby('resid')['coordination']
        .apply(lambda x: min(x, key=lambda c: priority.get(c, 99)))
        .to_dict()
        )
    df['coordination'] = df['resid'].map(resid_coord_map)
    return df

def filter_df(df):
    df = df[~df['coordination'].astype(str).isin(['unknown'])]
    df = df[(df['ion_q'] >= Q_ION_THRESHOLD) & (df['coord_nuc_qavg'] >= Q_NUC_THRESHOLD)]
    return df

def refine_df(df, min_cluster_size=8):
    resid_cluster_map = {
        resid: i for i, resid in enumerate(sorted(df["resid"].unique()))
    }
    df["resid_cluster"] = df["resid"].map(resid_cluster_map)
    refined_df, reservoir_df = split_outlier_ions_by_coordination(df)
    refined_df, reservoir_df = reassign_reservoir_ions(refined_df, reservoir_df)
    refined_df, reservoir_df = move_small_clusters_to_reservoir(refined_df, reservoir_df, min_cluster_size=min_cluster_size)
    refined_df = add_cluster_rmsf_column(refined_df)
    return refined_df

def cluster_df(df, path, files, min_cluster_size=8):
    clustered_rigid_df = apply_dbscan_to_rigid_df(df, eps=0.7, min_samples=min_cluster_size)
    rigid_df = clustered_rigid_df[clustered_rigid_df['dbscan_cluster'] != -1].copy()
    rigid_df = add_cluster_rmsf_column(rigid_df)

    chain_trajs = []
    pdb_files = sorted([f"{path}/{fi}" for fi in files])
    for pdb in pdb_files:
        traj = md.load(pdb)
        atom_indices = traj.topology.select('chainid == 0')#' and name P')
        if len(atom_indices) == 0:
            print(f"Warning: No atoms found for chain 0 in {pdb}")
            continue
        chain_traj = traj.atom_slice(atom_indices)
        chain_trajs.append(chain_traj)
    trajs = md.join(chain_trajs, check_topology=False)
    p_indices = np.arange(trajs.n_atoms)  # All atoms are P atoms
    rmsf_data = md.rmsf(trajs, reference=None, atom_indices=p_indices, precentered=True) * 10  # in Angstroms

    top = trajs.topology
    atom_to_resid = {i: top.atom(i).residue.resSeq for i in range(trajs.n_atoms)}

    # Create residue → RMSF mapping (take the average RMSF if multiple atoms per residue)
    res_rmsf_dict = defaultdict(list)
    for atom_idx, rmsf in enumerate(rmsf_data):
        resid = atom_to_resid[atom_idx]
        res_rmsf_dict[resid].append(rmsf)
    residue_rmsf = {resid: np.mean(vals) for resid, vals in res_rmsf_dict.items()}
    rigid_df["coord_residue_rmsf"] = rigid_df["coordinating_residues"].apply(
        lambda x: compute_group_residue_rmsf(x, residue_rmsf)
    )

    # Calculate cluster_rmsf using standard deviation (default)
    rigid_df = compute_cluster_rmsf(rigid_df)
    return rigid_df

def ion_scan(
        path, cmm, ion_chain='B', min_cluster_size=8,
        ion_selection=ION_SELECTION,
        reference_selection=ATOM_SELECTION,
        cutoff=CUTOFF):
    # Load CMM ions
    cmm_df = load_cmm(cmm)

    # Load PDBs/Ions
    ion_records, pdb_files = build_records(path, ion_chain, ion_selection, reference_selection, cutoff)

    # Build DF
    raw_df = build_df(ion_records, cmm_df)
    filt_df = filter_df(raw_df)
    ref_df = refine_df(filt_df, min_cluster_size)
    rigid_df = cluster_df(ref_df, path, pdb_files, min_cluster_size)

    # Output
    dfs = {
        "cmm": cmm_df,
        "coord": raw_df,
        "rigid":  rigid_df
    }
    return dfs

def split_outlier_ions_by_coordination(df):
    df = df.copy()

    # Convert coordinating_residues string to set of ints
    df['coord_set'] = df['coordinating_residues'].apply(
        lambda x: set(map(int, x.split(';'))) if isinstance(x, str) and x.strip() else set()
    )

    outlier_indices = []

    # Iterate over each cluster
    for cluster_id, group in df.groupby('resid_cluster'):
        for idx, row in group.iterrows():
            current_set = row['coord_set']
            # Compare with all other ions in the same cluster
            other_sets = group.drop(index=idx)['coord_set']

            has_overlap = any(len(current_set & other_set) > 0 for other_set in other_sets)
            if not has_overlap:
                outlier_indices.append(idx)

    # Move outliers to reservoir_df
    reservoir_df = df.loc[outlier_indices].drop(columns='coord_set').copy()
    refined_df = df.drop(index=outlier_indices).drop(columns='coord_set')

    return refined_df, reservoir_df

def reassign_reservoir_ions(refined_df, reservoir_df):
    # Prepare coordinate sets
    refined_df = refined_df.copy()
    reservoir_df = reservoir_df.copy()
    
    refined_df['coord_set'] = refined_df['coordinating_residues'].apply(
        lambda x: set(map(int, x.split(';'))) if isinstance(x, str) and x.strip() else set()
    )
    reservoir_df['coord_set'] = reservoir_df['coordinating_residues'].apply(
        lambda x: set(map(int, x.split(';'))) if isinstance(x, str) and x.strip() else set()
    )

    reassigned_rows = []
    keep_in_reservoir = []

    for idx, res_row in reservoir_df.iterrows():
        res_set = res_row['coord_set']
        best_cluster = None

        # Check against all clusters in refined_df
        for cluster_id, group in refined_df.groupby('resid_cluster'):
            overlap_found = any(
                len(res_set & ref_set) >= 2
                for ref_set in group['coord_set']
            )
            if overlap_found:
                best_cluster = cluster_id
                break  # assign to first matching cluster

        if best_cluster is not None:
            # Update the cluster label and store for reassignment
            new_row = res_row.copy()
            new_row['resid_cluster'] = best_cluster
            reassigned_rows.append(new_row)
        else:
            keep_in_reservoir.append(res_row)

    # Update DataFrames
    reassigned_df = pd.DataFrame(reassigned_rows).drop(columns='coord_set')
    reservoir_df_final = pd.DataFrame(keep_in_reservoir).drop(columns='coord_set')
    refined_df_final = pd.concat([refined_df.drop(columns='coord_set'), reassigned_df], ignore_index=True)

    return refined_df_final, reservoir_df_final

def move_small_clusters_to_reservoir(refined_df, reservoir_df, min_cluster_size=8):
    # Count ions per cluster
    cluster_sizes = refined_df['resid_cluster'].value_counts()

    # Identify small clusters
    small_clusters = cluster_sizes[cluster_sizes < min_cluster_size].index

    # Move them to reservoir
    small_cluster_df = refined_df[refined_df['resid_cluster'].isin(small_clusters)].copy()
    refined_df = refined_df[~refined_df['resid_cluster'].isin(small_clusters)].copy()

    # Concatenate to reservoir_df
    reservoir_df = pd.concat([reservoir_df, small_cluster_df], ignore_index=True)

    return refined_df, reservoir_df

def add_cluster_rmsf_column(df):
    df = df.copy()
    rmsf_per_cluster = {}

    # Compute RMSF for each cluster
    for cluster_id, group in df.groupby('resid_cluster'):
        coords = group[['x', 'y', 'z']].to_numpy()
        if len(coords) < 2:
            rmsf = np.nan  # or 0.0 if you prefer
        else:
            center = coords.mean(axis=0)
            displacements = np.linalg.norm(coords - center, axis=1)
            rmsf = displacements.std()
            #rmsf = np.sqrt((displacements ** 2).mean())  # RMSF from center

        rmsf_per_cluster[cluster_id] = rmsf

    # Map the computed RMSF back to each row
    df['cluster_rmsf'] = df['resid_cluster'].map(rmsf_per_cluster)
    return df

def apply_dbscan_to_rigid_df(df, eps=1.5, min_samples=2):
    """
    Apply DBSCAN clustering to rigid_df based on x, y, z coordinates.
    
    Parameters:
        df : pandas.DataFrame
            DataFrame containing columns 'x', 'y', 'z'
        eps : float
            Maximum distance between two samples for one to be considered as in the neighborhood of the other.
        min_samples : int
            Minimum number of samples in a neighborhood to form a cluster.

    Returns:
        clustered_df : pandas.DataFrame
            Copy of df with an added column 'dbscan_cluster' containing the cluster labels.
    """
    df = df.copy()

    # Extract coordinate array
    coords = df[['x', 'y', 'z']].values

    # Apply DBSCAN
    db = DBSCAN(eps=eps, min_samples=min_samples)
    labels = db.fit_predict(coords)

    # Add cluster labels to DataFrame
    df['dbscan_cluster'] = labels

    return df

def dbscan_cluster(df):
    df = df.copy()
    rmsf_per_cluster = {}

    # Compute RMSF for each cluster
    for cluster_id, group in df.groupby('dbscan_cluster'):
        coords = group[['x', 'y', 'z']].to_numpy()
        if len(coords) < 2:
            rmsf = np.nan  # or 0.0 if you prefer
        else:
            center = coords.mean(axis=0)
            displacements = np.linalg.norm(coords - center, axis=1)
            rmsf = displacements.std()
            #rmsf = np.sqrt((displacements ** 2).mean())  # RMSF from center

        rmsf_per_cluster[cluster_id] = rmsf

    # Map the computed RMSF back to each row
    df['cluster_rmsf'] = df['resid_cluster'].map(rmsf_per_cluster)
    return df

def compute_group_residue_rmsf(coord_str, rmsf_lookup):
    try:
        residues = [int(r) for r in coord_str.split(";")]
        values = [rmsf_lookup.get(r, np.nan) for r in residues]
        values = [v for v in values if not np.isnan(v)]
        return np.mean(values) if values else np.nan
    except:
        return np.nan

def compute_cluster_rmsf(df, coord_cols=['x', 'y', 'z'], method='std'):
    """
    Computes RMSF per DBSCAN cluster and adds a 'cluster_rmsf' column to the DataFrame.

    Parameters:
    -----------
    df : pd.DataFrame
        Must contain 'dbscan_cluster' and 3D coordinate columns.

    coord_cols : list
        Columns containing x, y, z coordinates (default: ['x', 'y', 'z'])

    method : str
        'std' for standard deviation (default), or 'rms' for root-mean-square displacement.

    Returns:
    --------
    pd.DataFrame with updated 'cluster_rmsf' column
    """
    df = df.copy()
    cluster_rmsf_map = {}

    for cluster_id, group in df.groupby("dbscan_cluster"):
        coords = group[coord_cols].values
        if len(coords) < 2:
            cluster_rmsf_map[cluster_id] = 0.0
        else:
            center = coords.mean(axis=0)
            displacements = np.linalg.norm(coords - center, axis=1)
            if method == 'rms':
                rmsf = np.sqrt((displacements ** 2).mean())
            else:
                rmsf = displacements.std()
            cluster_rmsf_map[cluster_id] = rmsf

    # Map to new column
    df['cluster_rmsf'] = df['dbscan_cluster'].map(cluster_rmsf_map)
    return df
