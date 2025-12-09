#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  5 16:39:23 2025

@author: frazaodesouzam2
"""

import mdtraj as md
import os
import numpy as np

def align_pdbs(pdb_files, ref_pdb, output_dir, core_residues):
    os.makedirs(output_dir, exist_ok=True)
    
    ref = md.load(ref_pdb)
    ref_core_indices = np.concatenate([
        ref.topology.select(f"chainid == 0 and name P and resid {r}") for r in core_residues
    ])
    
    for pdb in pdb_files:
        traj = md.load(pdb)
        core_indices = np.concatenate([
            traj.topology.select(f"chainid == 0 and name P and resid {r}") for r in core_residues
        ])
        
        # Check that core indices match
        if len(core_indices) != len(ref_core_indices):
            print(f"[WARN] Skipping {pdb}: core atom mismatch.")
            continue

        # Align using the core, but apply to all atoms
        traj.superpose(ref, atom_indices=core_indices, ref_atom_indices=ref_core_indices)
        
        # Save full aligned PDB
        out_name = os.path.basename(pdb).split(".")[0] + ".pdb"
        traj.save_pdb(os.path.join(output_dir, out_name))

        print(f"[INFO] Aligned and saved: {out_name}")
