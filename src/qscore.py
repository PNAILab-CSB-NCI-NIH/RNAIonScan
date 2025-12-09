#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 29 15:50:54 2025

@author: frazaodesouzam2
"""

import os
import subprocess
from pathlib import Path

# Set your root directory path
root_dir = Path('/Users/frazaodesouzam2/Documents/Project/FL_RnaseP/Conf_Space_CryoEM/Apo_1mM_Mg/QScore')
chimera_path ='/Applications/ChimeraX-1.9.app/Contents/MacOS/ChimeraX'
pdb_ext = '.pdb'
volume_exts = ['.mrc']
def score(root_dir, chimera_path, pdb_ext, volume_exts):
    # Loop through each V* folder
    paths = root_dir.glob('V*')
    qpaths = []
    for subdir in paths:
        if not subdir.is_dir():
            continue

        # Find PDB and volume files
        pdb_files = list(subdir.glob(f'*{pdb_ext}'))
        volume_files = []  # <-- DEFINE volume_files HERE
        for ext in volume_exts:
            volume_files.extend(subdir.glob(f'*{ext}'))

        # Now you can use volume_files!
        if len(pdb_files) == 1 and len(volume_files) == 1:
            pdb_file = pdb_files[0]
            volume_file = volume_files[0]
            output_file = subdir / f'{pdb_file.stem}_qscore.csv'

            script = f"""
            open {volume_file}
            open {pdb_file}
            qscore :1-417 toVolume #1 useGui false referenceGaussianSigma 0.6 outputFile {output_file}
            exit
            """

            script_file = subdir / 'qscore_script.cxc'
            with open(script_file, 'w') as f:
                f.write(script.strip())

            result = subprocess.run(
                [chimera_path, '--nogui', '--script', str(script_file)],
                capture_output=True, text=True)

            if os.path.exists(output_file):
                qpaths.append(output_file)
                print(f"[INFO] {subdir.name} qscore saved: {output_file}")
            else:
                print(f"[ERROR] Operation failed {subdir.name}:")
                print("STDOUT:", result.stdout)
                print("STDERR:", result.stderr)

            script_file.unlink()
        else:
            print(f"Warning: {subdir.name} does not contain exactly one PDB and one volume file")
    
    return qpaths