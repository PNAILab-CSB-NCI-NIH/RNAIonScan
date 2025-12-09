import os
import sys
import glob

SRC_DIR = os.path.dirname(__file__) + "/../src"
DATA_DIR = os.path.dirname(__file__) + "/../data"
sys.path.append(SRC_DIR)

from pathlib import Path
from align import align_pdbs
from qscore import score
from scan import ion_scan

PATH = "../"

def test_align():
    # Path settings
    pdb_pattern = f"{DATA_DIR}/V*/V*.cif"   # Pattern to obtain PDB files
    output_dir  = f"{DATA_DIR}/aligned"     # Output path for aligned PDBs
    reference   = f"{DATA_DIR}/V1/V1.cif"   # Alignment against this reference

    # Alignment settings
    # - Align PDBs through these core residues:
    core_residues = list(range(13, 16)) + list(range(20, 30)) + list(range(38, 85)) + \
                    list(range(97, 111)) + list(range(249, 262)) + list(range(273, 313)) + \
                    list(range(326, 341)) + list(range(386, 402))

    pdb_files = sorted(glob.glob(pdb_pattern))
    align_pdbs(
        pdb_files=pdb_files,
        ref_pdb=reference,
        output_dir=output_dir,
        core_residues=core_residues
    )

def test_QScore():
    # Path settings
    root_dir     = Path(f"{DATA_DIR}")
    chimera_path = '/Applications/ChimeraX_Daily.app/Contents/MacOS/ChimeraX'
    pdb_ext      = '.cif'
    volume_exts  = ['.mrc']

    if not os.path.exists(chimera_path): return

    qpaths = score(
        root_dir=root_dir,
        chimera_path=chimera_path,
        pdb_ext=pdb_ext,
        volume_exts=volume_exts
    )

    # Now, copy QScores to the prealigned folder
    target = f"{DATA_DIR}/aligned" # Output path for aligned PDBs
    for q in qpaths: os.system(f"cp {q} {target}")

def test_main():
    # Parameters
    path = f"{DATA_DIR}/aligned"
    CMM_results  = f"{DATA_DIR}/V*"
    ion_chain = "B"
    min_cluster_size = 2

    res = ion_scan(path, CMM_results, ion_chain, min_cluster_size)
    cmm = res["cmm"]
    df = res["rigid"]

    assert len(cmm) > 0
    assert len(df)  > 0
