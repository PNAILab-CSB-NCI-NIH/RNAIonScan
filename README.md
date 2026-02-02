# RNA Ion Scan: Ion Binding Mode Classifier for RNA

[![build](https://github.com/PNAILab-CSB-NCI-NIH/RNAIonScan/actions/workflows/ci.yml/badge.svg)](https://github.com/PNAILab-CSB-NCI-NIH/RNAIonScan/actions/workflows/ci.yml)  
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)  

RNAIonScan is a Python-based tool for identifying, classifying, and characterizing magnesium ion (Mg²⁺) behavior in RNA structural ensembles. The tool integrates Cryo-EM Q-scores, geometric annotation (CMM), and spatial clustering to identify stable (“rigid”) vs. dynamic (“flexible”) ions and to categorize coordination types (direct, hydrated, water-mediated).

This workflow is designed for large ensembles of aligned RNA structures obtained from cryo-EM, crystallography, modeling, or molecular simulations.
## Features
RNA Ion Scan automates the following tasks:
1. Mg²⁺ detection and coordination analysis: Identifies Mg²⁺ ions proximal to RNA atoms. Determines coordination residues within customizable distance cutoffs. Integrates Q-scores from cryo-EM maps to filter reliable ions. Merges ion annotations with CMM geometry classification.
2. Ion behavior classification: Distinguishes direct, hydrated, and water-mediated coordination. Resolves conflicts across frames to determine global coordination type.
3. Outlier detection and reclassification: Removes ions with incompatible coordination sets. Reassigns ions based on residue overlap. Moves small or noisy clusters to a reservoir.
4. Stability (RMSF) analysis: Computes RMSF per ion cluster. Distinguishes rigid vs. flexible coordination sites. Computes per-residue RMSF for coordinating nucleotides.
5. Spatial clustering using DBSCAN: Identifies robust ion-binding meta-sites across ensembles. Separates stable ion locations from transient ones.
6. Output generation: Generates CSV summaries of ion classifications.

## Installation

### 1. Clone the repository
```bash
git clone https://github.com/PNAILab-CSB-NCI-NIH/RNAIonScan.git
cd RNAIonScan
```

### 2. Create a Conda environment
```bash
conda create -n ion_scan python=3.11
```

### 3. Activate Environment
```bash
conda activate ion_scan
```

### 4. Install requirements with pip
```bash
pip install -r requirements
python -m ipykernel install --user --name ion_scan --display-name "RNAIonScan"
```

## Requirements

This package would require one to have the volumes and associated coordinates files (`.mrc` and `.pdb`). The files need to be organized as structured below (as shown in /data)
```
├── volumes/
│   ├── V1/
│   │   ├── V1.pdb
│   │   └── V1.mrc
│   ├── V2/
│   │   ├── V2.pdb
│   │   └── V2.mrc
│   ├── ...
│   ├── Vn/
│   │   ├── Vn.pdb
│   │   └── Vn.mrc
```

Once this is done, one can run the functions to align and compute the QScore. For running the RNAIonScan, one would also need to use [cmm-web-submitter](https://github.com/PNAILab-CSB-NCI-NIH/cmm-web-submitter) to extract the ion information. Note: The large mrc files are not included in the repository but are available at [Zenodo](https://zenodo.org/records/17832303).

## Usage

Once one have prepared the dataset as structured earlier, the following steps are requirements:
- a) Run [cmm-web-submitter](https://github.com/PNAILab-CSB-NCI-NIH/cmm-web-submitter) over this dataset to extract information from [CheckMyMetal](https://cmm.minorlab.org) website from all PDBs and corresponding density maps
- b) Align PDB files (see 1-pre_align.ipynb)
- c) Extract QScore (see 2-pre_qscore.ipynb)
- d) Run the pipeline (see 3-RNAIonScan.ipynb)

## Citing

If RNAIonScan helped your research, please cite:

Zenodo:
```bibtex
@software{ionscan2025,
  author = {Degenhardt, Maximilia F. S. and Degenhardt, Hermann F. and Wang, Yun-Xing},
  title = {RNAIonScan: Ion Binding Mode Classifier for RNA},
  year = {2025},
  publisher = {Zenodo},
  doi = {10.5281/zenodo.17832303},
  url = {https://doi.org/10.5281/zenodo.17832303},
}
```
