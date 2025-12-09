# Automated Web Submission Script for CMM Analysis

[![build](https://github.com/PNAILab-CSB-NCI-NIH/RNAIonScan/actions/workflows/ci.yml/badge.svg)](https://github.com/PNAILab-CSB-NCI-NIH/RNAIonScan/actions/workflows/ci.yml)  
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)  


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

### 3. Install requirements with pip

```bash
pip install -r requirements.txt
python -m ipykernel install --user --name ion_scan --display-name "RNAIonScan"
```

### 4. Install Playwright browser drivers

```bash
playwright install
```

### 5. Activate Environment
```bash
conda activate ion_scan
```

## 📄 License

MIT License. See `LICENSE` file for details.

## Citing

If CMM-Web-Submitter helped your research, please cite:

Zenodo:
```bibtex
@software{rnaionscan2025,
  author = {Degenhardt, Maximilia F. S. and Degenhardt, Hermann F. and Wang, Yun-Xing},
  title = {RNAIonScan: Ion Binding Mode Classifier for RNA},
  year = {2025},
  publisher = {Zenodo},
  doi = {10.5281/zenodo.17832303},
  url = {https://doi.org/10.5281/zenodo.17832303},
}
```
