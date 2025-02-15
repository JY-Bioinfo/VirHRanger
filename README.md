# **VirHRanger: Predicting animal host range from viral genomes and proteins using foundation models 🧬**



## Overview

This repository contains the code for implementing the ensemble model VirHRanger described in **"Developing Foundation Models for Predicting Viral Animal Host Range in Intelligent Surveillance"**, designed to predict the animal host range of viruses.

## Contents

- [Setup](#setup-)
- [HuggingFace](#huggingface-)
- [Input](#input-)
- [Options](#options-)
- [Usage](#usage-)
- [Citation](#citation-)
- [Troubleshooting](#troubleshooting-)
- [License](#license-)

## Setup ⚙️

### **Installation**

**1. Clone the Repository**:

   ```bash
   git clone https://github.com/JY-Bioinfo/VirHRanger
   cd VirHRanger
   ```

**2. Install Dependencies**

Create virtual environment and install required packages by running:

   ```bash
   conda env create -f environment.yml
   ```
Before each use, activate this environment using

```bash
   conda activate VirHRanger
   ```

**3. External Bioinformatics Tools**

Certain external bioinformatics tools are required depending on the features you wish to use:

- **Only use Genomic or Protein Foundation Model**: No external bioinformatics tools are needed.
- **Ensemble Model of Genomic Features**: The following bioinformatics tools are required:
  - BLAST
- **Ensemble Model of Genomic and Protein Features**: The following bioinformatics tools are required:
  - BLAST
  - CPB Machine
- **Ensemble Model of Genomic, Protein, and PPI Features**: The following bioinformatics tools are required:
  - BLAST
  - IUPred
  - SlimSuite
  - CPB Machine

Please note that the extraction of PPI (Protein-Protein Interaction) features is computationally intensive.

#### **BLAST**
We rely on BLAST (Basic Local Alignment Search Tool) to extract homology signals from viral genomes.

- **Setup Instructions**: To install BLAST locally, follow the instructions on NCBI
- Additionally, you may need to add BLAST executables to your system’s `PATH` to enable command-line usage.

#### **IUPred**
IUPred is a tool for predicting intrinsically disordered regions in proteins. Extraction of PPI features rely on the C or python2 version of IUPred (Requirements of SLimSuite).

- **Repository**: [IUPred Website](https://iupred2a.elte.hu/)
- **Setup Instructions**: Please download SLimSuite from the GitHub repository and put into the correct path.
- Additionally, you may need to add IUPred to your system’s `PATH` to enable command-line usage.

#### **SLimSuite**
SLimSuite is a software package of short linear motif (SLiM) discovery 

- **Repository**: [SlimSuite Repository](https://github.com/slimsuite/SLiMSuite)
- **Setup Instructions**: Please download SLimSuite from the GitHub repository and put into the correct path. For troubleshooting, please refer to the instruction of SLimProb on their official website (http://rest.slimsuite.unsw.edu.au/docs&page=module:slimprob)

#### **CPB Machine**
CPB Machine is a java program developed for extracting various features from viral coding sequences including nucleotide, dinucleotide and AA frequencies, to dinucleotide, codon and codon pair biases. These features are demonstrated useful in identifying the hosts of viruses in published research papers(DOI: 10.1126/science.aap9072)

- **Repository**: [CPB Machine Repository](https://github.com/rjorton/VirusFeatures)
- **Setup Instructions**: Please download the CPB_Machine.jar program from the GitHub repository and put into the correct path.


## HuggingFace 🤗

All models can be directly downloaded from the github repository, with the support from Git LFS. 
Due to the large file sizes of protein and genomic foundation models, We provide an alternative download option from the HuggingFace repository (https://huggingface.co/JYBioinfo/VirHRanger/).

## Input 🧬

We recommend to input complete viral genomes in GBFF format. (Represents nucleotide sequences, including metadata, annotation and the sequence itself) You can then process the input using the script
```bash
   bash Process_Input_from_Genbank.sh
   ```
If the ORF annotation is not available, we recommend to input nucleotide sequences in FASTA format and using the ensemble model of genomic features.

## Options 🔄

**VirHRanger** is an ensemble model using genomic, protein, and PPI Features. It comprises of six individual models, which can also be used stand alone. We also provide three ensemble model options depend on input limitations and computational budgets.

#### **Individual Models**
- Using viral genomes only
   - **Genomic Foundation Model (Genomic FM)**
   - **Homology-Based Model (Homology)**

- Using viral proteins only
   - **Protein Foundation Model (Protein FM)**
   - **Protein Trait-Based Model (Protein Traits)**

- Using viral genomes and proteins(or Open Reading Frame Annotation)
   - **Genomic Trait-Based Model (Genomic Traits)**
   - **PPI-Based Model (PPI)** 

#### **Ensemble Models**
- **Ensemble Model of Genomic Features**:
   - Genomic FM
   - Homology

- **Ensemble Model of Genomic and Protein Features**: 
   - Genomic FM
   - Protein FM
   - Homology
   - Genomic Traits
   - Protein Traits

- **(VirHRanger) Ensemble Model of Genomic, Protein, and PPI Features**
   - Genomic FM
   - Protein FM
   - Homology
   - Genomic Traits
   - Protein Traits
   - PPI

## Usage 🖥️

To run the algorithm, use the following command:

```bash
bash HostPredictor.sh \
      --input your_data \
      --output results
```

Make sure to replace `your_data` with the appropriate file or data you wish to process. The results will be saved in the `results` directory.

#### Example Command

```bash
bash HostPredictor_example.sh \
      -i Example_input/Bats/ \
      -o Out/ \
      -m Model/ \
      -c 16 \
      -e Genomic_Protein \
      -w whole \
      -g 0
```

#### Command Line Arguments

- **`-i`**: Input directory (folder of genomic or protein sequence data).
- **`-o`**: Output directory. The output will be named after the input folder name (results will be stored here).
- **`-m`**: Directory path containing the trained machine learning models.
- **`-c`**: Number of CPU cores to utilize (default is half of total cores).
- **`-e`**: Ensemble mode to use (options include Genomic, Genomic_Protein, Genomic_Protein_PPI).
- **`-f`**: Only effective when use one individual model (e.g., Genomic_FM, Protein_FM).
- **`-p`**: Processing mode of input files. (default is individual, meaning each viral record is stored in an individual folder. either "individual" or "batch").
- **`-w`**: Defines whether the entire process should run or if feature extraction can be skipped (options: "whole" or "only_predict").
- **`-g`**: GPU usage (set to a value between 0-7 for specific GPU or -1 for CPU).

---

## Citation 📜

If you have used VirHRanger in your research, please kindly cite the following publication:

@article {Guo2025.02.13.638012,
	author = {Guo, Jinyuan and Guo, Qian and Yin, Hengchuang and Han, Yilun and Geng, Peter X. and Hou, Jiaheng and Zhang, Haoyu and Tan, Jie and Li, Mo and Zhang, Yan and Jiang, Xiaoqing and Zhu, Huaiqiu},
	title = {Developing Foundation Models for Predicting Viral Animal Host Range in Intelligent Surveillance},
	elocation-id = {2025.02.13.638012},
	year = {2025},
	doi = {10.1101/2025.02.13.638012},
	publisher = {Cold Spring Harbor Laboratory},
	journal = {bioRxiv}
}


## Troubleshooting 🛠️

If you encounter issues:

- **Missing dependencies**: Ensure that all required dependencies (BLAST IUPred, SlimSuite, CPB Machine, and modified transformers package) are correctly installed or put into the right path.
- **CUDA errors**: Verify that your GPU setup is correctly configured. Set the correct GPU index in the command.

## License 📜

This project is licensed under the **Apache License 2.0**. See the [LICENSE](LICENSE) file for details.
