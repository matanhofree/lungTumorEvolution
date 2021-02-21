# Lung Tumor Evolution 
## "Emergence of a High-Plasticity Cell State during Lung Cancer Evolution"<br/>
## Nemanja Despot Marjanovic, Matan Hofree, Jason E. Chan, ..., Tyler Jacks, Aviv Regev, Tuomas Tammela, _Cell_, 2019<br/>
[10.1016/j.ccell.2020.06.012](https://doi.org/10.1016/j.ccell.2020.06.012)

The repository contains the princple computational methods and analysis scripts used in [1] 
Including code for the following analysis (Matlab):
* Basica single cell analysis workflow
    * Identify over-dispersed genes
    * NMF 
    * Graph cluster
    * tSNE and PHATE 2D representations   
* Genewise eCDF normalization
* NMI 
* NMF-WOT
* Single cell gene signatures 
* Consensus NMF and gsea of characteristic genes 

## Getting started
### Software requirements and dependencies 
* Matlab (code was tested with R2020a)
* R
* (Optional) Jupyter notebook with matlab kernel -- view and manipulate ".ipynb" notebooks
* (Optional) google cloud SDK

### Clone repository:
git clone https://github.com/matanhofree/lungTumorEvolution.git  
cd lungTumorEvolution  

### Download data:
* Login and download from the single cell portal:<br/>
https://singlecell.broadinstitute.org/single_cell/data/public/SCP971/emergence-of-a-high-plasticity-cell-state-during-lung-cancer-evolution-timecourse-smartseq2?filename=lungTE_data.tar.gz

* Download data from the public google bucket:<br/>
gsutil -u ADD GOOGLE USER PROJECT cp gs://regev-lung-tumorevo-mm-public/lungTE_data.tar.gz .<br/>
tar xzvf lungTE_data.tar.gz<br/>

### Main figure files:
* [Figure 1](https://github.com/matanhofree/lungTumorEvolution/blob/master/output/lung_tumor_evolution_Figure_1.ipynb) (Updated: 2/21/21)
* [Figure 2](https://github.com/matanhofree/lungTumorEvolution/blob/master/output/lung_tumor_evolution_Figure_2.ipynb) (Updated: 2/21/21)
* [Figure 3](https://github.com/matanhofree/lungTumorEvolution/blob/master/output/lung_tumor_evolution_Figure_3.ipynb) (Updated: 2/21/21)

#### Under construction (Last update 2/21/21)
