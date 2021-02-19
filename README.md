# Lung Tumor Evolution 
# "Emergence of a High-Plasticity Cell State during Lung Cancer Evolution", Nemanja Despot Marjanovic, Matan Hofree, Jason E. Chan, ..., Tyler Jacks, Aviv Regev, Tuomas Tammela, Cell, 2019

The repository contains code for reproducing the main manuscript figures. 
Including Matlab scripts for the following analysis:
* Basica single cell analysis workflow
    * Identify over-dispersed genes
    * NMF 
    * Graph cluster    
* NMI 
* eCDF normalization
* NMF-WOT
* ...

## Getting started
### Dependencies 
* Matlab (code was tested with R2020a)
* google cloud SDK
* (Optional) jupyter notebook with matlab kernel

### Dowload accompanying data
git clone https://github.com/matanhofree/lungTumorEvolution.git
cd lungTumorEvolution

cd data
gsutil -m cp -r gs://fc-0d580fef-66d6-47a8-afcf-45a2528719cb/lungTumorEvolution_data .
cd ..

### Main figure files:


# Under construction 2/19