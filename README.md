# HIV-phyloTSI-workshops
Materials for HIV-phyloTSI training workshops in Kigali and Gaborone.

The repository is structured as follows:
- *input*
- *output*
- *scripts*
- *Kigali* 
- *Gaborone*


## Prerequisites

To run the analyses, it is best to be on a Linux environment. 
It is possible to "simulate" such an environment on Windows, through WSL [tutorial](https://discourse.ubuntu.com/t/install-ubuntu-on-wsl2-on-windows-10/26269).
Once this is set up, it is possible to install the software required following the instructions [here](TODO), assuming R is installed.
In particular, the following pieces of software will be downloaded/installed:
- [phyloscanner](https://github.com/BDI-pathogens/phyloscanner)
- IQTREE 


## Running HIV-phyloTSI

Base on instructions given on the [HIV-phyloTSI github repository](https://github.com/BDI-pathogens/HIV-phyloTSI)
Specifically, need to:
    - clone the repo
    - install the python dependencies.
    - run the scripts

```{bash}
# move to the directory where you want to install HIV-phyloTSI
cd ~$HOME/git        # or any directory of your choice

# clone directories necessary to run the analysis
# BDIs code
git clone git@github.com/BDI-pathogens/HIV-phyloTSI.git
# workshop materials
git clone git@github.com/abriz97/HIV-phyloTSI-workshops.git

# install python dependencies through conda.
conda env create -f HIV-phyloTSI/hivphylotsi.yml
conda activate hivphylotsi

# run it on the input data 
python HIV-phyloTSI/HIV-phyloTSI.py \
    -d Model \
    -p HIV-phyloTSI-workshops/input/*PatStats.csv \ 
    -m HIV-phyloTSI-workshops/input/*maf.csv \ 
    -o HIV-phyloTSI-workshops/output/*.csv

# if more time, we can maybe discuss how to obtain the mafs?
```
