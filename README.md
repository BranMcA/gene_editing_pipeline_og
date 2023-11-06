# gene_editing_pipeline
This repo has a bundle of gene editing pipelines for nuclease editting, base editing , and epigenome editing


Nuclease Editting
1. Running env setting
   For historical reason, all sigularity container (*.sif) are stored on Target Science WIP (/hpc/grid/wip_drm_targetsciences/users/shangl02/github/ampliseq-baseEdit/singularity), copy the entire folder to your local branch root directory

2. How to run
main.sh experimentID fastqDataPath inputFormPath outputDirPath


Epigenome Editing
1. Running env setting
   MethylCTools is only compatible with Python2, thus build a python2 conda env and activate it
   For example:
    conda create --name py2env python=2.7
    conda activate py2env

2. How to run
   epiEdit/BisulfiteSeq_Pipeline.py inputFormPath outputDirPath fastqDataDirPath
 

