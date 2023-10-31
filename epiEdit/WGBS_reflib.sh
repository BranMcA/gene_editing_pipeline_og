ml ib htslib/1.13 bwa/0.7.17
#conda activate py2

sample=$1
refFasta=$2 #"/hpc/grid/wip_drm_targetsciences/projects/p073_GeneEditingSupport/epigenomeEditing/InternalDataset/Sidney/AS33 Fastq/AS33_methylC_lib/A33_01_amplicon.fa"  #"/hpc/grid/wip_drm_targetsciences/projects/Ensembl/release-106/homo_sapients/GRCh38/dna/chromosomes/Homo_sapiens.GRCh38.dna.allChromosomes.fa"
refLibDir=$3 # "/hpc/grid/wip_drm_targetsciences/projects/p073_GeneEditingSupport/epigenomeEditing/InternalDataset/Sidney/AS33 Fastq/AS33_methylC_lib"
refPosGz="${sample}.pos.gz"
refConvFa="${sample}.conv.fa"

methyC='/hpc/grid/wip_drm_targetsciences/users/shangl02/github/methylCtools/methylCtools'

## step 1. build ref lib
cd "$refLibDir"
python "${methyC}" fapos "${refFasta}" - | bgzip > "${refPosGz}"
tabix -s 1 -b 2 -e 2 "${refPosGz}"

python ${methyC} faconv "${refFasta}" "${refConvFa}"
bwa index -a bwtsw "${refConvFa}"