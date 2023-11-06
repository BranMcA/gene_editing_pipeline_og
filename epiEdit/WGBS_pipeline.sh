ml ib htslib/1.13 bwa/0.7.17 samtools/1.13
ml eb/2017  GCCcore/.5.4.0 picard/2.1.0-Java-1.8.0_92
#conda activate py2


sample=$1
rawDataDir=$2 #'/hpc/grid/wip_drm_targetsciences/projects/p073_GeneEditingSupport/epigenomeEditing/InternalDataset/Sidney/AS33 Fastq/AS33.fq'
outDir=$3 #'/hpc/grid/wip_drm_targetsciences/projects/p073_GeneEditingSupport/epigenomeEditing/InternalDataset/ProcessedData'
refLibDir=$4 #"${outDir}/refLib/${sample.reflib}" #/hpc/grid/wip_drm_targetsciences/projects/p073_GeneEditingSupport/epigenomeEditing/InternalDataset/Sidney/AS33 Fastq/AS33_methylC_lib" #'/hpc/grid/wip_drm_targetsciences/projects/Resource/methylCtools_reflib/GRCh38.ensembl106'

methyC='/hpc/grid/wip_drm_targetsciences/users/shangl02/github/methylCtools/methylCtools'
refConvFa="${sample}.conv.fa" #'human.GRCh38_ensemb106.conv.fa'
refPosGz="${sample}.pos.gz"
picard="/hpc/grid/wip_drm_targetsciences/users/shangl02/Software/picard/picard_2.25.0/picard.jar"

## Perform per sample
echo "fqconv step"

fqs=($(ls "${rawDataDir}" -1 | grep "^${sample}.*\(fastq\.gz\|fq\.gz\)$"))
# len=${#fqs[@]}
echo $sample ${fqs[@]}
if [ ${#fqs[@]} -ne 2 ]; then
    echo "Paired fastq files of ${sample} are not found in ${rawDataDir}. Now exit."
    exit 1
fi 
# python ${methyC} fqconv -1 "${rawDataDir}/${sample}_L001_R1_001.fastq.gz" -2 "${rawDataDir}/${sample}_L001_R2_001.fastq.gz" "${outDir}/${sample}.conv.fq"
python ${methyC} fqconv -1 "${rawDataDir}/${fqs[0]}" -2 "${rawDataDir}/${fqs[1]}" "${outDir}/${sample}.conv.fq"

############################################################################################################
####################### This is the old section to be removd after #########################################
##bsub -q long -M 131072 -R "rusage[mem=131072]" "
#bwa mem -p -M "${refLibDir}/${refConvFa}" "${outDir}/${sample}.conv.fq" >"${outDir}/${sample}.conv.sam";

# remove duplicate in SAM header
#DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
#echo $DIR
#python $DIR/RemoveSamDupHeader.py ${outDir}/${sample}.conv.sam;
############################################################################################################

echo "bwa step"
bwa mem -p -M "${refLibDir}/${refConvFa}" "${outDir}/${sample}.conv.fq" | samtools view -Sb - > "${outDir}/${sample}.conv.bam"

echo "bconv step"
python ${methyC} bconv "${outDir}/${sample}.conv.bam" "${outDir}/${sample}.conv.bconv.bam"
samtools sort "${outDir}/${sample}.conv.bconv.bam" -o "${outDir}/${sample}.conv.bconv.sorted.bam"
samtools index "${outDir}/${sample}.conv.bconv.sorted.bam"

echo "bcall step"
python ${methyC} bcall "${refLibDir}/${refPosGz}" "${outDir}/${sample}.conv.bconv.sorted.bam" - | bgzip > ${outDir}/${sample}.call.gz
tabix -s 1 -b 2 -e 2 ${outDir}/${sample}.call.gz 

echo "picard dedup"
java -jar ${picard}  MarkDuplicates I="${outDir}/${sample}.conv.bconv.sorted.bam" O="${outDir}/${sample}.conv.bconv.sorted.dedup.bam" M="${outDir}/${sample}.dedupMetric.txt"
samtools index "${outDir}/${sample}.conv.bconv.sorted.dedup.bam"

echo "bcall step (on dedup)"
python ${methyC} bcall "${refLibDir}/${refPosGz}" "${outDir}/${sample}.conv.bconv.sorted.dedup.bam" - | bgzip > ${outDir}/${sample}.dedup.call.gz
tabix -s 1 -b 2 -e 2 ${outDir}/${sample}.dedup.call.gz
