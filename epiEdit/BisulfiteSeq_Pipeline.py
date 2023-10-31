'''
This is the bisulfite-seq pipeline analyzed with methyLCtools.
use python2 which is required by methylC package
'''

import pandas as pd
import os
import sys
import subprocess

def run(*popenargs, **kwargs):
    input = kwargs.pop("input", None)
    check = kwargs.pop("handle", False)

    if input is not None:
        if 'stdin' in kwargs:
            raise ValueError('stdin and input arguments may not both be used.')
        kwargs['stdin'] = subprocess.PIPE

    process = subprocess.Popen(*popenargs, **kwargs)
    try:
        stdout, stderr = process.communicate(input)
    except:
        process.kill()
        process.wait()
        raise
    retcode = process.poll()
    if check and retcode:
        raise subprocess.CalledProcessError(
            retcode, process.args, output=stdout, stderr=stderr)
    return retcode, stdout, stderr

def prepareRefLib(sample, seq, dir):
    """
    create reflib of sample 
    """
    reflib = os.path.join(dir, "reflib", sample)
    if not os.path.exists(reflib):
        os.makedirs(reflib)

    ofile = os.path.join(reflib, sample+".amplicon.fa")
    o = open(ofile, "w") 
    o.write(">" + sample + "_amplicon\n" +seq+"\n")
    o.close()
    return ofile

def buildRefLib(sample, refFa, refLibDir):
    script=os.path.join(os.path.dirname(os.path.realpath(__file__)), 'WGBS_reflib.sh')
    #cmd="bash \"{}\" {} \"{}\" \"{}\"".format(script, sample, refFa, refLibDir)
    #print(cmd)
    run(["bash", script, sample, refFa, refLibDir])


def runMethylC(sample, fqdir, rootOutdir, reflibdir):
    script=os.path.join(os.path.dirname(os.path.realpath(__file__)), 'WGBS_pipeline.sh')
    outdir = os.path.join(rootOutdir, "methylC", sample)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    run(["bash", script, sample, fqdir, outdir, reflibdir])

def concatCall(rootOutdir, dedup):
    dir=os.path.join(rootOutdir, "methylC")
    ofile = os.path.join(dir, 'summary{0}.call'.format(".dedup" if dedup else ""))
    o=open(ofile, "w")
    o.write('chromosome\tcoordinate\tstrand\tcontext\tunconvertedCytosines\tconvertedCytosines\n')
    o.close()
    cmd='for i in {0}/*/*{1}call.gz; do zcat $i >> {2}; done'.format(dir, "dedup." if dedup else "[!d][!e][!d][!u][!p].", ofile)
    #print(cmd)
    subprocess.call(cmd, shell=True)
    #run(["eval"], cmd)
    return(ofile)

def addCol(file, formFile):
    # summaryFile = r'X:\projects\p073_GeneEditingSupport\epigenomeEditing\InternalDataset\ProcessedData\AS38\methylC\summary.call'
    # formFile = r'X:\projects\p073_GeneEditingSupport\epigenomeEditing\InternalDataset\Sidney\AS38 CRISPRoff mRNA delivery H2B D4+D28, CD81 D14+d28.csv'
    df_summary = pd.read_csv(file, sep='\t')
    df_summary['sampleID'] = df_summary['chromosome'].str.replace('_amplicon', '')
    df_form = pd.read_csv(formFile, usecols=['Sample Number', 'Sample Name'])
    merge_df = pd.merge(df_form, df_summary, left_on='Sample Number', right_on='sampleID', how='left')
    merge_df.drop('sampleID', axis=1, inplace=True)
    merge_df.to_csv(file, sep='\t', index=False)
    

# formFile="/hpc/grid/wip_drm_targetsciences/projects/p073_GeneEditingSupport/epigenomeEditing/InternalDataset/Sidney/AS48/AS48_EMseq_test.modified.csv"
# outdir='/hpc/grid/wip_drm_targetsciences/projects/p073_GeneEditingSupport/epigenomeEditing/InternalDataset/ProcessedData/AS48'
# rawDataDir='/hpc/grid/wip_drm_targetsciences/projects/p073_GeneEditingSupport/epigenomeEditing/InternalDataset/Sidney/AS48/fastq'

formFile = sys.argv[1] ##"/hpc/grid/wip_drm_targetsciences/projects/p073_GeneEditingSupport/epigenomeEditing/InternalDataset/Sidney/AS35 D4 + D7 mRNA delivery and AS33 re-run CRISPRoff BisSeq submission form.csv"
outdir = sys.argv[2] ##'/hpc/grid/wip_drm_targetsciences/projects/p073_GeneEditingSupport/epigenomeEditing/InternalDataset/ProcessedData/AS35.dedup'
rawDataDir = sys.argv[3] ##'/hpc/grid/wip_drm_targetsciences/projects/p073_GeneEditingSupport/epigenomeEditing/InternalDataset/Sidney/AS35 Fastq/AS35 fastq'

df = pd.read_csv(formFile)

print(df)
for index, row in df.iterrows():
    try :
        sample = df.loc[index,"Sample Number"]
        refseq = df.loc[index, "Full Amplicon"]

        print("prepare reflib directory of " + sample)
        # refFa = prepareRefLib(sample, refseq, outdir)

        print("build reflib of " + sample)
        # reflibDir=os.path.dirname(refFa)
        # buildRefLib(sample, refFa, reflibDir)

        print("run methyCtools of " + sample)
        # runMethylC(sample, rawDataDir, outdir, reflibDir)
    except(OSError,ValueError, TypeError):
        print("Fail to process row {0}".format(index))
        break

print("collect un-dedup call to generate summary")
summaryFile=concatCall(outdir, False)
addCol(summaryFile, formFile)

print("collect dedup call to generate summary")
summaryFile=concatCall(outdir, True)
addCol(summaryFile, formFile)
