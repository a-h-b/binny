import os
import sys
import shutil
import gzip
#import json
import yaml
import bz2
import re
from copy import deepcopy
import subprocess
import pandas as pd
from pathlib import Path
import urllib.request

def open_output(filename):
    return(open(OUTPUTDIR+'/'+filename, 'w+'))

# default executable for snakmake
shell.executable("bash")

# default configuration file
configfile:
    srcdir("config/config.default.yaml")

# some parameters
SRCDIR = srcdir("workflow/scripts")
BINDIR = srcdir("workflow/bin")
ENVDIR = srcdir("workflow/envs")

# get parameters from the config file
# output
if os.path.isabs(os.path.expandvars(config['outputdir'])):
    OUTPUTDIR = os.path.expandvars(config['outputdir'])
else:
    OUTPUTDIR = os.getcwd() + "/" + os.path.expandvars(config['outputdir'])

# input
if os.path.isabs(os.path.expandvars(config['raws']['Contigs'])):
    CONTIGS = os.path.expandvars(config['raws']['Contigs'])
else:
    CONTIGS = os.getcwd() + "/" + os.path.expandvars(config['raws']['Contigs'])
# Added depth file par to us instead of alignment
if config['raws']['Contig_depth']:
    if os.path.isabs(os.path.expandvars(config['raws']['Contig_depth'])):
        CONTIG_DEPTH = os.path.expandvars(config['raws']['Contig_depth'])
    else:
        CONTIG_DEPTH = os.getcwd() + "/" + os.path.expandvars(config['raws']['Contig_depth'])
else:
    CONTIG_DEPTH = None
    if os.path.isabs(os.path.expandvars(config['raws']['Alignment_metagenomics'])):
        MGaln = os.path.expandvars(config['raws']['Alignment_metagenomics'])
    else:
        MGaln = os.getcwd() + "/" + os.path.expandvars(config['raws']['Alignment_metagenomics'])

SAMPLE = config['sample']
if SAMPLE == "":
    SAMPLE = "_".join(OUTPUTDIR.split("/")[-2:])
SAMPLE = re.sub("_+","_",re.sub("[;|.-]","_",SAMPLE))
DBPATH = os.path.expandvars(config['db_path'])
if not os.path.isabs(DBPATH):
    DBPATH = os.getcwd() + "/" + DBPATH
if not os.path.exists(DBPATH):
    os.makedirs(DBPATH)
    urllib.request.urlretrieve("https://webdav-r3lab.uni.lu/public/R3lab/IMP/essential.hmm", DBPATH + "/essential.hmm")

# hardware parameters
MEMCORE = str(config['mem']['normal_mem_per_core_gb']) + "G"
BIGMEMCORE = str(config['mem']['big_mem_per_core_gb']) + "G"

# temporary directory will be stored inside the OUTPUTDIR directory
# unless a absolute path is set
TMPDIR = config['tmp_dir']
if not os.path.isabs(TMPDIR):
    TMPDIR = os.path.join(OUTPUTDIR, TMPDIR)
if not os.path.exists(TMPDIR):
    os.makedirs(TMPDIR)

# set working directory and dump output
workdir:
    OUTPUTDIR


def prepare_input_files(inputs, outputs):
    """
    Prepare file names from input into snakemake pipeline.
    """
    if len(inputs) != len(outputs):
        raise OSError("//Inputs and outputs are not of the same length: %s <> %s" % (', '.join(inputs), ', '.join(outputs)))
    for infilename, outfilename in zip(inputs, outputs):
        _, fname1 = os.path.split(infilename)
        _process_file(fname1, infilename, outfilename)

def _process_file(fname, inp, outfilename):
    """
    Write the input to the output. Handle raw, zip, or bzip input files.
    """
    print(inp, '=>', outfilename)
    import bz2
    # ungunzip
    if os.path.splitext(fname)[-1] in ['.gz', '.gzip']:
        with open(outfilename, 'wb') as whandle, gzip.open(inp, 'rb') as rhandle:
            shutil.copyfileobj(rhandle, whandle)
    # unbzip2
    elif os.path.splitext(fname)[-1] in ['.bz2', '.bzip2']:
        shell("bzip2 -dc {i} > {o}".format(i=inp, o=outfilename))
    # copy
    else:
        shutil.copy(inp, outfilename)

localrules: prepare_input_data, ALL, prepare_binny

rule ALL:
    input:
        "intermediary.tar.gz",
        "contigs2bin.tsv",
        directory("bins"),
        "contigs2bin_filtered.tsv"

yaml.add_representer(OrderedDict, lambda dumper, data: dumper.represent_mapping('tag:yaml.org,2002:map', data.items()))
yaml.add_representer(tuple, lambda dumper, data: dumper.represent_sequence('tag:yaml.org,2002:seq', data))
yaml.dump(config, open_output('binny.config.yaml'), allow_unicode=True,default_flow_style=False)


rule prepare_input_data:
    input:
        CONTIGS,
        CONTIG_DEPTH if CONTIG_DEPTH else MGaln
    output:
        "intermediary/assembly.fa",
        "intermediary/assembly.contig_depth.txt" if CONTIG_DEPTH else "reads.sorted.bam"
    threads: 1
    resources:
        runtime = "4:00:00",
        mem = MEMCORE,
    message: "Preparing input."
    run:
        prepare_input_files(input, output)

rule format_assembly:
    input:
        "intermediary/assembly.fa"
    output:
        "assembly.fa"
    threads: 1
    resources:
        runtime = "2:00:00",
        mem = MEMCORE,
    message: "Preparing assembly."
    conda: ENVDIR + "/IMP_fasta.yaml"
    shell:
       "fasta_formatter -i {input} -o {output} -w 80"

# contig depth
if not CONTIG_DEPTH:
    rule call_contig_depth:
        input:
            "reads.sorted.bam",
            "assembly.fa"
        output:
            "intermediary/assembly.contig_depth.txt"
        resources:
            runtime = "4:00:00",
            mem = BIGMEMCORE
        threads: 1
        conda: ENVDIR + "/IMP_mapping.yaml"
        log: "logs/analysis_call_contig_depth.log"
        message: "call_contig_depth: Getting data on assembly coverage with mg reads."
        shell:
            """
            echo "Running BEDTools for average depth in each position" >> {log}
            TMP_DEPTH=$(mktemp --tmpdir={TMPDIR} -t "depth_file_XXXXXX.txt")
            genomeCoverageBed -ibam {input[0]} | grep -v "genome" > $TMP_DEPTH
            echo "Depth calculation done" >> {log}

            ## This method of depth calculation was adapted and modified from the CONCOCT code
            perl {SRCDIR}/calcAvgCoverage.pl $TMP_DEPTH {input[1]} > {output}
            echo "Remove the temporary file" >> {log}
            rm $TMP_DEPTH
            """

#gene calling
rule annotate:
    input:
        'assembly.fa'
    output:
        "intermediary/annotation.filt.gff",
        "intermediary/prokka.faa",
        "intermediary/prokka.fna",
        "intermediary/prokka.ffn",
        "intermediary/prokka.fsa",
    threads: 8
    resources:
        runtime = "8:00:00",
        mem = MEMCORE
    log: "logs/analysis_annotate.log"
    conda: ENVDIR + "/IMP_annotation.yaml"
    message: "annotate: Running prokkaC."
    shell:
        """
        export PERL5LIB=$CONDA_PREFIX/lib/site_perl/5.26.2
        export LC_ALL=en_US.utf-8
        if [ ! -f $CONDA_PREFIX/db/hmm/HAMAP.hmm.h3m ]; then
          {BINDIR}/prokkaC --dbdir $CONDA_PREFIX/db --setupdb
        fi
	    {BINDIR}/prokkaC --dbdir $CONDA_PREFIX/db --force --outdir intermediary/ --prefix prokka --noanno --cpus {threads} --metagenome {input[0]} >> {log} 2>&1

	    # Prokka gives a gff file with a long header and with all the contigs at the bottom.  The command below removes the
        # And keeps only the gff table.

        LN=`grep -Hn "^>" intermediary/prokka.gff | head -n1 | cut -f2 -d ":" || if [[ $? -eq 141 ]]; then true; else exit $?; fi`
        LN1=1
        LN=$(($LN-$LN1))
        head -n $LN intermediary/prokka.gff | grep -v "^#" | sort | uniq | grep -v "^==" > {output[0]}
        """

rule cut_rRNA:
    input:
        "assembly.fa",
        "intermediary/annotation.filt.gff"
    output:
        "intermediary/assembly.cut.fa"
    log: "logs/binning_cut_rRNA.log"
    resources:
        runtime = "2:00:00",
        mem = MEMCORE
    threads: 1
    conda: ENVDIR + "/IMP_annotation.yaml"
    message: "cut_rRNA: Cutting contigs for vizbin."
    shell:
        """
        export PERL5LIB=$CONDA_PREFIX/lib/site_perl/5.26.2
        {SRCDIR}/fastaExtractCutRibosomal1000.pl -f {input[0]} -g {input[1]} -l {log} -o {output} -c {config[binning][vizbin][cutoff]}
        """


# essential genes
rule hmmer_essential:
    input:
        "intermediary/prokka.faa",
    output:
        "intermediary/prokka.faa.essential.hmmscan"
    params:
        dbs = DBPATH
    resources:
        runtime = "8:00:00",
        mem = MEMCORE
    conda: ENVDIR + "/IMP_annotation.yaml"
    threads: 12
    log: "logs/analysis_hmmer.essential.log"
    message: "hmmer: Running HMMER for essential."
    shell:
        """
        if [ ! -f {DBPATH}/essential.hmm.h3i ]; then
          hmmpress {DBPATH}/essential.hmm 2>> {log}
        fi
        hmmsearch --cpu {threads} --cut_tc --noali --notextw \
          --tblout {output} {params.dbs}/*.hmm {input} >/dev/null 2>> {log}
        """

rule makegff:
    input:
        "intermediary/prokka.faa.essential.hmmscan",
        "intermediary/annotation.filt.gff",
        "intermediary/prokka.faa"
    output:
        "intermediary/annotation.CDS.RNA.essential.gff"
    resources:
        runtime = "4:00:00",
        mem = MEMCORE
    threads: 1
    conda: ENVDIR + "/IMP_annotation.yaml"
    log: "logs/analysis_makegff.essential.log"
    message: "makegff: Adding hmmer results of essential genes to gff."
    shell:
        """
        export PERL5LIB=$CONDA_PREFIX/lib/site_perl/5.26.2
        {SRCDIR}/hmmscan_addBest2gff.pl -i {input[0]} -a {input[1]} \
         -n essential -o {output} -g $(grep ">" {input[2]} | wc -l) > {log} 2>&1
        """

rule mergegff:
    input:
        "intermediary/annotation.filt.gff",
        "intermediary/annotation.CDS.RNA.essential.gff"
    output:
        "intermediary/annotation_CDS_RNA_hmms.gff"
    resources:
        runtime = "2:00:00",
        mem = MEMCORE
    conda: ENVDIR + "/IMP_annotation.yaml"
    threads: 1
    log: "logs/analysis_mergegff.log"
    message: "mergegff: Merging gffs with hmmer results."
    shell:
        """
        export PERL5LIB=$CONDA_PREFIX/lib/site_perl/5.26.2
        {SRCDIR}/mergegffs.pl {output} {input} > {log} 2>&1
        """

# binning
rule vizbin:
    input:
        "intermediary/assembly.cut.fa"
    output:
        "vizbin.with-contig-names.points"
    resources:
        runtime = "12:00:00",
        mem = BIGMEMCORE
    threads: 1
    conda: ENVDIR + "/IMP_annotation.yaml"
    log: "logs/binning_vizbin.log"
    message: "vizbin: Running VizBin."
    shell:
        """
        TMP_VIZBIN=$(mktemp --tmpdir=intermediary -dt "VIZBIN_XXXXXX")
        java -jar {BINDIR}/VizBin-dist.jar \
         -a {config[binning][vizbin][dimension]} \
         -c {config[binning][vizbin][cutoff]} \
         -i {input} \
         -o $TMP_VIZBIN/data.points \
         -k {config[binning][vizbin][kmer]} \
         -p {config[binning][vizbin][perp]} > {log} 2>&1

         if [ -f $TMP_VIZBIN/data.points ]
           then
             paste <(grep "^>" {input} | sed -e 's/>//') \
              <(cat $TMP_VIZBIN/data.points | sed -e 's/,/\t/') > {output}
           fi
        rm -rf $TMP_VIZBIN
        """

rule prepare_binny:
    input:
       mgdepth='intermediary/assembly.contig_depth.txt',
       vizbin='vizbin.with-contig-names.points' ,
       gff='intermediary/annotation_CDS_RNA_hmms.gff'
    output:
       directory("intermediary/clusterFiles")
    message: "Prepare binny."
    shell:
       """
       mkdir -p {output} || echo "{output} exists"
       """

rule binny:
    input:
       outdir="intermediary/clusterFiles",
       mgdepth='intermediary/assembly.contig_depth.txt',
       vizbin="vizbin.with-contig-names.points",
       gff="intermediary/annotation_CDS_RNA_hmms.gff"
    output:
        expand("intermediary/reachabilityDistanceEstimates.{pk}.{nn}.tsv \
        intermediary/clusterFirstScan.{pk}.{nn}.tsv \
        intermediary/bimodalClusterCutoffs.{pk}.{nn}.tsv \
        intermediary/contigs2clusters.{pk}.{nn}.tsv \
        intermediary/contigs2clusters.{pk}.{nn}.RDS \
        intermediary/finalClusterMap.{pk}.{nn}.png \
        intermediary/clusteringWS.{pk}.{nn}.Rdata".split(),pk=config["binning"]["binny"]["pk"],nn=config["binning"]["binny"]["nn"])
    params:
        plot_functions = SRCDIR + "/IMP_plot_binny_functions.R",
        binnydir="intermediary/"
    resources:
        runtime = "12:00:00",
        mem = BIGMEMCORE
    threads: 1
    conda: ENVDIR + "/IMP_binning.yaml"
    log: "logs/binning_binny.log"
    message: "binny: Running Binny."
    script:
        SRCDIR + "/binny.R"

rule tar_binny_files:
    input:
        expand("intermediary/contigs2clusters.{pk}.{nn}.tsv \
        intermediary/finalClusterMap.{pk}.{nn}.png".split(),pk=config["binning"]["binny"]["pk"],nn=config["binning"]["binny"]["nn"])
    output:
        "intermediary.tar.gz",
        "contigs2bin.tsv",
        "finalClusterMap.png"
    threads: 1
    resources:
        runtime = "8:00:00",
        mem = MEMCORE
    params:
        intermediary = "intermediary/"
    log: "logs/binning_tar_binny_files.log"
    message: "tar_binny_files: Compressing intermediary files from binny."
    shell:
       """
       cp {input[0]} {output[1]}
       cp {input[1]} {output[2]}
       tar cvzf {output[0]} {params.intermediary} >> {log} 2>&1 && rm -r {params.intermediary} >> {log} 2>&1
       """

rule filter_output:
    input:
        "assembly.fa",
        "contigs2bin.tsv",
        OUTPUTDIR,
        config["binning"]["filtering"]["completeness"]),
        config["binning"]["filtering"]["purity"])
    output:
        directory("bins"),
        "contigs2bin_filtered.tsv"
    threads: 1
    resources:
        runtime = "0:30:00",
        mem = MEMCORE
    log: "logs/filter_output.log"
    message: "filtering output."
    shell:
       """
       ./{SRCDIR}/filter_binny_output.py {input[0]} {input[1]} {input[2]} {input[3]} {input[4]}
       """