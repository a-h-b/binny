import gzip
import os
import re
import shutil
import sys
import tarfile
import urllib.request

import pandas as pd
import yaml


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

if SRCDIR not in sys.path:
    sys.path.append(SRCDIR)
    import remove_unused_checkm_hmm_profiles as prepCheckM

# get parameters from the config file

# output
if os.path.isabs(os.path.expandvars(config['outputdir'])):
    OUTPUTDIR = os.path.expandvars(config['outputdir'])
else:
    OUTPUTDIR = os.path.join(os.getcwd() , os.path.expandvars(config['outputdir']))

# input
if os.path.isabs(os.path.expandvars(config['raws']['assembly'])):
    CONTIGS = os.path.expandvars(config['raws']['assembly'])
else:
    CONTIGS = os.path.join(os.getcwd(), os.path.expandvars(config['raws']['assembly']))
    
# Added depth file par to us instead of alignment
if config['raws']['contig_depth']:
    if os.path.isabs(os.path.expandvars(config['raws']['contig_depth'])):
        CONTIG_DEPTH = os.path.expandvars(config['raws']['contig_depth'])
    else:
        CONTIG_DEPTH = os.path.join(os.getcwd(), os.path.expandvars(config['raws']['contig_depth']))
else:
    CONTIG_DEPTH = None
    if os.path.isabs(os.path.expandvars(config['raws']['metagenomics_alignment'])):
        MGaln = os.path.expandvars(config['raws']['metagenomics_alignment'])
    else:
        MGaln = os.path.join(os.getcwd(), os.path.expandvars(config['raws']['metagenomics_alignment']))

# hardware parameters
MEMCORE = str(config['mem']['normal_mem_per_core_gb']) + "G"
if config['mem']['big_mem_avail']:
    BIGMEMCORE = str(config['mem']['big_mem_per_core_gb']) + "G"
else:
    BIGMEMCORE = False

def getThreads(max):
    if workflow.cores:
        realThreads = max if max <= workflow.cores else workflow.cores
    elif workflow.nodes:
        realThreads = max if max <= workflow.nodes else workflow.nodes
    else:
        realThreads = max
    return realThreads


SAMPLE = config['sample']
if SAMPLE == "":
    SAMPLE = "_".join(OUTPUTDIR.split("/")[-2:])
SAMPLE = re.sub("_+","_",re.sub("[;|.-]","_",SAMPLE))
DBPATH = os.path.expandvars(config['db_path'])
if not os.path.isabs(DBPATH):
    DBPATH = os.path.join(os.getcwd(), DBPATH)
if not os.path.exists(DBPATH):
    print("Setting up marker database")
    os.makedirs(DBPATH)
    urllib.request.urlretrieve("https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz", DBPATH + "/checkm_data_2015_01_16.tar.gz")
    print(os.path.join( DBPATH, "/checkm_data_2015_01_16.tar.gz"))
    checkm_tar = tarfile.open(os.path.join( DBPATH, "checkm_data_2015_01_16.tar.gz"))
    checkm_tar.extract("./taxon_marker_sets.tsv",DBPATH)
    checkm_tar.extract("./pfam/tigrfam2pfam.tsv",DBPATH)
    checkm_tar.extract("./hmms/checkm.hmm",DBPATH)
    markers_df = pd.read_csv(os.path.join(DBPATH, 'taxon_marker_sets.tsv'), sep='\t', skipinitialspace=True, header=None)
    markers_df = markers_df.sort_values(markers_df.columns[2])
    markers_df.to_csv(os.path.join(DBPATH, "taxon_marker_sets_lineage_sorted.tsv"), header=None, index=None, sep="\t")
    prepCheckM.remove_unused_checkm_hmm_profiles(os.path.join(DBPATH, "hmms/checkm.hmm"), os.path.join(DBPATH, 'taxon_marker_sets.tsv'), os.path.join(DBPATH, "pfam/tigrfam2pfam.tsv"), os.path.join(DBPATH, "hmms"))
    if os.path.exists(os.path.join(DBPATH, "checkm_data_2015_01_16.tar.gz")):
        os.remove(os.path.join(DBPATH, "checkm_data_2015_01_16.tar.gz"))
    if os.path.exists(os.path.join(DBPATH, "hmms/checkm.hmm")):
        os.remove(os.path.join(DBPATH, "hmms/checkm.hmm"))
    if os.path.exists(os.path.join(DBPATH, "taxon_marker_sets.tsv")) and os.path.exists(os.path.join(DBPATH, "taxon_marker_sets_lineage_sorted.tsv")):
        os.remove(os.path.join(DBPATH, "taxon_marker_sets.tsv"))

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

localrules: prepare_input_data, ALL


rule ALL:
    input:
        "bins/"

yaml.add_representer(OrderedDict, lambda dumper, data: dumper.represent_mapping('tag:yaml.org,2002:map', data.items()))
yaml.add_representer(tuple, lambda dumper, data: dumper.represent_sequence('tag:yaml.org,2002:seq', data))
yaml.dump(config, open_output('binny.config.yaml'), allow_unicode=True,default_flow_style=False)

rule prepare_input_data:
    input:
        CONTIGS,
        CONTIG_DEPTH if CONTIG_DEPTH else MGaln
    output:
        "intermediary/assembly.fa",
        "intermediary/assembly.contig_depth.txt" if CONTIG_DEPTH else "intermediary/reads.sorted.bam"
    threads: 1
    resources:
        runtime = "4:00:00",
        mem = MEMCORE
    message: "Preparing input."
    run:
        prepare_input_files(input, output)

rule format_assembly:
    input:
        "intermediary/assembly.fa"
    output:
        "intermediary/assembly.formatted.fa"
    threads: 1
    resources:
        runtime = "2:00:00",
        mem = MEMCORE
    message: "Preparing assembly."
    conda: 
       os.path.join(ENVDIR, "IMP_fasta_no_v.yaml")
    shell:
       "fasta_formatter -i {input} -o {output} -w 80"

# contig depth
if not CONTIG_DEPTH:
    rule call_contig_depth:
        input:
            "intermediary/reads.sorted.bam",
            "intermediary/assembly.formatted.fa"
        output:
            "intermediary/assembly.contig_depth.txt"
        resources:
            runtime = "4:00:00",
            mem = BIGMEMCORE if BIGMEMCORE else MEMCORE
        threads: 
            getThreads(2) if BIGMEMCORE else getThreads(8)
        conda: 
            os.path.join(ENVDIR, "IMP_mapping.yaml")
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
        'intermediary/assembly.formatted.fa'
    output:
        "intermediary/annotation.filt.gff",
        "intermediary/prokka.faa",
        "intermediary/prokka.fna",
        "intermediary/prokka.ffn",
        "intermediary/prokka.fsa",
    threads: 
        # getThreads(20)
        workflow.cores
    resources:
        runtime = "120:00:00",
        mem = MEMCORE
    log: "logs/analysis_annotate.log"
    benchmark: "logs/analysis_annotate_benchmark.txt"
    conda: 
        os.path.join(ENVDIR, "IMP_annotation.yaml")
    message: "annotate: Running prokkaP."
    shell:
        """
        export PERL5LIB=$CONDA_PREFIX/lib/site_perl/5.26.2
        export LC_ALL=en_US.utf-8
	{BINDIR}/prokkaP --dbdir $CONDA_PREFIX/db --force --outdir intermediary/ --tmpdir {TMPDIR} --prefix prokka --noanno --cpus {threads} --metagenome {input[0]} >> {log} 2>&1
        
	# Prokka gives a gff file with a long header and with all the contigs at the bottom.  The command below keeps only the gff table.

        LN=`grep -Hn "^>" intermediary/prokka.gff | head -n1 | cut -f2 -d ":" || if [[ $? -eq 141 ]]; then true; else exit $?; fi`
        LN1=1
        LN=$(($LN-$LN1))
        head -n $LN intermediary/prokka.gff | grep -v "^#" | sort | uniq | grep -v "^==" > {output[0]}
        """

# Find markers on contigs
rule mantis_checkm_marker_sets:
    input:
        "intermediary/prokka.faa"
    output:
        "intermediary/mantis_out/output_annotation.tsv",
        "intermediary/mantis_out/integrated_annotation.tsv",
        "intermediary/mantis_out/consensus_annotation.tsv"
    resources:
        runtime = "48:00:00",
        mem = MEMCORE
    conda: 
        os.path.join(BINDIR, "mantis/mantis_env.yml")
    threads: 
        # getThreads(20)
        # getThreads(14)
        workflow.cores
    log: "logs/analysis_checkm_markers.log"
    benchmark: "logs/analysis_checkm_markers_benchmark.txt"
    message: "MANTIS: Running MANTIS with CheckM marker sets."
    shell:
        """
        if [ -d intermediary/mantis_out ]; then rm intermediary/mantis_out/* || true ; fi >> {log} 2>&1
        python {BINDIR}/mantis/ run_mantis -i {input[0]} -da heuristic -mc {BINDIR}/mantis/MANTIS.config \
                                           -o intermediary/mantis_out -c {threads} -et 1e-3 >> {log} 2>&1
        """

rule binny:
    input:
        mgdepth='intermediary/assembly.contig_depth.txt',
        raw_gff='intermediary/annotation.filt.gff',
        assembly="intermediary/assembly.formatted.fa",
        hmm_markers="intermediary/mantis_out/consensus_annotation.tsv"
    output:
        directory("bins")
    params:
        sample=SAMPLE,
        py_functions = SRCDIR + "/binny_functions.py",
        binnydir="intermediary/",
        t2p=DBPATH + "/pfam/tigrfam2pfam.tsv",
        marker_sets=DBPATH + "/taxon_marker_sets_lineage_sorted.tsv",
        completeness=config["binning"]["binny"]["bin_quality"]["completeness"],
        purity=config["binning"]["binny"]["bin_quality"]["purity"],
        kmers=config["binning"]["binny"]["kmers"],
        cutoff=config["binning"]["binny"]["cutoff"],
        cutoff_marker=config["binning"]["binny"]["cutoff_marker"],
        max_n_contigs=config["binning"]["binny"]["max_n_contigs"],
        distance_metric=config["binning"]["binny"]["distance_metric"],
        max_embedding_tries=config["binning"]["binny"]["embedding"]["max_iterations"],
        tsne_early_exag_iterations=config["binning"]["binny"]["embedding"]["tsne_early_exag_iterations"],
        tsne_main_iterations=config["binning"]["binny"]["embedding"]["tsne_main_iterations"],
        include_depth_initial=config["binning"]["binny"]["clustering"]["include_depth_initial"],
        include_depth_main=config["binning"]["binny"]["clustering"]["include_depth_main"],
        hdbscan_min_samples=config["binning"]["binny"]["clustering"]["hdbscan_min_samples"],
        hdbscan_epsilon=config["binning"]["binny"]["clustering"]["hdbscan_epsilon"],
        gff="intermediary/annotation_CDS_RNA_hmms_checkm.gff",
    resources:
        runtime = "12:00:00",
        mem = BIGMEMCORE if BIGMEMCORE else MEMCORE
    threads: 
        # getThreads(2) if BIGMEMCORE else getThreads(20)
        workflow.cores
    conda: 
        os.path.join(ENVDIR, "py_binny_linux.yaml")
    log: "logs/binning_binny.log"
    benchmark: "logs/binning_binny_benchmark.txt"
    message: "binny: Running Python Binny."
    script:
        os.path.join(SRCDIR, "binny_main.py")

