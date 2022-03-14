import gzip
import os
import glob
import re
import shutil
import sys
import tarfile
import urllib.request

import pandas as pd
import yaml


def open_output(filename):
    return open(OUTPUTDIR+'/'+filename, 'w+')

# default executable for snakmake
shell.executable("bash")

# default configuration file
configfile:
    srcdir("config/config.default.yaml")

# some parameters
SRCDIR = srcdir("workflow/scripts")
BINDIR = srcdir("workflow/bin")
ENVDIR = srcdir("workflow/envs")
CONDA_DIR = srcdir("conda")

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
        CONTIG_DEPTH = os.path.join(os.getcwd(),os.path.expandvars(config['raws']['contig_depth']))
else:
    CONTIG_DEPTH = None
    if all([os.path.isabs(path) for path in glob.glob(config['raws']['metagenomics_alignment'])]):
        MGaln = [os.path.expandvars(path) for path in glob.glob(config['raws']['metagenomics_alignment'])]
    else:
        MGaln = [os.path.join(os.getcwd(), os.path.expandvars(path)) for path in glob.glob(config['raws']['metagenomics_alignment'])]
    # print(MGaln)
    # print(config['raws']['metagenomics_alignment'])
    # print(glob.glob(config['raws']['metagenomics_alignment']))
    # Get filenames of all bam files without extension, even if the name contains '.'
    mappings_ids = ['.'.join(bam.split('/')[-1].split('.')[:-1]) for bam in MGaln]
    # print(mappings_ids)
    # ????:
    # Note that if a rule has multiple output files, Snakemake requires them to all have exactly the same wildcards.
    # Otherwise, it could happen that two jobs running the same rule in parallel want to write to the same file.
    # and
    # The best solution is to have a dictionary that translates a sample id to the inconsistently named files and
    # use a function (see Functions as Input Files) to provide an input file ...
    garbage_dict_so_snakemake_gets_it = {map_id: 'sample_%06.d' % (index + 1) for index, map_id in enumerate(mappings_ids)}
    # print(garbage_dict_so_snakemake_gets_it)

# Use existing env for Prokka if specified
if config['prokka_env'] and config['prokka_env'].split('.')[-1] in ['yaml', 'yml']:
    if os.path.isabs(config['prokka_env']):
        PROKKA_ENV = os.path.expandvars(config['prokka_env'])
    else:
        PROKKA_ENV = os.path.join(os.getcwd(), config['prokka_env'])
    print(PROKKA_ENV)
elif config['prokka_env']:
    PROKKA_ENV = config['prokka_env']
    print('named', PROKKA_ENV)
else:
    PROKKA_ENV = None
# Use existing env for Mantis if specified
if config['mantis_env'] and config['mantis_env'].split('.')[-1] in ['yaml', 'yml']:
    if os.path.isabs(config['mantis_env']):
        MANTIS_ENV = os.path.expandvars(config['mantis_env'])
    else:
        MANTIS_ENV = os.path.join(os.getcwd(), config['mantis_env'])
    print(MANTIS_ENV)
elif config['mantis_env']:
    MANTIS_ENV = config['mantis_env']
    print('named', MANTIS_ENV)
else:
    MANTIS_ENV = None

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
    urllib.request.urlretrieve("https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz", os.path.join(DBPATH, "checkm_data_2015_01_16.tar.gz"))
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
    print("Initializing conda environments.")

# temporary directory will be stored inside the OUTPUTDIR directory
# unless an absolute path is set
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
        # "intermediary/assembly.contig_depth.txt" if CONTIG_DEPTH else expand("intermediary/reads_{mappings_id}_sorted.bam", mappings_id=mappings_ids)
        "intermediary/assembly.contig_depth.txt" if CONTIG_DEPTH else ["intermediary/reads_{0}_sorted.bam".format(garbage_dict_so_snakemake_gets_it[mappings_id]) for mappings_id in mappings_ids]
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
    threads:
        workflow.cores
    resources:
        runtime = "2:00:00",
        mem = MEMCORE
    message: "Preparing assembly."
    conda:
       os.path.join(ENVDIR, "fasta_processing.yaml")
    shell:
       "seqkit seq {input} -o {output} -w 80"

# print([f"intermediary/reads_{mappings_id}_sorted.bam" for mappings_id in mappings_ids])
# print([f"intermediary/assembly_contig_depth_{mappings_id}.txt" for mappings_id in mappings_ids])
# print(expand(["intermediary/reads_{mappings_id}_sorted.bam", "intermediary/assembly.formatted.fa"], mappings_id=mappings_ids))

# contig depth
if not CONTIG_DEPTH:
    rule call_contig_depth:
        input:
            # "intermediary/reads_{sample}_sorted.bam"
            # [f"intermediary/reads_{mappings_id}_sorted.bam" for mappings_id in mappings_ids]
            # expand(["intermediary/reads_{mappings_id}_sorted.bam"], mappings_id=mappings_ids)
            assembly = "intermediary/assembly.formatted.fa",
            mapping=lambda wildcards: "intermediary/reads_{0}_sorted.bam".format(garbage_dict_so_snakemake_gets_it[wildcards.sample])
        output:
            "intermediary/assembly_contig_depth_{sample}.txt"
            # [f"intermediary/assembly_contig_depth_{mappings_id}.txt" for mappings_id in mappings_ids]
            # expand("intermediary/assembly_contig_depth_{mappings_id}.txt", mappings_id=mappings_ids)
        resources:
            runtime = "4:00:00",
            mem = BIGMEMCORE if BIGMEMCORE else MEMCORE
        threads:
            # getThreads(2) if BIGMEMCORE else getThreads(8)
            max(1, int(workflow.cores / len(mappings_ids))) if mappings_ids else 1
        conda:
            os.path.join(ENVDIR, "mapping.yaml")
        log: "logs/analysis_call_contig_depth_{sample}.log"
        message: "call_contig_depth: Getting data on assembly coverage with mg reads."
        shell:
            """
            echo "Running BEDTools for average depth in each position" >> {log}
            TMP_DEPTH=$(mktemp --tmpdir={TMPDIR} "depth_file_XXXXXXXXXXXXXXXX.txt")
            genomeCoverageBed -ibam {input.mapping} | grep -v "genome" > $TMP_DEPTH
            echo "Depth calculation done" >> {log}

            ## This method of depth calculation was adapted and modified from the CONCOCT code
            echo "Getting average contig depth." >> {log}
            perl {SRCDIR}/calcAvgCoverage.pl $TMP_DEPTH {input.assembly} > {output} && \
            echo "Done. Removing the temporary file" >> {log} 
            rm $TMP_DEPTH
            """

    rule merge_contig_depths:
        input:
            [f"intermediary/assembly_contig_depth_{mappings_id}.txt" for mappings_id in mappings_ids]
        output:
            "intermediary/assembly.contig_depth.txt"
        resources:
            runtime = "1:00:00",
            mem = MEMCORE
        threads:
            getThreads(1)
        conda:
            os.path.join(ENVDIR, "mapping.yaml")
        log: "logs/merge_contig_depth.log"
        message: "Merging depth files."
        shell:
            """
            first_file=true
            for file in {input}; do
              if [[ $first_file == 'true' ]]; then
                # echo "First file."
                cp $file {output}
                first_file=false
              else
                # echo "File $COUNTER"
                TMP_DEPTH=$(mktemp --tmpdir={TMPDIR} "tmp_XXXXXXXXXXX.tsv")
                paste {output} <( cut -f 2 $file) > $TMP_DEPTH \
                      && mv $TMP_DEPTH {output}
              fi
            done
            rm intermediary/assembly_contig_depth_*.txt
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
        "intermediary/prokka.fsa"
    threads:
        # getThreads(20)
        workflow.cores
    resources:
        runtime = "120:00:00",
        mem = MEMCORE
    log: "logs/analysis_annotate.log"
    benchmark: "logs/analysis_annotate_benchmark.txt"
    conda:
        PROKKA_ENV if PROKKA_ENV else os.path.join(ENVDIR, "prokka.yaml")
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
    params:
        binny_cfg=srcdir("config/binny_mantis.cfg")
    resources:
        runtime = "48:00:00",
        mem = MEMCORE
    conda:
        MANTIS_ENV if MANTIS_ENV else os.path.join(ENVDIR, "mantis.yaml")
    threads:
        workflow.cores
    log: "logs/analysis_checkm_markers.log"
    benchmark: "logs/analysis_checkm_markers_benchmark.txt"
    message: "MANTIS: Running MANTIS with CheckM marker sets."
    shell:
        """
        if [ -d intermediary/mantis_out ]; then rm intermediary/mantis_out/* || true ; fi >> {log} 2>&1
        mantis run -i {input[0]} \
                   -da heuristic \
                   -mc {params.binny_cfg} \
                   -o intermediary/mantis_out \
                   -c {threads} \
                   --no_taxonomy \
                   -et 1e-3 >> {log} 2>&1
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
        gff="intermediary/annotation_CDS_RNA_hmms_checkm.gff",
        min_completeness=config["binning"]["binny"]["bin_quality"]["min_completeness"],
        start_completeness=config["binning"]["binny"]["bin_quality"]["start_completeness"],
        purity=config["binning"]["binny"]["bin_quality"]["purity"],
        kmers=config["binning"]["binny"]["kmers"],
        min_cutoff=config["binning"]["binny"]["min_cont_length_cutoff"],
        max_cutoff=config["binning"]["binny"]["max_cont_length_cutoff"],
        min_cutoff_marker=config["binning"]["binny"]["min_cont_length_cutoff_marker"],
        max_cutoff_marker=config["binning"]["binny"]["max_cont_length_cutoff_marker"],
        nx_val=config["binning"]["binny"]["NX_value"],
        max_n_contigs=config["binning"]["binny"]["max_n_contigs"],
        distance_metric=config["binning"]["binny"]["distance_metric"],
        max_embedding_tries=config["binning"]["binny"]["embedding"]["max_iterations"],
        include_depth_initial=config["binning"]["binny"]["clustering"]["include_depth_initial"],
        include_depth_main=config["binning"]["binny"]["clustering"]["include_depth_main"],
        hdbscan_min_samples_range=config["binning"]["binny"]["clustering"]["hdbscan_min_samples_range"],
        hdbscan_epsilon_range=config["binning"]["binny"]["clustering"]["hdbscan_epsilon_range"]
    resources:
        runtime = "12:00:00",
        mem = BIGMEMCORE if BIGMEMCORE else MEMCORE
    threads:
        # getThreads(2) if BIGMEMCORE else getThreads(20)
        workflow.cores
    conda:
        os.path.join(ENVDIR, "binny_linux.yaml")
    log: "logs/binning_binny.log"
    benchmark: "logs/binning_binny_benchmark.txt"
    message: "binny: Running Python Binny."
    script:
        os.path.join(SRCDIR, "binny_main.py")
