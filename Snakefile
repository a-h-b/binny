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


# default executable for snakmake
shell.executable("bash")


#functions
def open_output(filename):
    return open(OUTPUTDIR + '/' + filename, 'w+')

def getThreads(max):
    if hasattr(workflow, 'cores'):
        realThreads = max if max <= workflow.cores else workflow.cores
    elif hasattr(workflow, 'nodes'):
        realThreads = max if max <= workflow.nodes else workflow.nodes
    else:
        realThreads = max
    return realThreads

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



# default configuration file
configfile:
    srcdir("config/config.default.yaml")

# some parameters
SRCDIR = srcdir("workflow/scripts")
BINDIR = srcdir("workflow/bin")
ENVDIR = srcdir("workflow/envs")
CONDA_DIR = config['conda_source']

if SRCDIR not in sys.path:
    sys.path.append(SRCDIR)
    import remove_unused_checkm_hmm_profiles as prepCheckM

# get parameters from the config file

# output
if os.path.isabs(os.path.expandvars(config['outputdir'])):
    OUTPUTDIR = os.path.expandvars(config['outputdir'])
else:
    OUTPUTDIR = os.path.join(os.getcwd() , os.path.expandvars(config['outputdir']))

if not os.path.exists(OUTPUTDIR):
            os.makedirs(OUTPUTDIR)

# input
if os.path.isabs(os.path.expandvars(config['raws']['assembly'])):
    CONTIGS = os.path.expandvars(config['raws']['assembly'])
else:
    CONTIGS = os.path.join(os.getcwd(), os.path.expandvars(config['raws']['assembly']))

# Added depth file par to use instead of alignment
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
    # Get filenames of all bam files without extension, even if the name contains '.'
    mappings_ids = ['.'.join(bam.split('/')[-1].split('.')[:-1]) for bam in MGaln]
    # print(mappings_ids)
    # ????:
    # Note that if a rule has multiple output files, Snakemake requires them to all have exactly the same wildcards.
    # Otherwise, it could happen that two jobs running the same rule in parallel want to write to the same file.
    # and
    # The best solution is to have a dictionary that translates a sample id to the inconsistently named files and
    # use a function (see Functions as Input Files) to provide an input file ...
    sample_id_map_dict = {map_id: 'sample_%06.d' % (index + 1) for index, map_id in enumerate(mappings_ids)}

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
if config['mem']['big_mem_avail'] > 0:
    BIGMEMCORE = str(config['mem']['big_mem_per_core_gb']) + "G"
    BIGMEMS = config['mem']['big_mem_avail']
else:
    BIGMEMCORE = False
    

#clean up sample name
SAMPLE = config['sample']
if SAMPLE == "":
    SAMPLE = "_".join(OUTPUTDIR.split("/")[-2:])
SAMPLE = re.sub("_+","_",re.sub("[;|.-]","_",SAMPLE))


#set up DBs, if necessary
if not config['db_path']:
    DBPATH = srcdir('database')
else:
    DBPATH = os.path.expandvars(config['db_path'])
    if not os.path.isabs(DBPATH):
        DBPATH = os.path.join(os.getcwd(), DBPATH)
if not os.path.exists(os.path.join(DBPATH, "taxon_marker_sets_lineage_sorted.tsv")):
    print("Setting up marker database")
    if not os.path.exists(DBPATH):
        os.makedirs(DBPATH)
    if not os.path.exists(os.path.join(DBPATH, "checkm_data_2015_01_16.tar.gz")):	
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

# dump config
yaml.add_representer(OrderedDict, lambda dumper, data: dumper.represent_mapping('tag:yaml.org,2002:map', data.items()))
yaml.add_representer(tuple, lambda dumper, data: dumper.represent_sequence('tag:yaml.org,2002:seq', data))
yaml.dump(config, open_output('binny.config.yaml'), allow_unicode=True, default_flow_style=False)

workdir:
    OUTPUTDIR

onsuccess:
    shell("mkdir -p job.errs.outs &>> logs/cleanup.log; ( mv binny*stdout job.errs.outs || touch job.errs.outs ) &>> logs/cleanup.log") 


# Snakemake workflow

localrules: prepare_input_data, ALL


rule ALL:
    input:
        os.path.join(OUTPUTDIR, 'binny.done')


rule prepare_input_data:
    input:
        CONTIGS,
        CONTIG_DEPTH if CONTIG_DEPTH else MGaln
    output:
        os.path.join(OUTPUTDIR, "intermediary/assembly.fa"),
        os.path.join(OUTPUTDIR, "intermediary/assembly.contig_depth.txt") if CONTIG_DEPTH
        else [os.path.join(OUTPUTDIR, "intermediary/reads_{0}_sorted.bam".format(sample_id_map_dict[mappings_id]))
              for mappings_id in mappings_ids]
    threads: 1 
#    resources:
#        runtime = "4:00:00",
#        mem = MEMCORE
    message:
        "Preparing input."
    run:
        prepare_input_files(input, output)

rule format_assembly:
    input:
        os.path.join(OUTPUTDIR, "intermediary/assembly.fa")
    output:
        os.path.join(OUTPUTDIR, "intermediary/assembly.formatted.fa")
    params:
        min_length=int(config["min_cont_length_cutoff_marker"]) - 1
    threads: 
        getThreads(1)
    resources:
        runtime = "2:00:00",
        mem = MEMCORE
    message:
        "Preparing assembly."
    conda:
       os.path.join(ENVDIR, "fasta_processing.yaml")
    shell:
       """
       seqkit seq {input} -o {output} -w 80 -m 499 \
       && \
       rm -f {input}
       """

# contig depth
if not CONTIG_DEPTH:
    rule call_contig_depth:
        input:
            assembly=os.path.join(OUTPUTDIR, "intermediary/assembly.formatted.fa"),
            mapping=lambda wildcards: os.path.join(OUTPUTDIR,
                "intermediary/reads_{0}_sorted.bam".format(sample_id_map_dict[wildcards.sample]))
        output:
            os.path.join(OUTPUTDIR, "intermediary/assembly_contig_depth_{sample}.txt")
        resources:
            runtime = "4:00:00",
            mem = BIGMEMCORE if BIGMEMCORE else MEMCORE
        threads:
            1
        conda:
            os.path.join(ENVDIR, "mapping.yaml")
        log:
            os.path.join(OUTPUTDIR, "logs/analysis_call_contig_depth_{sample}.log")
        message:
            os.path.join(OUTPUTDIR, "call_contig_depth: Getting data on assembly coverage with mg reads.")
        shell:
            """
            echo "Running BEDTools for average depth in each position" >> {log}
            export TMPDIR={TMPDIR}
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
            [os.path.join(OUTPUTDIR, f"intermediary/assembly_contig_depth_{mappings_id}.txt")
             for mappings_id in mappings_ids]
        output:
            os.path.join(OUTPUTDIR, "intermediary/assembly.contig_depth.txt")
        params:
            int_dir=os.path.join(OUTPUTDIR, "intermediary")
        resources:
            runtime = "1:00:00",
            mem = MEMCORE
        threads:
            getThreads(1)
        conda:
            os.path.join(ENVDIR, "mapping.yaml")
        log:
            os.path.join(OUTPUTDIR, "logs/merge_contig_depth.log")
        message:
            "Merging depth files."
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
                export TMPDIR={TMPDIR}
                TMP_DEPTH=$(mktemp --tmpdir={TMPDIR} "tmp_XXXXXXXXXXX.tsv")
                paste {output} <( cut -f 2 $file) > $TMP_DEPTH \
                      && mv $TMP_DEPTH {output}
              fi
            done
            rm {params.int_dir}/assembly_contig_depth_*.txt
            """

#gene calling
rule annotate:
    input:
        assembly=os.path.join(OUTPUTDIR, 'intermediary/assembly.formatted.fa')
    output:
        gff=os.path.join(OUTPUTDIR, "intermediary/prokka.gff"),
        gff_filt=os.path.join(OUTPUTDIR, "intermediary/annotation.filt.gff"),
        faa=os.path.join(OUTPUTDIR, "intermediary/prokka.faa"),
    params:
        int_dir=os.path.join(OUTPUTDIR, "intermediary")
    threads:
        getThreads(24)
    resources:
        runtime = "48:00:00",
        mem = MEMCORE
    log:
        os.path.join(OUTPUTDIR, "logs/analysis_annotate.log")
    benchmark:
        os.path.join(OUTPUTDIR, "logs/analysis_annotate_benchmark.txt")
    conda:
        PROKKA_ENV if PROKKA_ENV else os.path.join(ENVDIR, "prokka.yaml")
    message:
        "annotate: Running prokkaP."
    shell:
        """
        export PERL5LIB=$CONDA_PREFIX/lib/site_perl/5.26.2
        export LC_ALL=en_US.utf-8
        export TMPDIR={TMPDIR}
	    {BINDIR}/prokkaP --dbdir $CONDA_PREFIX/db --force --outdir {params.int_dir} --tmpdir {TMPDIR} --prefix prokka \
	                     --noanno --cpus {threads} --metagenome {input.assembly} >> {log} 2>&1
        
	    # Prokka gives a gff file with a long header and with all the contigs at the bottom.  The command below keeps
	    # only the gff table.

        LN=`grep -Hn "^>" {output.gff} | head -n1 | cut -f2 -d ":" || if [[ $? -eq 141 ]]; then true; else exit $?; fi`
        LN1=1
        LN=$(($LN-$LN1))
        head -n $LN {output.gff} | grep -v "^#" | sort | uniq | grep -v "^==" > {output.gff_filt}
        """

# Find markers on contigs
rule mantis_checkm_marker_sets:
    input:
        proteins=os.path.join(OUTPUTDIR, "intermediary/prokka.faa")
    output:
        os.path.join(OUTPUTDIR, "intermediary/mantis_out/output_annotation.tsv"),
        os.path.join(OUTPUTDIR, "intermediary/mantis_out/integrated_annotation.tsv"),
        os.path.join(OUTPUTDIR, "intermediary/mantis_out/consensus_annotation.tsv"),
        out_dir=directory(os.path.join(OUTPUTDIR, "intermediary/mantis_out"))
    params:
        binny_cfg=srcdir("config/binny_mantis.cfg")
    resources:
        runtime = "48:00:00",
        mem = MEMCORE
    conda:
        MANTIS_ENV if MANTIS_ENV else os.path.join(ENVDIR, "mantis.yaml")
    threads:
        getThreads(80)
    log:
        os.path.join(OUTPUTDIR, "logs/analysis_checkm_markers.log")
    benchmark:
        os.path.join(OUTPUTDIR, "logs/analysis_checkm_markers_benchmark.txt")
    message:
        "MANTIS: Running MANTIS with CheckM marker sets."
    shell:
        """
        if [ -d {output.out_dir} ]; then rm {output.out_dir}/* || true ; fi >> {log} 2>&1
        mantis run -i {input.proteins} \
                   -da heuristic \
                   -mc {params.binny_cfg} \
                   -o {output.out_dir} \
                   -c {threads} \
                   --no_taxonomy \
                   -et 1e-3 >> {log} 2>&1
        """

rule binny:
    input:
        mgdepth=os.path.join(OUTPUTDIR, 'intermediary/assembly.contig_depth.txt'),
        raw_gff=os.path.join(OUTPUTDIR, 'intermediary/annotation.filt.gff'),
        assembly=os.path.join(OUTPUTDIR, "intermediary/assembly.formatted.fa"),
        hmm_markers=os.path.join(OUTPUTDIR, "intermediary/mantis_out/consensus_annotation.tsv")
    output:
        os.path.join(OUTPUTDIR, 'binny.done')
    params:
        binny_out=OUTPUTDIR,
        sample=SAMPLE,
        py_functions=SRCDIR + "/binny_functions.py",
        t2p=DBPATH + "/pfam/tigrfam2pfam.tsv",
        marker_sets=DBPATH + "/taxon_marker_sets_lineage_sorted.tsv",
        gff=os.path.join(OUTPUTDIR, "intermediary/annotation_CDS_RNA_hmms_checkm.gff"),
        min_completeness=config["bin_quality"]["min_completeness"],
        start_completeness=config["bin_quality"]["start_completeness"],
        purity=config["bin_quality"]["purity"],
        kmers=config["kmers"],
        mask_disruptive_sequences=config["mask_disruptive_sequences"],
        extract_scmags=config["extract_scmags"],
        coassembly_mode=config["coassembly_mode"],
        min_cutoff=config["min_cont_length_cutoff"],
        max_cutoff=config["max_cont_length_cutoff"],
        min_cutoff_marker=config["min_cont_length_cutoff_marker"],
        max_cutoff_marker=config["max_cont_length_cutoff_marker"],
        nx_val=config["NX_value"],
        max_n_contigs=config["max_n_contigs"],
        max_marker_lineage_depth_lvl=config["max_marker_lineage_depth_lvl"],
        distance_metric=config["distance_metric"],
        max_embedding_tries=config["embedding"]["max_iterations"],
        include_depth_initial=config["clustering"]["include_depth_initial"],
        include_depth_main=config["clustering"]["include_depth_main"],
        hdbscan_min_samples_range=config["clustering"]["hdbscan_min_samples_range"],
        hdbscan_epsilon_range=config["clustering"]["hdbscan_epsilon_range"],
        write_contig_data=config["write_contig_data"]
    resources:
        runtime = "12:00:00",
        mem = BIGMEMCORE if BIGMEMCORE else MEMCORE
    threads:
        getThreads(BIGMEMS) if BIGMEMCORE else getThreads(80)
    conda:
        os.path.join(ENVDIR, "binny_linux.yaml")
    log:
        os.path.join(OUTPUTDIR, "logs/binning_binny.log")
    benchmark:
        os.path.join(OUTPUTDIR, "logs/binning_binny_benchmark.txt")
    message:
            "binny: Running Python Binny."
    script:
        os.path.join(SRCDIR, "binny_main.py")
