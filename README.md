# binny

## Installing binny
For binny to work, you need [conda](https://www.anaconda.com/). 

1) Clone this repository to your disk:
```
git clone https://github.com/a-h-b/binny.git
```
Change into the binny directory:
```
cd binny
```
At this point, you have all the scripts you need to run the workflow using snakemake. If you want to use the **comfortable binny wrapper**, follow the points 2-6. 

2) Adjust the file VARIABLE_CONFIG to your requirements (have a tab between the variable name and your setting):
* SNAKEMAKE_VIA_CONDA - set this to true, if you don't have snakemake in your path and want to install it via conda. Leave empty, if you don't need an additional snakemake.
* SNAKEMAKE_EXTRA_ARGUMENTS - if you want to pass additional arguments to snakemake, put them here (e.g. --latency-wait=320 for slower file systems). Leave empty usually.
* LOADING_MODULES - insert a bash command to load modules, if you need them to run conda. Leave empty, if you don't need to load a module.
* SUBMIT_COMMAND - insert the bash command you'll usually use to submit a job to your cluster to run on a single cpu for a few days. You only need this, if you want to have the snakemake top instance running in a submitted job. You alternatively have the option to run it on the frontend via tmux. Leave empty, if you want to use this version and have [tmux](https://github.com/tmux/tmux/wiki) installed.
* BIND_JOBS_TO_MAIN - if you use the option to run the snakemake top instance in a submitted job and need to bind the other jobs to the same node, you can set this option to true. See FAQ below for more details.
* NODENAME_VAR - if you use the BIND_JOBS_TO_MAIN option, you need to let dadasnake know, how to access the node name (e.g.SLURMD_NODENAME on slurm).
* SCHEDULER - insert the name of the scheduler you want to use (currently `slurm`, `slurm_simple` or `sge`). This determines the cluster config given to snakemake, e.g. the cluster config file for slurm is config/slurm.config.yaml . Also check that the settings in this file is correct. If you have a different system, contact us ( https://github.com/a-h-b/binny/issues ).
* MAX_THREADS - set this to the maximum number of cores you want to be using in a run. If you don't set this, the default will be 50 (which is more than will be used). Users can override this setting at runtime.
* NORMAL_MEM_EACH - set the size of the RAM of one core of your normal copute nodes (e.g. 8G). If you're not planning to use binny to submit to a cluster, you don't need to set this. 
* BIGMEM_MEM_EACH - set the size of the RAM of one core of your bigmem (or highmem) compute nodes. If you're not planning to use binny to submit to a cluster or don't have separate bigmem nodes, you don't need to set this.


3) Decide how you want to run binny, if you let it submit jobs to the cluster:
Only do one of the two:
* if you want to submit the process running snakemake to the cluster:
```
cp runscripts/binny_submit.sh binny
chmod 755 binny
```
* if you want to keep the process running snakemake on the frontend using tmux:
```
cp runscripts/binny_tmux.sh binny
chmod 755 binny
```

4) **optional**: Install snakemake via conda:
If you want to use snakemake via conda (and you've set SNAKEMAKE_VIA_CONDA to true), install the environment, as [recommended by Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html):
```
conda install -c conda-forge mamba
if ! [ -x "$(command -v unzip)" ]; then
  mamba create --prefix $PWD/conda/snakemake_env snakemake mamba unzip -c conda-forge -c bioconda
else
  mamba create --prefix $PWD/conda/snakemake_env snakemake mamba -c conda-forge -c bioconda
fi
```
*Note*: binny will not work with snakemake version 6.3.0 .

5) **optional**: Set permissions / PATH:
binny is meant to be used by multiple users. Set the permissions accordingly. I'd suggest:
* to have read access for all files for the users plus 
* execution rights for the binny file and the .sh scripts in the subfolder runscripts
* read, write and execution rights for the conda subfolder 
* to add the binny directory to your path. 
* It can also be useful to make the VARIABLE_CONFIG file not-writable, because you will always need it. The same goes for config.default.yaml once you've set the paths to the databases you want to use (see below).

6) Initialize conda environments:
This run sets up the conda environments, including [Mantis](https://academic.oup.com/gigascience/article/10/6/giab042/6291114), that will be usable by all users: 
```
./binny -i config/config.init.yaml 
```
This step will take several minutes to an hour. It will also create a folder with the name "database". It contains the database of unique genes from [CheckM](https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz). You can move this elsewhere and specify the path in the config, if you wish.

7) **Optional** test run:
You should be able to test run by 
```
./binny -l -n "TESTRUN" -r config/config.test.yaml
```
If all goes well, binny will run in the current session, load the conda environments, and make and fill a directory called testoutput. A completed run contains a file "contigs2bin.tsv". 
If you don't want to see binny's guts at this point, you can also run this with the -c or -f settings to submit to your cluster or start a tmux session (see How to run binny below). 

## How to run binny
To run the binny, you need a config file, plus data: 
* The config file (in yaml format) is read by Snakemake to determine the inputs, arguments and outputs. 
* You need contigs in a fasta file and the alignments of metagenomic reads in bam format, or the contigs plus a tab-separated depth file. The paths have to be set in the config file. 

### Using the binny wrapper
As shown in the installation description above, binny can be run in a single step, by calling the binny executable. Since most of the configuration is done via the config file, the options are very limited. You can either:
* -c run (submit to a cluster) binny and make a report (-r), or
* -l run (in the current terminal) binny and make a report (-r), or
* -f run (in a tmux session on the frontend) binny *only available in the tmux installation* and make a report (-r), or
* just make a report (-r), or 
* run a dryrun (-d), or 
* unlock a working directory, if a run was killed (-u)
* initialize the conda environmnets only (-i) - you should only need this during the installation. 
It is strongly recommended to **first run a dryrun on a new configuration**, which will tell you within a few seconds and without submission to a cluster whether your chosen steps work together, the input files are where you want them, and your sample file is formatted correctly. In all cases you need the config file as the last argument. 
```
binny -d -r config.yaml
```
You can also set the number of cpus to maximally run at the same time with -t. The defaults (1 for local/frontend runs and 50 for clusters) are reasonable for many settings and if you don't know what this means, you probably don't have to worry. But you may want to increase the numbers for larger datasets or bigger infrastructure, or decrease the numbers to match your environment's constraints.
You can add a name for your main job (-n NAME), e.g.:
```
binny -c -n RUNNAME -r config.yaml
```
Note that spaces in RUNNAME are not allowed and dots will be replaced by underscores.

If you use the tmux version, you can see the tmux process running by typing `tmux ls`. You can also see the progress by checking the stdandard error file `tail RUNNAME_XXXXXXXXXX.stderr`.

Depending on your dataset and settings and your cluster's scheduler, the workflow will take a few minutes to hours to finish. 

### Running snakemake manually
Once metagenomic data and the config file are present, the workflow can be started from the binny directory by the snakemake command:
```
snakemake -s Snakefile --configfile /PATH/TO/YOUR/CONFIGFILE --use-conda
```
If you're using a computing cluster, add your cluster's submission command and the number of jobs you want to maximally run at the same time, e.g.:
```
snakemake -j 50 -s Snakefile --cluster "qsub -l h_rt={resources.runtime},h_vmem=8G -pe smp {threads} -cwd" --configfile /PATH/TO/YOUR/CONFIGFILE --use-conda 
```
This will submit most steps as their own job to your cluster's queue. The same can be achieved with a [cluster configuration](https://snakemake.readthedocs.io/en/stable/executing/cluster-cloud.html#cluster-execution):
```
snakemake -j 50 -s Snakefile --cluster-config PATH/TO/SCHEDULER.config.yaml --cluster "{cluster.call} {cluster.runtime}{resources.runtime} {cluster.mem_per_cpu}{resources.mem} {cluster.threads}{threads} {cluster.partition}" --configfile /PATH/TO/YOUR/CONFIGFILE --use-conda
```
If you want to share the conda installation with colleagues, use the `--conda-prefix` argument of Snakemake
```
snakemake -j 50 -s Snakefile --cluster-config PATH/TO/SCHEDULER.config.yaml --cluster "{cluster.call} {cluster.runtime}{params.runtime} {cluster.mem_per_cpu}{resources.mem} {cluster.threads}{threads} {cluster.partition}" --use-conda --conda-prefix /PATH/TO/YOUR/COMMON/CONDA/DIRECTORY
```
Depending on your dataset and settings, and your cluster's queue, the workflow will take a few minutes to days to finish.

