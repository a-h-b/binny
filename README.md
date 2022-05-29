[![DOI](https://zenodo.org/badge/327396590.svg)](https://zenodo.org/badge/latestdoi/327396590)



# binny

## Quickstart
Here is a quick guide on the installation and test run of binny. Please check out the longer description below to set up binny on a cluster environment.
Make sure you have conda (and optionally, recommended mamba) installed.

1) Clone this repository with git
```
git clone https://github.com/a-h-b/binny.git
cd binny
```

2) Optional configuration 
You can set the following fields in the config file `config/config.init.yaml`, or adjust the paths mentioned below. Or you can keep the default settings: By default, all dependencies will be put into `conda` in the `binny` directory.

a) Optional: Change config file with path for conda environments.
```
my_conda_path="absolute/path/to/conda/dir" #adjust path here
sed -i -e "s|conda_source: \"\"|conda_source: \"${my_conda_path}\"|g" config/config.*.yaml
```

b) Optional: If you already have environments for Snakemake, Prokka and/or Mantis, you can add them to the config file and binny will use them, e.g.:
(for Snakemake, if you have it already installed in your PATH set the param to `"in_path"`. By default, a new Snakemake env will be created in the binny dir)
```
my_prokka_env="absolute/path/to/prokka/env.yaml" or "my_named_prokka_env" #choose path/env here
my_mantis_env="absolute/path/to/mantis/env.yml" or "my_named_mantis_env" #choose path/env here
my_snakemake_env="absolute/path/to/snakemake/env" or "my_named_snakemake_env" or "in_path" #choose path/env here
sed -i -e "s|prokka_env: \"\"|prokka_env: \"${my_prokka_env}\"|g" \
       -e "s|mantis_env: \"\"|mantis_env: \"${my_mantis_env}\"|g" \
       -e "s|snakemake_env: \"\"|snakemake_env: \"${my_snakemake_env}\"|g" config/config.*.yaml
```

3) Install the snakemake and conda environments, and databases
```
./binny -i config/config.init.yaml 
```

4) Test run
```
./binny -l -n "TESTRUN" -r config/config.test.yaml
```
If all goes well, binny will run in the current session, load the conda environments, and make and fill a directory called `test_output`. A completed run should contain four fasta files with one bin each in `test_output/bins`. 
If you don't want to see binny's guts at this point, you can also run this with the -c or -f settings to submit to your cluster or start a tmux session (see How to run binny below). 

Please see the comments in `config/config.default.yaml` and the longer descriptions below to set up your own runs.

## Adjusting the VARIABLE_CONFIG
* Make sure to separate variable name and your setting with a tab
* SNAKEMAKE_VIA_CONDA - set this to true, if you don't have snakemake in your path and want to install it via conda. Leave empty, if you don't need an additional snakemake.
* SNAKEMAKE_EXTRA_ARGUMENTS - if you want to pass additional arguments to snakemake, put them here (e.g. --latency-wait=320 for slower file systems). Leave empty usually.
* LOADING_MODULES - insert a bash command to load modules, if you need them to run conda. Leave empty, if you don't need to load a module.
* SUBMIT_COMMAND - insert the bash command you'll usually use to submit a job to your cluster to run on a single cpu for a few days. You only need this, if you want to have the snakemake top instance running in a submitted job (`-c` option without `-x`). You alternatively have the option to run the snakemake top instance on the frontend via tmux (`-c -x`). Leave empty, if you want to use [tmux](https://github.com/tmux/tmux/wiki) and have it installed.
* BIND_JOBS_TO_MAIN - if you use the option to run the snakemake top instance in a submitted job and need to bind the other jobs to the same node, you can set this option to true. See FAQ below for more details.
* NODENAME_VAR - if you use the BIND_JOBS_TO_MAIN option, you need to let dadasnake know, how to access the node name (e.g.SLURMD_NODENAME on slurm).
* SCHEDULER - insert the name of the scheduler you want to use (currently `slurm`, `slurm_simple` or `sge`). This determines the cluster config given to snakemake, e.g. the cluster config file for slurm is config/slurm.config.yaml . Also check that the settings in this file are correct. If you have a different system, contact us ( https://github.com/a-h-b/binny/issues ).
* MAX_THREADS - set this to the maximum number of cores you want to be using in a run. If you don't set this, the default will be 50 (which is more than will be used). Users can override this setting at runtime.



## How to run binny
To run the binny, you need a config file, plus data: 
* The config file (in yaml format) is read by Snakemake to determine the inputs, arguments and outputs. 
* You need contigs in a fasta file and the alignments of metagenomic reads in bam format, or the contigs plus a tab-separated depth file. The paths have to be set in the config file. 

### Using the binny wrapper
As shown in the installation description above, binny can be run in a single step, by calling the binny executable. Since most of the configuration is done via the config file, the options are very limited. You can choose:
* `-c` : run binny *submitted to a cluster* and make a report (`-r`), or
* `-c -x` : run binny *in a tmux session, submitting jobs to a cluster* and make a report (`-r`), or
* `-l` : run binny *in the current terminal* and make a report (`-r`), or
* `-f` : run binny *in a tmux session on the frontend* (same as `-l -x`) and make a report (-r), or
* `-r` : just make a *report*, or 
* `-d` : run a *dryrun*, or 
* `-u` : *unlock* the output directory named in the config file (**only necessary, if a run was killed**)
* `-i` : *initialize the conda environments only* (**you should only need this during the installation**) 
It is strongly recommended to **first run a dryrun on a new configuration**, which will tell you within a few seconds and without submission to a cluster whether your chosen steps work together, the input files are where you want them, and your sample file is formatted correctly. In all cases you need the config file as the last argument. 
```
binny -d config.yaml
```
You can also set the number of cpus to maximally run at the same time with `-t`. The defaults (1 for local/frontend runs and 50 for clusters) are reasonable for many settings and if you don't know what this means, you probably don't have to worry. But you may want to decrease the numbers to match your environment's constraints.
You can add a name for your main job (`-n NAME`), e.g.:
```
binny -c -n MYBINNYrun -r config.yaml
```
Note that spaces in NAME are not allowed and dots will be replaced by underscores.

If you use the tmux flag, you can see the tmux process running by typing `tmux ls`. You can also see the progress by checking the stdandard error file `tail NAME_XXXXXXXXXX.stderr`.

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
