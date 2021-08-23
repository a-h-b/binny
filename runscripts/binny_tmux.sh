#! /bin/bash -i

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
VARCONFIG=$DIR/VARIABLE_CONFIG

while IFS=$'\t' read var val; do unset $var ; declare $var="$val" ; done < $VARCONFIG

if [ -z "$MAX_THREADS" ]; then
    MAX_THREADS=50
fi

usage() {
    echo "Usage: $0 [-u|d|c|f|l|i] [-t number] [-r] [-n name] /absolute_path/to/config_file " 1>&2
    echo "       -n <name for main job>, only works with -c and -f" 1>&2
    echo "       -r if set, a report is generated (it's recommended to run -c, -f and -l with -r)" 1>&2
    echo "       -d if set, a dryrun is performed" 1>&2
    echo "       -c if set, the whole thing is submitted to the cluster" 1>&2
    echo "       -i if set, only the conda environments will be installed, if they don't exist" 1>&2
    echo "       -u if set, the working directory will be unlocked (only necessary for crash/kill recovery)" 1>&2
    echo "       -f if set, the main snakemake thread and indivdual rules are run in a tmux session" 1>&2
    echo "       -l if set, the main snakemake thread and indivdual rules are run in the current terminal session" 1>&2
    echo "       -t <max_threads> maximum number of cpus to use for all rules at a time. Defaults to $MAX_THREADS for -c, and to 1 for -l and -f. No effect on -r, -d or -u only." 1>&2

}

while getopts n:t:udflcrhi flag
do
    case $flag in
        i)
            INITIAL=true;;
        u)
            UNLOCK=true;;
        d)
            DRYRUN=true;;
        c)
            CLUSTER=true;;
        n)
            JNAME=$OPTARG;;
        r)
            REPORT=true;;
        f)
            FRONTEND=true;;
        l)
            LAPTOP=true;;
        t)
            THREADS=$OPTARG;;
        h)
            usage
            exit;;
        *)  
            echo "Unimplemented option: -$OPTARG" >&2 
            usage
            exit 1;;
        :) 
            echo "Missing option argument for -$OPTARG" >&2 
            usage
            exit 1;;
        ?)
            usage
            exit
             ;;
    esac
done

shift $((OPTIND-1))

if [ -z "$1" ]; then
    echo "missing input"
    usage
    exit 1
else
    CONFIGFILE=$1
fi


#if the file cannot be found
if [[ !  -e "$1" ]]; then
   echo "Configfile "$1" was not found."
   echo "Provide full path."
   exit 1
fi

if [ "$SNAKEMAKE_VIA_CONDA" = true ]; then
   CONDA_START="conda activate $DIR/conda/snakemake_env"
   CONDA_END="conda deactivate"
   CONDA_END_t="conda deactivate;"
else
   CONDA_START=""
   CONDA_END=""
   CONDA_END_t=""
fi

START_TIME=`date +%s`
NAMEHASH=`echo $START_TIME| cksum | awk '{print $1}'`
if [ -z "$JNAME" ]; then
    JNAME="binny_${NAMEHASH}"
else
    JNAME="${JNAME}_${NAMEHASH}"
fi

if [ "$UNLOCK" = true ]; then
    echo "Unlocking working directory."
    eval $LOADING_MODULES
    eval $CONDA_START
    snakemake $SNAKEMAKE_EXTRA_ARGUMENTS --cores 1 -s $DIR/Snakefile --unlock --configfile $CONFIGFILE
    eval $CONDA_END
elif [ "$DRYRUN" = true ]; then
    echo "Dryrun."
    eval $LOADING_MODULES
    eval $CONDA_START
    snakemake $SNAKEMAKE_EXTRA_ARGUMENTS --cores 1 -s $DIR/Snakefile --dryrun --config sessionName=$JNAME --configfile $CONFIGFILE
    eval $CONDA_END
elif [ "$INITIAL" = true ]; then
    eval $LOADING_MODULES
    eval $CONDA_START
    echo 'Getting MANTIS'
    # mkdir -r ${DIR}/workflow/bin/mantis
    # git clone https://github.com/PedroMTQ/mantis.git ${DIR}/workflow/bin/mantis
    curl -L https://github.com/PedroMTQ/mantis/archive/master.zip --output $DIR/workflow/bin/mantis.zip
    unzip -q $DIR/workflow/bin/mantis.zip -d $DIR/workflow/bin/ && mv $DIR/workflow/bin/mantis-master $DIR/workflow/bin/mantis && rm $DIR/workflow/bin/mantis.zip
    snakemake $SNAKEMAKE_EXTRA_ARGUMENTS --verbose --cores 1 -s $DIR/Snakefile --conda-create-envs-only --use-conda --conda-prefix $DIR/conda --local-cores 1 --configfile $CONFIGFILE
    DB_PATH=`grep "db_path:" $CONFIGFILE | cut -f 2 -d " "`
    temp="${DB_PATH%\"}"
    DB_PATH="${temp#\"}"
    echo $DB_PATH
    if [[ ! "$DB_PATH" = /* ]]
      then
      DB_PATH=${DIR}/$DB_PATH
    fi
    sed -i -e "s|\#nog_hmm_folder\=|nog_hmm_folder=NA|g" \
           -e "s|\#pfam_hmm_folder\=|pfam_hmm_folder=NA|g" \
           -e "s|\#kofam_hmm_folder\=|kofam_hmm_folder=NA|g" \
           -e "s|\#tigrfam_hmm_folder\=|tigrfam_hmm_folder=NA|g" \
           -e "s|\#ncbi_hmm_folder\=|ncbi_hmm_folder=NA|g" \
           -e "s|\#ncbi_dmp_path_folder\=|ncbi_dmp_path_folder=NA|g" \
                      -e "s|\#custom_ref\=path\/to\/hmm/custom1\.hmm|custom_ref=${DIR}/database/hmms/checkm_tf/checkm_filtered_tf.hmm\ncheckm_filtered_tf_weight=0.5\ncustom_ref=${DIR}/database/hmms/checkm_pf/checkm_filtered_pf.hmm\ncheckm_filtered_pf_weight=1|g" \
           ${DIR}/workflow/bin/mantis/MANTIS.config
    for i in ${DIR}/conda/*.yaml; do
      env_name=$(head -n 1 ${i} | cut -d' ' -f2)
      if [[ ${env_name} == 'mantis_env' ]]; then
        mantis_env=$(basename -s .yaml ${i})
      fi
    done
    echo "Setting up Mantis with the CheckM databases"
    conda activate ${DIR}/conda/${mantis_env}
    # Make sure a compiler for cython is available
    if ! [ -x "$(command -v gcc)" ]; then
      conda install -c conda-forge gcc_linux-64 --yes
      $CONDA_PREFIX/etc/conda/activate.d/activate-binutils_linux-64.sh
      $CONDA_PREFIX/etc/conda/activate.d/activate-gcc_linux-64.sh
    fi
    python ${DIR}/workflow/bin/mantis/ setup_databases --chunk_size 1200
    hmmpress database/hmms/checkm_tf/checkm_filtered_tf.hmm
    hmmpress database/hmms/checkm_pf/checkm_filtered_pf.hmm
    python ${DIR}/workflow/bin/mantis/ check_installation
    conda deactivate
    echo "Done."
    exit 0
elif [ "$CLUSTER" = true ]; then
    if [ -z "$THREADS" ]; then
        THREADS=$MAX_THREADS
    fi
    echo "Submitting workflow to cluster - snakemake stays on the frontend, via tmux."
    tmux new -s $JNAME -d
    tmux send-keys -t $JNAME "$LOADING_MODULES >> $JNAME.stdout 2>> $JNAME.stderr" C-m
    tmux send-keys -t $JNAME "$CONDA_START >> $JNAME.stdout 2>> $JNAME.stderr" C-m
    if [ "$REPORT" = true ]; then
    tmux send-keys -t $JNAME "snakemake $SNAKEMAKE_EXTRA_ARGUMENTS --cores $THREADS --jobs $THREADS -s $DIR/Snakefile --keep-going --local-cores 1 --cluster-config $DIR/config/$SCHEDULER.config.yaml --cluster \"{cluster.call} {cluster.runtime}{resources.runtime} {cluster.mem_per_cpu}{resources.mem} {cluster.nodes} {cluster.qos} {cluster.threads}{threads} {cluster.partition} {cluster.stdout}\" --configfile $CONFIGFILE --config sessionName=$JNAME --use-conda --conda-prefix $DIR/conda >> $JNAME.stdout 2>> $JNAME.stderr; snakemake $SNAKEMAKE_EXTRA_ARGUMENTS --cores 1 -s $DIR/Snakefile --configfile $CONFIGFILE --use-conda --conda-prefix $DIR/conda --report report.html >> $JNAME.stdout 2>> $JNAME.stderr; $CONDA_END_t tmux kill-session" C-m
    else
    tmux send-keys -t $JNAME "snakemake $SNAKEMAKE_EXTRA_ARGUMENTS --cores $THREADS --jobs $THREADS -s $DIR/Snakefile --keep-going --local-cores 1 --cluster-config $DIR/config/$SCHEDULER.config.yaml --cluster \"{cluster.call} {cluster.runtime}{resources.runtime} {cluster.mem_per_cpu}{resources.mem} {cluster.threads}{threads} {cluster.nodes} {cluster.qos} {cluster.partition} {cluster.stdout}\" --configfile $CONFIGFILE --config sessionName=$JNAME --use-conda --conda-prefix $DIR/conda >> $JNAME.stdout 2>> $JNAME.stderr; $CONDA_END_t tmux kill-session" C-m
    fi
elif [ "$FRONTEND" = true ]; then
    echo "Running workflow on frontend - don't use this setting except with small datasets and with no more than one run at a time."
    if [ -z "$THREADS" ]; then
        THREADS=1
    fi
    tmux new -s $JNAME -d
    tmux send-keys -t $JNAME "$LOADING_MODULES >> $JNAME.stdout 2>> $JNAME.stderr" C-m
    tmux send-keys -t $JNAME "$CONDA_START >> $JNAME.stdout 2>> $JNAME.stderr" C-m
    if [ "$REPORT" = true ]; then
        tmux send-keys -t $JNAME "snakemake $SNAKEMAKE_EXTRA_ARGUMENTS --cores $THREADS -s $DIR/Snakefile --keep-going --configfile $CONFIGFILE --config sessionName=$JNAME --use-conda --conda-prefix $DIR/conda >> $JNAME.stdout 2>> $JNAME.stderr; snakemake $SNAKEMAKE_EXTRA_ARGUMENTS --cores 1 -s $DIR/Snakefile --configfile $CONFIGFILE --use-conda --conda-prefix $DIR/conda --report report.html >> $JNAME.stdout 2>> $JNAME.stderr; $CONDA_END_t tmux kill-session" C-m
    else
        tmux send-keys -t $JNAME "snakemake $SNAKEMAKE_EXTRA_ARGUMENTS --cores $THREADS -s $DIR/Snakefile --keep-going --configfile $CONFIGFILE --config sessionName=$JNAME --use-conda --conda-prefix $DIR/conda >> $JNAME.stdout 2>> $JNAME.stderr; $CONDA_END_t tmux kill-session" C-m
    fi
elif [ "$LAPTOP" = true ]; then
    echo "Running workflow in current session - don't use this setting except with small datasets and databases."
    JNAME=${JNAME//./_}
    if [ -z "$THREADS" ]; then
        THREADS=1
    fi
    eval $LOADING_MODULES
    eval $CONDA_START
    if [ "$REPORT" = true ]; then
        snakemake $SNAKEMAKE_EXTRA_ARGUMENTS --cores $THREADS -s $DIR/Snakefile --keep-going --configfile $CONFIGFILE --config sessionName=$JNAME --use-conda --conda-prefix $DIR/conda 
        snakemake $SNAKEMAKE_EXTRA_ARGUMENTS --cores $THREADS -s $DIR/Snakefile --configfile $CONFIGFILE --use-conda --conda-prefix $DIR/conda --report report.html 
        eval $CONDA_END
    else
        snakemake $SNAKEMAKE_EXTRA_ARGUMENTS --cores $THREADS -s $DIR/Snakefile --keep-going --configfile $CONFIGFILE --config sessionName=$JNAME --use-conda --conda-prefix $DIR/conda 
        eval $CONDA_END
    fi    
elif [ "$REPORT" = true ]; then
    echo "Writing report."
    eval $LOADING_MODULES
    eval $CONDA_START
    snakemake $SNAKEMAKE_EXTRA_ARGUMENTS --cores 1 -s $DIR/Snakefile --report report.html --configfile $CONFIGFILE --use-conda --conda-prefix $DIR/conda
    eval $CONDA_END
else
    echo "Nothing was done, please give -u, -d, -r, -c, -f, -i, or -l to start anything."
fi


