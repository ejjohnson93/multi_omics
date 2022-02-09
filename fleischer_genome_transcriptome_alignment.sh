#!/bin/bash -l
# Use the current working DATA_INectory
#SBATCH -D ./
# Use the current environment for this job.
#SBATCH --export=ALL
# Define job name
#SBATCH -J alignments
# Define a standard output file. When the job is running, %u will be replaced by user name,
# %N will be replaced by the name of the node that runs the batch script, and %j will be replaced by job id number.
#SBATCH -o Alignments_breast_cancer.%u.%N.%j.out
# Define a standard error file
#SBATCH -e Alignments_breast_cancer.%u.%N.%j.err
# Request the partition
#SBATCH -p nodes
# Request the number of nodes
#SBATCH -N 2
# Specify the number of tasks per node
#SBATCH --ntasks-per-node=8
# Specify the number of tasks
##SBATCH --ntasks=16 #this option is not set as we have already set --ntasks-per-node
# Request the number of cpu per task
#SBATCH --cpus-per-task=5
# This asks for 3 days
#SBATCH -t 3-00:00:00
# Specify memory per core
#SBATCH --mem-per-cpu=9000M
# Insert your own username to get e-mail notifications
#SBATCH --mail-user=ejohn16@liverpool.ac.uk
# Notify user by email when certain event types occur
#SBATCH --mail-type=ALL
#
# Load your own modules
module purge
module load apps/anaconda3/2019.10-bioconda
module load apps/samtools/1.9
# List all modules
module list
#
#
echo =========================================================
echo SLURM job: submitted  date = `date`
date_start=`date +%s`

echo Executable file:
echo MPI parallel job.
echo -------------
echo Job output begins
echo -----------------
echo

hostname

echo "Print the following environmental variables:"
echo "Job name                     : $SLURM_JOB_NAME"
echo "Job ID                       : $SLURM_JOB_ID"
echo "Job user                     : $SLURM_JOB_USER"
echo "Job array index              : $SLURM_ARRAY_TASK_ID"
echo "Submit directory             : $SLURM_SUBMIT_DIR"
echo "Temporary directory          : $TMPDIR"
echo "Submit host                  : $SLURM_SUBMIT_HOST"
echo "Queue/Partition name         : $SLURM_JOB_PARTITION"
echo "Node list                    : $SLURM_JOB_NODELIST"
echo "Hostname of 1st node         : $HOSTNAME"
echo "Number of nodes allocated    : $SLURM_JOB_NUM_NODES or $SLURM_NNODES"
echo "Number of processes          : $SLURM_NTASKS"
echo "Number of processes per node : $SLURM_TASKS_PER_NODE"
echo "Requested tasks per node     : $SLURM_NTASKS_PER_NODE"
echo "Requested CPUs per task      : $SLURM_CPUS_PER_TASK"
echo "Scheduling priority          : $SLURM_PRIO_PROCESS"


echo "Running parallel job:"

### Script to carry out qc and alignments for supplied datasets, indexes and base directory provided here
# fastq files need to end in _1 or _1 & _2 to designate whether they are single or paired-end reads

BASEDIR=$SLURM_SUBMIT_DIR
HISAT2_INDEX=$SLURM_SUBMIT_DIR/genomes/human/hg19/hg19_hisat2_index/hg19
SALMON_INDEX=$SLURM_SUBMIT_DIR/genomes/human/hg19/hg19_salmon_index/hg19

# Name of directory to analyse provided here

i=skin_fleischer_PRJNA454681

# Assigns variables and creates directories for data

DATA_IN=${BASEDIR}/${i}
HISAT2_OUTPUT=${DATA_IN}/hisat2
SALMON_OUTPUT=${DATA_IN}/salmon_quant

mkdir -p ${HISAT2_OUTPUT} ${SALMON_OUTPUT}

cd ${DATA_IN}

# Alignment, HISAT2 and salmon

#If files ending in _2 exist in the target directory then the reads are paired...

if ls ${DATA_IN}/*_2.* 1> /dev/null 2>&1; then
    echo "Sequences are paired-end, beginning alignments with HISAT2 and salmon..."
    echo
    
    # HISAT2 alignment

    for k in ${DATA_IN}/*_1.fastq.gz; do

    SEQ_ID=$(basename ${k} _1.fastq.gz)
    echo "Aligning reads for ${SEQ_ID} with HISAT2..."
    hisat2 -x ${HISAT2_INDEX} \
        -p $SLURM_CPUS_PER_TASK \
        -1 "${DATA_IN}/${SEQ_ID}_1.fastq.gz" \
        -2 "${DATA_IN}/${SEQ_ID}_2.fastq.gz" \
        | samtools view -b | samtools sort > ${HISAT2_OUTPUT}/${SEQ_ID}_hisat2.bam

    done
    
    # salmon alignment 

    for j in ${DATA_IN}/*_1.fastq.gz; do

    SEQ_ID=$(basename ${j} _1.fastq.gz)
    echo "Salmon quasi-alignment for ${SEQ_ID}..."
        salmon quant -i ${SALMON_INDEX} \
            -l A \
            -1 "${DATA_IN}/${SEQ_ID}_1.fastq.gz" \
            -2 "${DATA_IN}/${SEQ_ID}_2.fastq.gz" \
            -o "${SALMON_OUTPUT}/${SEQ_ID}_quant" \
            -p 4

    done
    
    #If there were no files ending in _2 (previous step) but files ending in _1 exist in the target directory then the reads are single read...

elif ls ${DATA_IN}/*_1.* 1> /dev/null 2>&1; then
    echo "Reads are single-read, beginning alignment with HISAT2..."
    echo

    for k in ${DATA_IN}/*_1.fastq.gz; do

    SEQ_ID=$(basename ${k} _1.fastq.gz)
    echo "Aligning reads for ${SEQ_ID} with HISAT2..."
    hisat2 -x ${HISAT2_INDEX} \
        -p $SLURM_CPUS_PER_TASK \
        -U "${DATA_IN}/${SEQ_ID}_1.fastq.gz" \
        | samtools view -b | samtools sort > ${HISAT2_OUTPUT}/${SEQ_ID}_hisat2.bam

    done

    for j in ${DATA_IN}/*_1.fastq.gz; do

    SEQ_ID=$(basename ${j} _1.fastq.gz)
    echo "Salmon quasi-alignment for ${SEQ_ID}..."
        salmon quant -i ${SALMON_INDEX} \
            -l A \
            -r "${DATA_IN}/${SEQ_ID}_1.fastq.gz" \
            -o "${SALMON_OUTPUT}/${SEQ_ID}_quant" \
            -p 4

    done

else 

    # No files ending in _1 or _2 suggest the files are incorrectly named or the target directory is wrong, script determines no reads available..

    echo "No reads available..."

fi 


ret=$?


echo
echo ---------------
echo Job output ends
date_end=`date +%s`
seconds=$((date_end-date_start))
minutes=$((seconds/60))
seconds=$((seconds-60*minutes))
hours=$((minutes/60))
minutes=$((minutes-60*hours))
echo =========================================================
echo SLURM job: finished   date = `date`
echo Total run time : $hours Hours $minutes Minutes $seconds Seconds
echo =========================================================
exit $ret
