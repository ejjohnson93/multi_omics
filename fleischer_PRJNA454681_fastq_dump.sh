#!/bin/bash -l
# Use the current working directory
#SBATCH -D ./
# Use the current environment for this job.
#SBATCH --export=ALL
# Define job name
#SBATCH -J fastq_dump
# Define a standard output file. When the job is running, %u will be replaced by user name,
# %N will be replaced by the name of the node that runs the batch script, and %j will be replaced by job id number.
#SBATCH -o fastq_dump.%u.%N.%j.out
# Define a standard error file
#SBATCH -e fastq_dump.%u.%N.%j.err
# Request the partition
#SBATCH -p nodes
# Request the number of nodes
#SBATCH -N 1
# Request the number of cores
#SBATCH -n 8
# Specify time limit in format a-bb:cc:dd, where a is days, b is hours, c is minutes, and d is seconds.
#SBATCH -t 3-00:00:00
# Request the memory on the node or request memory per core
# PLEASE don't set the memory option as we should use the default memory which is based on the number of cores
##SBATCH --mem-per-cpu=9000M
# Insert your own username to get e-mail notifications (note: keep just one "#" before SBATCH)
#SBATCH --mail-user=ejohn16@liverpool.ac.uk
# Notify user by email when certain event types occur
#SBATCH --mail-type=ALL
#
#
# Load your own modules
module purge
module load apps/sratoolkit/2.9.6-1
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

# If you use all of the cores with mpirun, you do not need
# to specify how many MPI processes to use - that is the default
# the ret flag is the return code, so you can spot easily if your code failed.
#mpirun  $EXEC


# download fastq files using fastq-dump
# PRJNA454681
# https://www.ncbi.nlm.nih.gov/sra?linkname=bioproject_sra_all&from_uid=454681

# Save temp files to temp directory instead of home directory
echo '/repository/user/main/public/root = "/tmp/users/ejohn16/"' > $HOME/.ncbi/user-settings.mkfg

OUTDIR=$SLURM_SUBMIT_DIR/skin_fleischer_PRJNA454681
mkdir -p ${OUTDIR}
#SRR7093809 - SRR7093951

echo "Downloading files for bioproject: PRJNA454681. Fleischer, ageing skin dataset."

for i in $(seq 809 951)
do
	echo "SRR7093${i}"

	fastq-dump --split-files \
		   --outdir ${OUTDIR} \
		   SRR7093${i}

	gzip ${OUTDIR}/*.fastq

done

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
