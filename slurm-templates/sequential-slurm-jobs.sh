module load slurm

# ====================================================================================
#
#  Simple example for scheduling jobs sequentially
#
# ====================================================================================

# Run first job and get its ID for dependencies when running job 2.
job_1_message=$(sbatch example-job.sbatch)
job_1_id=$(echo $job_1_message | grep -oh "[0-9]*$")

# Change into the proper directory to run job 2. If directory doesn't
# matter or job 2 is in the current directory, you don't need the cd command.
cd path/to/second/directory

# Queue job 2, dependent on the completion of job 1. The only differece here
# is the inclusion of --depend=afterok.
job_2_message=$(sbatch --depend=afterok:${job_1_id} second-example-job.sbatch)
job_2_id=$(echo $job_2_message | grep -oh "[0-9]*$")

# Change into job 3 directory. (only if necessary)
cd path/to/third/directory

# Queue job 3, dependent on the completion of job 2.
job_3_message=$(sbatch --depend=afterok:${job_2_id} third-example-job.sbatch)
job_3_id=$(echo $job_3_message | grep -oh "[0-9]*$")

# To schedule more than 3 jobs, just repeat this pattern, replacing the variable names where necessary.


# ===================================================================================
# 
#  This example includes an array job.
# 
# ===================================================================================

# Run first job and get its ID for dependencies when running job 2.
job_1_message=$(sbatch example-job.sbatch)
job_1_id=$(echo $job_1_message | grep -oh "[0-9]*$")

# Our second job includes an array job. Here, we are using some logic to grab all the
# relevant files for our array, and we are putting the files into a metadata file for
# easy access by the array job.
cd path/to/directory
ls ../01_fastp/raw/ -p | grep -v / | grep ".fastq.gz" | grep "R1" > example-file-list.txt
# Then get the number of array jobs to run:
num_jobs=$(cat example-file-list.txt | wc -l)

# Now, we schedule the array job:
array_message=$(sbatch --depend=afterok:${job_1_id} --array=0-$((num_jobs - 1)) example-array-job.sbatch)
array_id=$(echo $array_message | grep -oh "[0-9]*$")

# When all array jobs finish, run the bam loop
cd path/to/other/directory
job_3_message=$(sbatch --depend=afterok:${array_id} example-job-3.sbatch)
job_3_id=$(echo $job_3_message | grep -oh "[0-9]*$")


# ================================================================================================
#
#  This example schedules two jobs simultaneously, and a third that runs after they both complete.
#
# ================================================================================================

# Schedule first job
job_1_message=$(sbatch example-job.sbatch)
job_1_id=$(echo $job_1_message | grep -oh "[0-9]*$")

# Schedule second job, changing directories (only if needed)
cd path/to/directory
job_2_message=$(sbatch second-example-job.sbatch)
job_2_id=$(echo $job_2_message | grep -oh "[0-9]*$")

# Schedule third job, depending on the completion of both jobs 1 and 2.
cd path/to/other/directory
job_3_message=$(sbatch --depend=afterok:${job_1_id}:${job_2_id} third-example-job.sbatch)
job_3_id=$(echo $job_3_message | grep -oh "[0-9]*$")

# You can imagine adding another job that jobs 1 and 2 both depend on by using
# the --depend=afterok flag in the sbatch commands for both jobs 1 and 2.

