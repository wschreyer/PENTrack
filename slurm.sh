#!/bin/sh

# this is an example script showing how to run PENTrack on a computing cluster using the SLURM task scheduler

#SBATCH --array=1-250  # run multiple simulations in parallel
#SBATCH --time=4:00:00 # set maximum time simulations are allowed to run (should be as close as possible to the actual running time reported by sacct command)
###SBATCH --mem=4000M  # reserve additional memory (e.g. when using 3D maps of field, only reserve as much as you really need as reported by sacct command)
###SBATCH --account=   # set an account, e.g. in case you have an account with priority access to resources


# combine Job ID and array task ID to a unique ID that is passed to PENTrack
ID=$SLURM_ARRAY_TASK_ID  # Slurm array task index
JOB=$SLURM_ARRAY_JOB_ID  # Slurm job ID

# run PENTrack with the created ID, the default config, and the default output directory
./PENTrack $JOB$ID ./in/config.in ./out/
