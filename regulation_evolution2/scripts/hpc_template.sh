#!/bin/bash

#! Name of the job:
##SBATCH -J PROCESS_NAME

#! Which project should be charged:
#SBATCH -A VROOMANS-SL2-CPU
#! Partition to run (icelake or cclake, remember to also change the module we load to match the partition)
#SBATCH -p icelake
#SBATCH -c 1
#SBATCH --time=36:00:00
#! Replace --time with this line when testing
##SBATCH --qos=intr
#SBATCH --nodes=1
#SBATCH --ntasks=1

. /etc/profile.d/modules.sh
#! Uncomment the following lines to load a specific set of modules
## module purge
#! Or rhel7/default-ccl if using cclake
## module load rhel8/default-icl

# Command goes here
COMMAND_TO_RUN
