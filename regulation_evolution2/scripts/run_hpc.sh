#!/bin/bash

#! Use from project directory
#! USAGE: sbatch --array=x-y ./script/run-hpc.sh PATH_TO_PARAMETER_FILE [N_RESTARTS] [OPTIONS]
#! If you change the account or partition they also have to be changed for the restart commands!

#! Name of the job:
#SBATCH -J cell_evolution

#! Which project should be charged:
#SBATCH -A VROOMANS-SL2-CPU
#! Partition to run (icelake or cclake,
#! remember to also change in the restart command and match the modules we load)
#SBATCH -p icelake
#SBATCH -c 1
#SBATCH --time=36:00:00
#! Replace --time with this line when testing
##SBATCH --qos=intr
#SBATCH --nodes=1
#SBATCH --ntasks=1

parfile="$1"
outrunsdir="$2"
restarts=${3:-0}
options=("${@:4}")
#! Assumes script execution from project directory (change this to point to project dir otherwise)
projectdir="$SLURM_SUBMIT_DIR"
#! Colortable file path
colortablefile="$projectdir/data/colortable.ctb"
#! Default directory for each run (depends on the iteration of the job array)
#! If the job is submitted without the array argument, a single run will be created in "runs"
rundir="$outrunsdir/$SLURM_ARRAY_TASK_ID"
#! Default path to the log file of the executable
logfile="$rundir/log.txt"

. /etc/profile.d/modules.sh
#! Uncomment the following lines to load a specific set of modules
## module purge
#! Or rhel7/default-ccl if using cclake
## module load rhel8/default-icl

CMD="$projectdir/bin/cell_evolution $parfile -name $rundir/ -colortablefile $colortablefile ${options[*]}"

mkdir -p "$rundir"
cp "$parfile" "$rundir"
echo "Command used: $CMD" > "$logfile"
$CMD >> "$logfile" 2>&1 &

jobid="$SLURM_ARRAY_JOB_ID"_"$SLURM_ARRAY_TASK_ID"
for ((i=1; i <= restarts; i++))
do
    jobid=$(sbatch -p "$SLURM_JOB_PARTITION" -A VROOMANS-SL2-CPU --parsable --dependency=afterany:"$jobid" -J rest_"$SLURM_ARRAY_TASK_ID"_"$i" --time=36:00:00 "$projectdir"/scripts/restart_run.sh "$rundir" "$parfile" "" "" "" "${options[@]}")
done
wait
