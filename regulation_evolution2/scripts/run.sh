#!/bin/bash

#! Use from project directory
#! USAGE: bash scripts/run.sh PARFILE RUNSDIR FIRSTRUN LASTRUN PATH_TO_PARAMETER_FILE [OPTIONS]

parfile="$1"
outrunsdir="$2"
firstrun="$3"
lastrun="$4"
options=("${@:5}")

#! Assumes script execution from project directory (change this to point to project dir otherwise)
projectdir="./regulation_evolution2/"
#! Colortable file path
colortablefile="$projectdir/data/colortable.ctb"

for i in $(seq "$firstrun" "$lastrun")
do
  #! Default directory for each run (depends on the iteration of the job array)
  rundir="$outrunsdir/$i"
  #! Default path to the log file of the executable
  logfile="$rundir/log.txt"

  CMD="$projectdir/bin/cell_evolution $parfile -name $rundir/ -colortablefile $colortablefile ${options[*]}"

  mkdir -p "$rundir"
  cp "$parfile" "$rundir"
  echo "Command used: $CMD" > "$logfile"
  $CMD >> "$logfile" 2>&1 &
  echo "PID $!: $CMD"
done
wait
