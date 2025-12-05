#!/bin/bash
set -e

# RUN FROM PROJECT DIRECTORY

rundir=$1  # Run from which to extract the ancestor
parfile=$2
imgfile=$3
latticesize=$4
outdir=${5:-"$rundir/ancestor/"}
food=${6:-"-1"}
options=("${@:7}")

echo "REMEMBER TO PASS '-noevolreg' FOR THE RUNS IF REQUIRED!!!" 1>&2

projectdir="."
ancfile="$outdir/ancestor_template.csv"
cellsfile="$outdir/cells.csv"
lattfile="$outdir/latt.csv"
logfile="$outdir/log.txt"
mkdir -p "$outdir"
cp "$parfile" "$outdir"
python -m scripts ancestor_template "$rundir/celldata_/" "$ancfile" -f "$food"
python -m scripts make_competition -s "$latticesize" "$ancfile" "$imgfile" "$cellsfile" "$lattfile"

echo Running ancestor simulation, output will be redirected to "$logfile"
CMD="$projectdir/bin/cell_evolution $parfile -name $outdir/ -colortablefile $projectdir/data/colortable.ctb -celldatafile $cellsfile -latticefile $lattfile ${options[*]}"

echo "Command used: $CMD" > "$logfile"
$CMD >> "$logfile" 2>&1
