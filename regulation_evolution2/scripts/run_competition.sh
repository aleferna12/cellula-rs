#!/bin/bash
set -e

templatecellfile1=$1
templatecellfile2=$2
templatefile=$3
cellcompfile=$4
lattcompfile=$5
parfile=$6
rundir=$7
imagefile=${8:-data/competition/scatter.png}
food=${9:-"-1"}
options=("${@:10}")

python -m scripts make_templates "$templatefile" -i "$templatecellfile1" "$templatecellfile2" -f "$food"
python -m scripts make_competition -s 2000 "$templatefile" "$imagefile" "$cellcompfile" "$lattcompfile"

mkdir -p "$rundir"
CMD="./bin/cell_evolution $parfile -name $rundir/ -colortablefile ./data/colortable.ctb -noevolreg -groupextinction -celldatafile $cellcompfile -latticefile $lattcompfile ${options[*]}"
echo "Command used: $CMD" > "$rundir/log.txt"
$CMD >> "$rundir/log.txt" 2>&1 &
echo "PID $!: $CMD"
