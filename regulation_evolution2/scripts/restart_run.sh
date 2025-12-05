#!/bin/bash
set -e

# Restarts a simulation from its latest backup

rundir=$(realpath "$1")
parfile=$(realpath "$2")
outputrundir=${3:-"$rundir"}
binfile=${4:-"$(dirname "$parfile")/../bin/cell_evolution"}
# Log file (new info is appended at the end)
logfile=${5:-"$outputrundir/log.txt"}
options=("${@:6}")

findfiles() {
  find "$1"/t*.csv -type f -printf "%f\n" | sort
}

lastfile() {
  find "$1"/t*.csv | sort | tail -1
}

celldir="$rundir/celldata_"
lattdir="$rundir/lattice_"
fooddir="$rundir/fooddata_"
cellfiles=$(findfiles "$celldir")
lattfiles=$(findfiles "$lattdir")
foodfiles=$(findfiles "$fooddir")

l1=$(comm -12 <(printf '%s\n' "${cellfiles}") <(printf '%s\n' "${lattfiles}"))
l2=$(comm -12 <(printf '%s\n' "${l1}") <(printf '%s\n' "${foodfiles}"))
# Grabs the string starting from the last 't' and adds a t before it
# did it like this bc using a '\n' character seems not tot work...
lastfile="t${l2##*t}"

mkdir -p "$outputrundir"
colortablefile=$(dirname "$parfile")/colortable.ctb
CMD="$binfile $parfile -name $outputrundir/ -celldatafile $celldir/$lastfile -latticefile $lattdir/$lastfile -fooddatafile $fooddir/$lastfile -colortablefile $colortablefile -existing_dirs ${options[*]}"

echo "Restarting simulation from files labeled: $lastfile" >> "$logfile"
echo "Command used: $CMD" >> "$logfile"
$CMD >> "$logfile" 2>&1
