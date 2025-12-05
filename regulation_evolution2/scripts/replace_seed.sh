#!/bin/bash
set -e

# Replaces the seed passed in the par file with the seed actually used in the simulation (printed in the log file)

logfile="$1"
parfile="$2"

actual_seed=$(sed -rn 's/Seed for the random generator is: ?([0-9]+)/\1/p' "$logfile")
used_pattern='(rseed ?= ?)([0-9-]+)'
used_seed=$(sed -rn "s/$used_pattern/\2/p" "$parfile")

if [ "$used_seed" == '-1' ]; then
  sed -ri "s/$used_pattern/\1$actual_seed/" "$parfile"
elif [ "$used_seed" == '' ]; then
  echo "rseed = $actual_seed" >> "$parfile"
else
  if [ "$used_seed" != "$actual_seed" ]; then
    echo "Error: seed in parameter file was set as '$used_seed', but seed reported in log file is '$actual_seed'. Are these files from the same simulation?" 1>&2
    exit 1
  fi
fi
