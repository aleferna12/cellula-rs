#!/bin/bash

rundir=$1
prefix=$2  # Prefix of image files like "food", "tau" etc
outdir=${3:-$rundir}
framerate=${4:-20}

ffmpeg -framerate "$framerate" -pattern_type glob -i "$rundir/movie_/$prefix*.png" -c:v libx264 -sws_flags lanczos -pix_fmt yuv420p -crf 18 "$outdir/$(basename "$rundir").mp4"
