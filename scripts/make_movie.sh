#!/bin/bash

rundir=$1
outfile=$2
framerate=$3
scale=${4:-1}

if [ "$scale" -eq 1 ]; then
    vf=()
else
    vf=(-vf "scale=iw*$scale:ih*$scale:flags=neighbor")
fi

ffmpeg -framerate "$framerate" -pattern_type glob -i "$rundir/images/*.webp" \
"${vf[@]}" -c:v libx264 -crf 18 -pix_fmt yuv444p "$outfile"