#!/bin/bash

input=$1
output=$2
scale=$3

ffmpeg -i "$input" vf "scale=iw*$scale:ih*$scale:flags=neighbor" \
-c:v libx264 -crf 18 -pix_fmt yuv444p "$output"