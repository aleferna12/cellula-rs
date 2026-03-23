#!/usr/bin/bash

nruns=$1
# Comma separated list
max_pops=$2
# Comma separated list
lat_sizes=$3

echo compiling model
cargo build -p benchmark --profile fastest

echo generating morpheus XML files
IFS=',' read -ra lat_sizes_list <<< "$lat_sizes"
for size in "${lat_sizes_list[@]}"; do
  sed "s/value=\"lat_size lat_size 0\"/value=\"$size $size 0\"/" benchmark/fixtures/morpheus.xml > benchmark/fixtures/dyn/morpheus_"$size".xml
done

echo running benchmark
hyperfine \
  --runs "$nruns" \
  --parameter-list max_pop "$max_pops" \
  --parameter-list size "$lat_sizes" \
  --export-markdown benchmark/results.md \
  --export-json benchmark/results.json \
  'morpheus --file benchmark/fixtures/dyn/morpheus_{size}.xml --outdir benchmark/out/morpheus/{max_pop}/{size} --set max_population={max_pop}' \
  'target/fastest/benchmark {max_pop} {size}'

echo plotting results
python benchmark/plot.py
