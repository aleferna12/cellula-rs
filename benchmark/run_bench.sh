#!/usr/bin/bash

cargo build -p benchmark --profile fastest

hyperfine \
  --warmup 1 \
  --runs 3 \
  --export-markdown benchmark/results.md \
  --export-json benchmark/results.json \
  'target/fastest/benchmark' \
  'morpheus --file benchmark/fixtures/morpheus.xml'