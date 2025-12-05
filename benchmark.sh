#!/usr/bin/bash
set -e

hyperfine --export-markdown perf.md "bash regulation_evolution2/scripts/run.sh regulation_evolution2/data/perf.par regulation_evolution2/ 0 0" \
  "cargo run --profile fastest -- run perf.toml" --runs 10

echo Finished runs. Results are in "perf.md"