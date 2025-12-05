#!/usr/bin/bash
set -e

cargo build --profile fastest
cd regulation_evolution2;
cmake -DCMAKE_BUILD_TYPE=Release .
make

echo Finished compilation
