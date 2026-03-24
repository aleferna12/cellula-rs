| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `morpheus --file benchmark/fixtures/dyn/morpheus_500.xml --outdir benchmark/out/morpheus/32/500 --set max_population=32` | 2.194 ± 0.085 | 2.141 | 2.548 | 8.45 ± 0.34 |
| `target/fastest/benchmark 32 500` | 0.260 ± 0.003 | 0.255 | 0.267 | 1.00 |
| `morpheus --file benchmark/fixtures/dyn/morpheus_500.xml --outdir benchmark/out/morpheus/64/500 --set max_population=64` | 2.970 ± 0.029 | 2.924 | 3.051 | 11.44 ± 0.18 |
| `target/fastest/benchmark 64 500` | 0.494 ± 0.004 | 0.488 | 0.502 | 1.90 ± 0.03 |
| `morpheus --file benchmark/fixtures/dyn/morpheus_500.xml --outdir benchmark/out/morpheus/128/500 --set max_population=128` | 4.398 ± 0.062 | 4.337 | 4.571 | 16.93 ± 0.31 |
| `target/fastest/benchmark 128 500` | 0.949 ± 0.010 | 0.929 | 0.971 | 3.65 ± 0.06 |
| `morpheus --file benchmark/fixtures/dyn/morpheus_500.xml --outdir benchmark/out/morpheus/256/500 --set max_population=256` | 7.232 ± 0.174 | 7.047 | 7.599 | 27.84 ± 0.75 |
| `target/fastest/benchmark 256 500` | 1.827 ± 0.022 | 1.801 | 1.894 | 7.03 ± 0.12 |
| `morpheus --file benchmark/fixtures/dyn/morpheus_500.xml --outdir benchmark/out/morpheus/512/500 --set max_population=512` | 13.987 ± 0.762 | 12.691 | 14.983 | 53.86 ± 3.01 |
| `target/fastest/benchmark 512 500` | 3.628 ± 0.121 | 3.526 | 3.997 | 13.97 ± 0.50 |
| `morpheus --file benchmark/fixtures/dyn/morpheus_500.xml --outdir benchmark/out/morpheus/1024/500 --set max_population=1024` | 28.672 ± 1.912 | 25.826 | 32.755 | 110.40 ± 7.48 |
| `target/fastest/benchmark 1024 500` | 7.494 ± 0.180 | 7.270 | 7.843 | 28.85 ± 0.77 |
