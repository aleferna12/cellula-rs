| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `morpheus --file benchmark/fixtures/dyn/morpheus_1000.xml --outdir benchmark/out/morpheus/64/1000 --set max_population=64` | 2.511 ± 0.070 | 2.435 | 2.629 | 5.01 ± 1.93 |
| `target/fastest/benchmark 64 1000` | 0.501 ± 0.192 | 0.439 | 1.049 | 1.00 |
| `morpheus --file benchmark/fixtures/dyn/morpheus_1000.xml --outdir benchmark/out/morpheus/128/1000 --set max_population=128` | 4.267 ± 0.135 | 4.095 | 4.571 | 8.51 ± 3.28 |
| `target/fastest/benchmark 128 1000` | 0.828 ± 0.005 | 0.820 | 0.836 | 1.65 ± 0.63 |
| `morpheus --file benchmark/fixtures/dyn/morpheus_1000.xml --outdir benchmark/out/morpheus/256/1000 --set max_population=256` | 7.604 ± 0.150 | 7.398 | 7.827 | 15.17 ± 5.83 |
| `target/fastest/benchmark 256 1000` | 1.581 ± 0.012 | 1.567 | 1.612 | 3.15 ± 1.21 |
| `morpheus --file benchmark/fixtures/dyn/morpheus_1000.xml --outdir benchmark/out/morpheus/512/1000 --set max_population=512` | 14.097 ± 0.363 | 13.570 | 14.839 | 28.12 ± 10.82 |
| `target/fastest/benchmark 512 1000` | 3.108 ± 0.046 | 3.043 | 3.211 | 6.20 ± 2.38 |
| `morpheus --file benchmark/fixtures/dyn/morpheus_1000.xml --outdir benchmark/out/morpheus/1024/1000 --set max_population=1024` | 63.682 ± 40.304 | 26.599 | 129.021 | 127.01 ± 94.02 |
| `target/fastest/benchmark 1024 1000` | 7.796 ± 1.042 | 6.902 | 10.662 | 15.55 ± 6.32 |
| `morpheus --file benchmark/fixtures/dyn/morpheus_1000.xml --outdir benchmark/out/morpheus/2048/1000 --set max_population=2048` | 149.659 ± 17.502 | 128.871 | 173.127 | 298.49 ± 119.79 |
| `target/fastest/benchmark 2048 1000` | 15.328 ± 0.876 | 14.470 | 17.017 | 30.57 ± 11.87 |
