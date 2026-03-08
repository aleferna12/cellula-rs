This is the template repo for the [cellulars](https://github.com/aleferna12/cellula-rs)
project.

To reduce boilerplate when starting a new cellulars project, 
this repo comes with several quality of life features (which are described below).

### Installation

Either click "Use this template" in the top-right corner of the page, 
or clone this repo directly with:

```commandline
git clone https://github.com/aleferna12/cellulars-template <my_project>
```

After, cd into the local repository and run the tests with:

```commandline
cargo test
```

### Running the model through the command line interface

Two commands are provided out-of-the-box:

#### 1. Run a new simulation

```commandline
cargo run -- run [-l LAYOUT_PNG] [-t TEMPLATE_PARQUET] <CONFIG_TOML>
```

Starts a simulation from a `CONFIG_TOML` file with parameters. 
A default config file with fully documented parameters is provided at 
[`config/1_cell.toml`](https://github.com/aleferna12/cellulars-template/blob/main/config/1_cell.toml). 

Optionally a `LAYOUT_PNG` file can be used to specify cell positions and a 
`TEMPLATE_PARQUET` file can be used to initialize cells with specific parameters.

For example, passing a template file with a migrating cell and a dividing cell together with this layout:

<img alt="squares in B/W" src="https://github.com/aleferna12/cellulars-template/blob/main/tests/fixtures/squares_layout.png" width="200" align="left">
squares_layout.png
<br clear="all">
<br/>

initializes a simulation like this:

<img alt="squares were replaced by cells" src="https://github.com/aleferna12/cellulars-template/blob/main/tests/out/layout_template/images/0000000000.webp" width="200" align="left">
Time step 0<br>
red = migrating<br>
blue = dividing<br>
black = solid object
<br clear="all">

<img alt="cells grew" src="https://github.com/aleferna12/cellulars-template/blob/main/tests/out/layout_template/images/0000000064.webp" width="200" align="left">
Time step 64
<br clear="all">

#### 2. Resume a simulation from backup

```commandline
cargo run -- resume [-t TIME_STEP] [-c CONFIG_FILE] <SIMULATION_DIRECTORY>
```

Restarts a simulation stored at `SIMULATION_DIRECTORY` from the last backup available (or from a specific `TIME_STEP`).
By default, uses the parameters saved in the simulation directory, but those can be overwritten with `CONFIG_FILE`.

For file specifications and more detailed explanations of the arguments, pass the `--help` flag to any of the commands.
For example:

```commandline
cargo run -- run --help
```

### Overriding parameters in the CLI

The parameters specified in the TOML configuration file can be overwritten with environment variables.
The syntax for the environment variables is:

```
CPM__<HEADER1>__<HEADER2>...<HEADER_N>__<PARAMETER_NAME>
```

with double underscores separating the header levels and single underscores replacing spaces.

For example, to overwrite the parameter specifying the number of starting cells to use two cells instead of one, 
you would do:

```bash
CPM__CELL__STARTING_CELLS=2 cargo run -- run config/1_cell.toml
```

on Linux/MacOS, or:

```powershell
$env:CPM__CELL__STARTING_CELLS=2; cargo run -- run config/1_cell.toml
```

on PowerShell.

List parameters should be comma-separated, for example:

```bash
CPM__IO__PLOT__ORDER=spin,border,center cargo run -- run config/1_cell.toml
```

### Examples

Examples for how to use cellulars can be found in the 
[examples folder](https://github.com/aleferna12/cellula-rs/tree/master/cellulars/examples).

### Profiles

Three profiles are available:

1. dev (default)
2. release
3. fastest

It's recommended to run the code with at least `release` for better performance:

```commandline
cargo run --profile release -- run config/1_cell.toml
```

Using `fastest` can further enhance performance, but at a pretty large cost in compilation time. 

### Benchmarks

#### 1. Model benchmarks

There are several [config files](https://github.com/aleferna12/cellulars-template/tree/main/benches/fixtures) 
aimed at benchmarking the model included in this project.

To run a benchmark, use:

```commandline
cargo bench --bench model_bench --profile fastest --no-default-features -- <BENCH_FILE>/<DURATION>mcs
```

which will run a criterion benchmark of the model with parameters `BENCH_FILE` for `DURATION` time steps 
(where `DURATION` is one of 1, 1000, or 10000).

You can also run all benchmarks by omitting the arguments after `--`, 
but be aware this will take a second to run even if using the `fastest` profile.

#### 2. IO benchmarks

Plots can be benchmarked with:

```commandline
cargo bench --bench io_bench --profile fastest --no-default-features
```

which is useful to identify plots driving up execution time of the IO loop.