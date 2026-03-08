This is the template repo for the [cellulars](https://github.com/aleferna12/cellula-rs)
project.

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

<img src="https://github.com/aleferna12/cellulars-template/blob/main/tests/fixtures/squares_layout.png" width="200">
squares_layout.png
<div style="clear: both;"></div>

initializes a simulation like this:

<img src="https://github.com/aleferna12/cellulars-template/blob/main/tests/out/layout_template/images/0000000000.webp" width="200">
time step 0
<div style="clear: both;"></div>

<img src="https://github.com/aleferna12/cellulars-template/blob/main/tests/out/layout_template/images/0000000064.webp" width="200">
time step 64
<div style="clear: both;"></div>


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