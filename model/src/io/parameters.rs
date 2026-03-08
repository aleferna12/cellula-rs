//! Contains logic related to the simulation parameters.

#![allow(missing_docs)]

use cellulars::constants::{CellIndex, FloatType};
use serde::{Deserialize, Serialize};
use std::path::Path;
use strum_macros::EnumIter;

// When you add parameters, dont forget to document them (and their defaults)
/// Parameters for the model.
///
/// Documentation for each parameter is in `examples/1_cell.toml`
#[derive(Serialize, Deserialize, Clone)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
pub struct Parameters {
    pub general: GeneralParameters,
    pub pond: PondParameters,
    pub cell: CellParameters,
    pub potts: PottsParameters,
    pub io: IoParameters
}

impl Parameters {
    /// Parses parameters from a config file at `path` + env. variables.
    pub fn parse(path: impl AsRef<Path>) -> anyhow::Result<Parameters> {
        let path = path.as_ref();
        log::info!("Reading parameters from {} and environmental variables", path.display());
        let params: Parameters = config::Config::builder()
            .add_source(
                config::File::from(path)
            ).add_source(
                // Converts an env CPM_TIME_STEPS to time-steps
                config::Environment::default()
                    .prefix("CPM")
                    .prefix_separator("__")
                    .separator("__")
                    .convert_case(config::Case::Kebab)
                    .list_separator(",")
                    .with_list_parse_key("io.plot.order")
                    .try_parsing(true)
            ).build()?
            .try_deserialize()?;
        params.check_conflicts()?;
        Ok(params)
    }

    /// Checks for conflicting parameters choices and panics if any are found.
    pub fn check_conflicts(&self) -> anyhow::Result<()> {
        #[cfg(not(feature = "fixed-boundary"))]
        if self.pond.enclose && self.pond.neigh_r > 1 {
            anyhow::bail!(
                "`enclose` can only be used with `neigh-r=1`. \
                 If you need an enclosed pond with larger neighborhoods, enable the `fixed_boundary` feature."
            );
        }
        #[cfg(feature = "fixed-boundary")]
        if !self.pond.enclose {
            anyhow::bail!("`enclose` must be `true` when the `fixed_boundary` feature is enabled")
        }
        Ok(())
    }
}

/// General simulation parameters.
#[derive(Serialize, Deserialize, Clone)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
pub struct GeneralParameters {
    pub time_steps: u32,
    pub seed: Option<u64>
}

/// Parameters determining how a pond is created (see [`pond`](crate::pond));
#[derive(Serialize, Deserialize, Clone)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
pub struct PondParameters {
    pub width: usize,
    pub height: usize,
    #[serde(default = "param_defaults::false_flag")]
    pub enclose: bool,
    pub neigh_r: u8,
}

/// Parameters specifying how cells are created and behave (see [`cell`](crate::my_cell)).
#[derive(Serialize, Deserialize, Clone)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
pub struct CellParameters {
    pub starting_cells: CellIndex,
    pub max_cells: CellIndex,
    pub search_radius: FloatType,
    pub starting_area: u32,
    pub target_area: u32,
    pub div_area: u32,
    #[serde(default = "param_defaults::true_flag")]
    pub divide: bool,
    #[serde(default = "param_defaults::true_flag")]
    pub migrate: bool,
    pub update_period: u32,
}

/// Parameters for the cellular automata update algorithm (see [`potts`](crate::my_potts)).
#[derive(Serialize, Deserialize, Clone)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
pub struct PottsParameters {
    pub boltz_t: FloatType,
    pub size_lambda: FloatType,
    pub chemotaxis_mu: FloatType,
    pub adhesion: AdhesionParameters
}

/// Parameters used in cell adhesion (see [`cellulars::static_adhesion`]).
#[derive(Serialize, Deserialize, Clone)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
pub struct AdhesionParameters {
    pub cell_energy: FloatType,
    pub medium_energy: FloatType,
    pub solid_energy: FloatType,
}

/// Parameters used to control IO operations (see [`io_manager`](crate::io::io_manager)).
#[derive(Serialize, Deserialize, Clone)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
pub struct IoParameters {
    pub outdir: String,
    #[serde(default = "param_defaults::false_flag")]
    pub replace_outdir: bool,
    pub image_period: u32,
    pub info_period: u32,
    pub data: DataParameters,
    pub plot: PlotParameters,
    pub movie: Option<MovieParameters>,
}

/// Parameters used to determine how and when to save data (see [`io_manager`](crate::io::io_manager)).
#[derive(Serialize, Deserialize, Clone)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
pub struct DataParameters {
    pub cells_period: u32,
    pub lattice_period: u32
}

/// Parameters used to display the real-time movie of the simulation.
///
/// Omitting these from the configuration file disables the movie window (same as setting `show` = False).
/// The `movie` feature flag must be on for the movie to be displayed.
#[derive(Serialize, Deserialize, Clone)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
pub struct MovieParameters {
    #[serde(default = "param_defaults::false_flag")]
    pub show: bool,
    pub width: usize,
    pub height: usize,
    pub frame_period: u32
}

/// Parameters using for plotting.
// We flatten the parameters here to allow order to be an env variable
#[derive(Serialize, Deserialize, Clone)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
pub struct PlotParameters {
    pub order: Box<[PlotType]>,
    pub solid_color: String,
    pub medium_color: Option<String>,
    pub center_color: String,
    pub chem_center_color: String,
    pub border_color: String,
    pub area_min_color: String,
    pub area_max_color: String,
    pub chem_min_color: String,
    pub chem_max_color: String,
    pub migrating_color: String,
    pub dividing_color: String,
    pub division_axis_color: String,
    pub division_axis_length: i32,
}


/// A type of plot.
#[derive(Serialize, Deserialize, Clone, EnumIter, Debug)]
#[serde(rename_all = "kebab-case")]
pub enum PlotType {
    /// Cell spin.
    Spin,
    /// Cell center.
    Center,
    /// Cell perceived chemical center.
    ChemCenter,
    /// Cell border.
    Border,
    /// Cell type.
    CellType,
    /// Cell area.
    Area,
    /// Background chemical.
    Chem,
    /// Division axis of cell.
    DivisionAxis
}

// This is a workaround while https://github.com/serde-rs/serde/issues/368 is pending
mod param_defaults {
    pub fn false_flag() -> bool { false }
    pub fn true_flag() -> bool { true }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse() -> anyhow::Result<()> {
        Parameters::parse("config/1_cell.toml")?;
        Ok(())
    }
}
