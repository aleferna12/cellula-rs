use clap::Parser;
use evo_cpm::model::Model;
use evo_cpm::parameters::Parameters;
use evo_cpm::utils::TEST_SEED;
// TODO: profile inlining

fn main() {
    let mut model = Model::new(Parameters::parse_from(["", "--seed", &TEST_SEED.to_string()]));
    model.welcome();
    model.setup();
    model.run(100)
}
