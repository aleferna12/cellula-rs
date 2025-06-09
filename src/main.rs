use clap::Parser;
use evo_cpm::model::Model;
use evo_cpm::parameters::Parameters;

fn main() {
    let model = Model::new(Parameters::parse());
    model.welcome();
}
