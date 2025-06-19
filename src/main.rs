use clap::Parser;
use evo_cpm::model::Model;
use evo_cpm::parameters::Parameters;
use evo_cpm::io;
// TODO: profile inlining core functions

fn main() {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();
    let parameters = Parameters::parse();
    io::welcome(&parameters);
    let mut model = Model::new(parameters);
    model.setup();
    io::simulation_frame(&model.env).save("./test1.png").unwrap();
    model.run(1000);
    io::simulation_frame(&model.env).save("./test2.png").unwrap()
}
