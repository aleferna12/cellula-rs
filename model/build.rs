fn main() {
    // Tell Cargo where the linker can find your shared library
    println!("cargo:rustc-link-search=native=model/lib");

    // Tell Cargo which library to link
    println!("cargo:rustc-link-lib=dylib=kinect");

    // Tell Cargo where to find lib at runtime
    println!("cargo:rustc-link-arg=-Wl,-rpath,model/lib");
}