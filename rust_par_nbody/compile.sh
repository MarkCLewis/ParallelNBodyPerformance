# This isn't being used, but I'm recording it here so we can potentially
# add flags to the Cargo build.
rustc -C opt-level=3 -C target-cpu=ivybridge -C codegen-units=1  nbody.rs -o nbody.rust-3.rust_run
