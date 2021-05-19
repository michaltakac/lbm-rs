# Lattice-Boltzmann Method on GPUs in Rust

This project is a fork of [gnzlbg's](https://github.com/gnzlbg/) implementation of [Lattice Boltzmann solver](https://github.com/gnzlbg/lbm-rs). I liked the general structure of his code.

This implementation uses ArrayFire library for computations on GPUs. The VTK writer code was removed for and instead real-time visualization with ArrayFire's Forge is used.
