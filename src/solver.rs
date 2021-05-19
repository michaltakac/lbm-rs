use crate::boundary;
use crate::grid::{StructuredGrid};
use crate::traits::{Distribution, Physics};
use crate::FloatNum;
use arrayfire::*;
// use time;

/// Lattice-Boltzmann Solver state
pub struct Solver<P: Physics, D: Distribution> {
    distribution: D,
    grid: StructuredGrid<P::Distribution>,
    pub bcs: boundary::Handler<P::Distribution>,
    physics: P,
    f: Array<FloatNum>, // Distribution functions
    f_hlp: Array<FloatNum>,
}

impl<P: Physics, D: Distribution> Solver<P, D> {
    /// Create a new solver from a `grid` and `physics`.
    pub fn new(distribution: D, grid: StructuredGrid<D>, bcs: boundary::Handler<D>, physics: P) -> Self {
        Self {
            distribution,
            grid,
            bcs,
            physics,
            f: constant(
                0.0 as FloatNum,
                dim4!(D::size() as u64 * P::Distribution::size() as u64),
            ),
            f_hlp: constant(
                0.0 as FloatNum,
                dim4!(D::size() as u64 * P::Distribution::size() as u64),
            ),
        }
    }

    /// Initialize distributions
    pub fn initialize(&mut self) {
        // Start in equilibrium state
        let u_sq: Array<FloatNum> =
            flat(&(pow(&self.physics.ux, &(2.0 as FloatNum), false) + pow(&&self.physics.uy, &(2.0 as FloatNum), false)));
        let eu: Array<FloatNum> = flat(
            &(&mul(&transpose(&self.distribution.ex(), false), &flat(&self.physics.ux), true)
                + &mul(&transpose(&self.distribution.ey(), false), &flat(&self.physics.uy), true)),
        );

        self.f = flat(&mul(&transpose(&self.distribution.weights(), false), &flat(&self.physics.density), true))
            * ((1.0 as FloatNum)
                + (3.0 as FloatNum) * &eu
                + (4.5 as FloatNum) * (&pow(&eu, &(2.0 as FloatNum), false))
                - (1.5 as FloatNum) * (&tile(&flat(&u_sq), dim4!(9))));
        
        self.f_hlp = self.f
    }

    /// Mutable reference to the distribution function
    fn f_mut(&mut self) -> &mut Array<FloatNum> {
        &mut self.f
    }

    /// Reference to the distribution function
    fn f_ref(&self) -> &Array<FloatNum> {
        &self.f
    }

    /// Streaming step
    fn streaming(&mut self) {
        let mut f_hlp = std::mem::replace(&mut self.f_hlp, self.f.clone());
        self.f_hlp = view!(f_hlp[P::Distribution::neighbors_index()]);
    }

    /// Collision step
    fn collision(&mut self) {
        let mut f = std::mem::replace(&mut self.f, Default::default());
        
        self.physics.collision(f, self.grid.dimensions);
    }

    /// Applies boundary conditions
    fn apply_boundary_conditions(&mut self) {
        let f_hlp = self.f_hlp.copy();
        let mut f = std::mem::replace(&mut self.f, f_hlp);

        self.f = f;
        self.bcs.apply(self.f, self.f_hlp, self.physics);
        // f.par_chunks_mut(P::Distribution::size())
        //     .zip(self.grid.par_ids())
        //     .for_each(|(f, c)| {
        //         let r = self.bcs.apply(
        //             &f,
        //             &self.f_hlp,
        //             |v, n: P::Distribution| v[n.value()],
        //             |v, n| v[Self::f_idx(c, n)],
        //             self.grid.x(c),
        //         );

        //         if let Some(r) = r {
        //             for n in P::Distribution::all() {
        //                 f[n.value()] = r.as_ref()[n.value()];
        //             }
        //         }
        //     });
    }

    /// Executes `n_it` iterations writing output every `n_out` iterations.
    pub fn run(&mut self, iter: usize, write_output: bool) {
        use time::Duration;
        let d = Duration::span(|| {
            let d = Duration::span(|| self.streaming());
            if write_output {
                self.substep("streaming", d);
            }

            let d = Duration::span(|| self.collision());
            if write_output {
                self.substep("collision", d);
            }

            let d = Duration::span(|| self.apply_boundary_conditions());
            if write_output {
                self.substep("bcs", d);
            }

            // if write_output {
            //     let d = Duration::span(|| self.visualize());
            //     self.substep("visualize", d);
            // }
        });
        if write_output {
            self.step(iter, d);
        }
    }

    /// Prints line info of a whole iteration step
    fn step(&self, n_it: usize, duration: time::Duration) {
        println!(
            "#{} | duration: {} ms",
            n_it,
            duration.num_milliseconds()
        );
    }
    /// Prints line info of an iteration sub-step
    fn substep(&self, name: &str, duration: time::Duration) {
        println!(
            "# [{}] | duration: {} \u{03BC}s",
            name,
            duration.num_microseconds().unwrap()
        );
    }
    /// Writes the solution to a VTK file.
    fn visualize(&self, iter: usize) {
        unimplemented!();
    }
}
