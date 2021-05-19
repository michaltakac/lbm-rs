use std;
use arrayfire::*;
use crate::FloatNum;
use crate::grid;
use crate::traits::{Distribution, Collision, Physics};
use crate::distribution::D2Q9;

/// Navier-Stokes distributions:
// pub trait NSDistribution {
//     #[inline(always)]
//     fn density(f: &Array<FloatNum>, dims: Dim4) -> Array<FloatNum> {
//         let rho = sum(f, 1);
//         moddims(&rho, dims)
//     }

//     // #[inline(always)]
//     // fn pressure(f: &Array<FloatNum>) -> FloatNum {
//     //     Self::density(f) * Self::c_squ()
//     // }

//     #[inline(always)]
//     fn velocities(f: &Array<FloatNum>, density: &Array<FloatNum>, dims: Dim4) -> (Array<FloatNum>, Array<FloatNum>) {
//         let fex = mul(&transpose(&Self::ex(), false), f, true);
//         let fey = mul(&transpose(&Self::ey(), false), f, true);

//         let ux = moddims(&(sum(&fex, 1) / density), dims);
//         let uy = moddims(&(sum(&fey, 1) / density), dims);
//         (ux, uy)
//     }
// }

// impl NSDistribution for D2Q9 {}
// impl Distribution for D3Q27 {}

/// Single relaxation time (SRT) algorithm
#[derive(Copy, Clone)]
pub struct SingleRelaxationTime {
    pub omega: FloatNum,
    pub re: FloatNum,
    pub nu: FloatNum,
    pub tau: FloatNum,
}

impl<D: Distribution> Collision<D> for SingleRelaxationTime {
    // TODO: implement updating of other parameters, not just omega
    #[inline(always)]
    fn set_omega(&mut self, new_omega: FloatNum) {
        self.omega = new_omega
    }
    #[inline(always)]
    fn collision(&self, f_hlp: &Array<FloatNum>, dist: D, dims: Dim4) -> &H {
        let f_2d = moddims(&f_hlp, dim4!(dims.elements(), D::size()));

        let mut density = D::density(&f_2d, dims);
        let (ux, uy) = D::velocities(&f_2d, &density, dims: Dim4);

        let mut f = f_hlp.clone();
        {
            // Collision
            let u_sq = flat(&(pow(&ux, &(2.0 as FloatNum), false) + pow(&uy, &(2.0 as FloatNum), false)));
            let eu = flat(
                &(&mul(&transpose(&dist.ex(), false), &flat(&ux), true)
                    + &mul(&transpose(&dist.ey(), false), &flat(&uy), true)),
            );
            let feq = flat(&mul(&transpose(&dist.weights(), false), &flat(&density), true))
                * ((1.0 as FloatNum)
                    + (3.0 as FloatNum) * &eu
                    + (4.5 as FloatNum) * (&pow(&eu, &(2.0 as FloatNum), false))
                    - (1.5 as FloatNum) * (&tile(&flat(&u_sq), dim4!(9))));

            // Relaxation step
            f = self.omega * &feq + (1.0 - self.omega) * f_hlp;
        }

        f
    }
}

#[derive(Clone)]
pub struct NavierStokes<D: Distribution, C: Collision<D>> {
    dims: Dim4,
    pub inflow_density: FloatNum,
    pub inflow_accel: FloatNum,
    pub ux: Array<FloatNum>,
    pub uy: Array<FloatNum>,
    pub uz: Array<FloatNum>,
    pub density: Array<FloatNum>,
    collision: C,
    __dist: std::marker::PhantomData<D>,
}

impl<D: Distribution, C: Collision<D>> NavierStokes<D, C> {
    pub fn new(density: FloatNum, accel: FloatNum, dims: Dim4, col: C) -> Self {
        Self {
            dims,
            inflow_density: density,
            inflow_accel: accel,
            ux: constant::<FloatNum>(accel, dims),
            uy: constant::<FloatNum>(0.0, dims),
            uz: constant::<FloatNum>(0.0, dims),
            density: constant::<FloatNum>(density, dims),
            collision: col,
            __dist: std::marker::PhantomData {},
        }
    }
    #[inline(always)]
    pub fn density(&self, f: &Array<FloatNum>, dims: Dim4) -> Array<FloatNum> {
        let rho = sum(f, 1);
        moddims(&rho, self.dims)
    }
    #[inline(always)]
    pub fn velocities(&mut self, solid_nodes: &Array<bool>, f: &Array<FloatNum>, density: &Array<FloatNum>) -> (Array<FloatNum>, Array<FloatNum>) {
        let (ux, uy) = D::velocities(f, density, self.dims);
        
        let zeros = constant(0.0 as FloatNum, self.dims);
        let ux = select(&zeros, solid_nodes, &ux);
        let uy = select(&zeros, solid_nodes, &uy);

        (ux, uy)
    }
}

impl<D: Distribution, C: Collision<D>> Physics
    for NavierStokes<D, C> {
    type Distribution = D;
    #[inline(always)]
    fn collision<H>(&self, f_hlp: &H, dims: Dim4) {
        self.collision.collision(f_hlp, dims);
    }

    fn visualize(&self) {
        unimplemented!();
    }
}