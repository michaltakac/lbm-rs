use std;
use arrayfire::*;
use crate::FloatNum;
use crate::traits::Distribution;

#[derive(Clone)]
#[repr(C)]
pub struct D2Q9 {
    dims: Dim4,
    ex: Array<FloatNum>,
    ey: Array<FloatNum>,
    weights: Array<FloatNum>
}

impl Default for D2Q9 {
    fn default() -> D2Q9 {
        let dims = dim4!(9);
        let t1: FloatNum = 4. / 9.;
        let t2: FloatNum = 1. / 9.;
        let t3: FloatNum = 1. / 36.;
        D2Q9 {
            dims,
            ex: Array::<FloatNum>::new(&[0., 1., 0., -1., 0., 1., -1., -1., 1.], dims),
            ey: Array::<FloatNum>::new(&[0., 0., 1., 0., -1., 1., 1., -1., -1.], dims),
            weights: Array::new(&[t1, t2, t2, t2, t2, t3, t3, t3, t3], dims)
        }
    }
}

impl D2Q9 {
    pub fn new() -> Self {
        Self { ..Default::default() }
    }
}

impl Distribution for D2Q9 {
    #[inline(always)]
    fn c_squ() -> FloatNum {
        1. / 3.
    }
    #[inline(always)]
    fn size() -> u64 {
        9
    }
    fn dims(&self) -> Dim4 {
        self.dims
    }
    fn ex(&self) -> Array<FloatNum> {
        self.ey
    }
    fn ey(&self) -> Array<FloatNum> {
        self.ey
    }
    fn weights(&self) -> Array<FloatNum> {
        self.weights
    }
    #[inline(always)]
    fn index() -> Array<u32> {
        (range::<u32>(dim4!(1, 8), 1) + 1) * Self::size() as u32
    }
    #[inline(always)]
    fn opposite_index() -> Array<u32> {
        (range::<u32>(dim4!(1, 8), 1) + 1) * Self::size() as u32
    }
}
