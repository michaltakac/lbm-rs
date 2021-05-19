use arrayfire::*;
use crate::FloatNum;
use crate::geometry;
use crate::grid;

pub trait Simulation {
    fn set_omega(&mut self);
    fn set_initial_conditions(&mut self);
}

pub trait Geometry {
    fn generate(&self, domain: Array<FloatNum>) -> Array<FloatNum>;
}

pub trait Distribution: Sized + Clone + Sync + Send {
    fn c_squ() -> FloatNum;
    fn dims(&self) -> Dim4;
    fn ex(&self) -> Array<FloatNum>;
    fn ey(&self) -> Array<FloatNum>;
    fn weights(&self) -> Array<FloatNum>;
    fn index() -> Array<u32>;
    fn opposite_index() -> Array<u32>;
    fn size() -> u64;
}

pub trait Collision<D: Distribution>: Clone + Sync + Send {
    fn set_omega(&mut self, new_omega: FloatNum);
    fn collision(&self, f_hlp: &Array<FloatNum>, dims: Dim4);
}

pub trait Physics: Clone + Sync + Send {
    type Distribution: Distribution;
    fn collision<FH>(
        &self,
        f_h: &FH,
    ) -> Array<FloatNum>;

    fn visualize(&self);
}