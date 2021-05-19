use crate::streaming::{stream_2d, stream_3d};
use crate::traits::Distribution;
use crate::FloatNum;
use arrayfire::*;

/// Rectangular grid of up to three dimensions.
#[derive(Clone, Debug)]
pub struct StructuredGrid<D: Distribution> {
    pub dimensions: Dim4,
    pub x: u64,
    pub y: u64,
    pub z: u64,
    __dist: std::marker::PhantomData<D>
}

impl<D: Distribution> StructuredGrid<D> {
    #[inline(always)]
    pub fn new(x: u64, y: u64, z: u64) -> Self {
        let dimensions = match z {
            0 => dim4!(x, y),
            1 => dim4!(x, y),
            _ => dim4!(x, y, z),
        };
        Self {
            x,
            y,
            z,
            dimensions,
            __dist: std::marker::PhantomData {}
        }
    }

    #[inline(always)]
    pub fn size(&self) -> u64 {
        match self.dimensions.ndims() {
            2 => self.x * self.y,
            3 => self.x * self.y * self.z,
        }
    }

    #[inline(always)]
    pub fn main_index(&self) -> Array<FloatNum> {
        moddims(
            &range::<FloatNum>(dim4!(self.size() * D::size()), 0),
            dim4!(self.x, self.y, D::size()),
        )
    }

    /// Returns the neighboring indices
    #[inline(always)]
    pub fn neighbors_index(&self) -> Array<FloatNum> {
        match self.dimensions.ndims() {
            2 => flat(&stream_2d(&self.main_index())),
            3 => flat(&stream_3d(&self.main_index())),
        }
    }
}
