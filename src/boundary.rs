use arrayfire::*;
use std::collections::HashMap;
use crate::grid::StructuredGrid;
use crate::traits::{Geometry, Distribution};
use crate::geometry::Circle;
use crate::FloatNum;

#[derive(Copy, Clone, PartialEq, PartialOrd)]
pub enum Type {
    BounceBack,
    Inflow(FloatNum, FloatNum),
}

pub trait AnyCondition: Send + Sync {
    #[inline(always)]
    fn condition(&self) -> Type;
    #[inline(always)]
    fn contains(&self, x: grid::X) -> bool;
}

pub struct Condition<T: Geometry + Send + Sync> {
    condition: Type,
    geometry: T,
}

impl<T: Geometry + Send + Sync> Condition<T> {
    pub fn new(c: Type, g: T) -> Condition<T> {
        Condition {
            condition: c,
            geometry: g,
        }
    }
}

impl<T: Geometry + Send + Sync> AnyCondition for Condition<T> {
    #[inline(always)]
    fn condition(&self) -> Type {
        self.condition
    }
    #[inline(always)]
    fn contains(&self, x: grid::X) -> bool {
        self.geometry.contains(x)
    }
}

pub struct Handler<D: Distribution> {
    grid: StructuredGrid<D>,
    occupied_nodes: Array<FloatNum>,
    boundary_conditions: HashMap<&'static str, Box<dyn AnyCondition>>,
    to_reflect: Array<u32>,
    reflected: Array<u32>,
}

impl<D: Distribution> Handler<D> {
    pub fn new(grid: StructuredGrid<D>) -> Self {
        Self {
            grid,
            boundary_conditions: HashMap::default(),
            occupied_nodes: constant::<FloatNum>(0.0, grid.dimensions),
            to_reflect: constant::<u32>(0, grid.dimensions),
            reflected: constant::<u32>(0, grid.dimensions)
        }
    }
    pub fn add(&mut self, label: &'static str, bc: Box<dyn AnyCondition>) {
        self.boundary_conditions.insert(label, bc);
    }

    pub fn update_bounceback_indices(&mut self) {
        // matrix offset of each Occupied Node
        let on = locate(&self.occupied_nodes);

        // Bounceback indexes
        let ci = D::index();
        let nbi = D::opposite_index();

        self.to_reflect = flat(&tile(&on, dim4!(ci.elements() as u64)))
            + flat(&tile(&ci, dim4!(on.elements() as u64)));
        self.reflected = flat(&tile(&on, dim4!(nbi.elements() as u64)))
            + flat(&tile(&nbi, dim4!(on.elements() as u64)));
    }

    #[inline(always)]
    pub fn solid_boundary(&self, x: grid::X) -> bool {
        for bc in &self.boundary_conditions {
            if bc.contains(x) && bc.condition() == Type::BounceBack {
                return true;
            }
        }
        false
    }

    #[inline(always)]
    pub fn idx(&self, x: grid::X) -> Option<usize> {
        for (idx, bc) in self.boundary_conditions.iter().enumerate() {
            if bc.contains(x) {
                return Some(idx);
            }
        }
        None
    }

    #[inline(always)]
    pub fn apply<F, H>(
        &self,
        f: &Array<FloatNum>,
        f_hlp: &Array<FloatNum>,
    ) -> Array<FloatNum> {
        let mut r = f.clone();

        for bc in &self.boundary_conditions {
            match bc.condition() {
                Type::BounceBack => {
                    let mut idxrs = Indexer::default();
                    idxrs.set_index(&self.to_reflect, 0, None);
                    let bouncedback = index_gen(self.f, idxrs);
                    eval!(f[self.reflected] = bouncedback);
                }
                Type::Inflow(density, accel) => {
                   
                }
            }
        }
        r
    }
}