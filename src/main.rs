pub mod lbm;
pub mod grid;
pub mod boundary;
pub mod distribution;
pub mod geometry;
pub mod physics;
pub mod traits;
pub mod streaming;
pub mod solver;
pub mod lib;

use arrayfire::*;
use self::solver::Solver;
use self::lib::FloatNum;

// Configure the numerical method
type Distribution = distribution::D2Q9;
type Collision = physics::ns::SingleRelaxationTime;
type Physics = physics::NavierStokes<Distribution, Collision>;

fn main() {
    let distribution = distribution::d2q9::D2Q9::new();
    // Initialize the grid, physical parameters, and solver:
    let grid = grid::StructuredGrid::<Distribution>::new(300, 150, 0);

    // Add initial static boundary conditions
    let bcs = boundary::Handler::<Distribution>::new(grid);
    {
        // Cylinder:
        {
            let cyl = Box::new(boundary::Condition::new(
                boundary::Type::BounceBack,
                geometry::Circle::new(grid.x, grid.y),
            ));
            bcs.add("cylinder", cyl);
        }
        // Bottom channel wall:
        {
            let bottom_wall = Box::new(boundary::Condition::new(
                boundary::Type::BounceBack,
                geometry::Plane::new((0, 1), (0, 0)),
            ));
            bcs.add("bottom_wall", bottom_wall);
        }
        // Top channel wall:
        {
            let top_wall = Box::new(boundary::Condition::new(
                boundary::Type::BounceBack,
                geometry::Plane::new((0, -1), (0, grid.y - 1)),
            ));
            bcs.add("top_wall", top_wall);
        }
    }

    // Physical parameters.
    let rho0: FloatNum = 1.0;
    // Reynolds number
    let re: FloatNum = 220.0;
    // Lattice speed
    let u_accel: FloatNum = 0.1;

    let obstacle_x: u64 = grid.x / 5 + 1; // x location of the cylinder
    let obstacle_y: u64 = grid.y / 2 + grid.y / 30; // y location of the cylinder
    let obstacle_r: u64 = grid.y / 10 + 1; // radius of the cylinder

    // Kinematic viscosity
    let nu: FloatNum = u_accel * 2.0 * obstacle_r as FloatNum / re;
    // Relaxation time
    let tau: FloatNum = (3.0 as FloatNum) * nu + (0.5 as FloatNum);
    // Relaxation parameter
    let omega: FloatNum = (1.0 as FloatNum) / tau;

    // Initialize physics
    let physics: Physics = Physics::new(rho0, u_accel, grid.dimensions, Collision { omega, re, nu, tau });

    // Initialize physical boundary conditions
    {
        // Periodic forced inflow:
        {
            let bc = Box::new(boundary::Condition::new(
                boundary::Type::Inflow(
                    physics.inflow_density,
                    physics.inflow_accel,
                ),
                geometry::Plane::new((1, 0), (0, 0)),
            ));
            bcs.add("left_wall_inflow", bc);
        }
    }

    let mut s = Solver::<Physics, Distribution>::new(distribution, grid, bcs, physics);

    // Initialize distribution functions
    s.initialize();

    // Simulation loop
    let mut n_it = 1000;
    let n_out = 100;
    assert!(n_it > 0);
    let mut iter = 0;

    loop {
        let write_output = n_out > 0 && iter % n_out == 0;
        s.run(iter, write_output);
        n_it -= 1;

        if n_it == 0 {
            break;
        }
    }

}