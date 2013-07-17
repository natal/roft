use std::vec;
use std::num::Zero;
use OpenCL::hl::*;
use OpenCL::vector::Vector;
use nalgebra::traits::norm::Norm;

use rs2cl::nalgebra2cl::CLVec3f64;

pub struct ConstraintsGeometry
{
  stiffness:   f64,
  rest_length: f64,
  impulse:     f64,
  rb1:         uint,
  rb2:         uint
}

pub struct SoftBodyGpu
{
  ext_forces:  CLVec3f64,
  constraints: ~[ConstraintsGeometry],

  // point masses
  positions:  ~[CLVec3f64],
  velocities: ~[CLVec3f64],
  masses:     ~[f64],
}

impl SoftBodyGpu
{
  pub fn from_mesh(vbuf:      ~[CLVec3f64],
                   edges:     ~[(uint, uint)],
                   invmasses: ~[f64],
                   stiffness: ~[f64]) -> SoftBodyGpu
  {
    assert!(vbuf.len() == invmasses.len(),
            "Vertex buffer and mass informations must have the same size.");

    assert!(stiffness.len() == edges.len(),
            "Edge buffer and stiffness informations must have the same size.");

    // create constraints
    let mut constraints = ~[];

    for edges.iter().zip(stiffness.iter()).advance |(&(v1, v2), s)|
    {
      constraints.push(ConstraintsGeometry {
        stiffness:   s.clone(),
        rest_length: (vbuf[v1] - vbuf[v2]).norm(),
        impulse:     Zero::zero(),
        rb1:         v1,
        rb2:         v2
      });
    }

    SoftBodyGpu {
      constraints: constraints,
      ext_forces:  Zero::zero(),
      positions:   vbuf,
      velocities:  vec::from_elem(invmasses.len(), Zero::zero()),
      masses:      invmasses
    }
  }
}

impl SoftBodyGpu
{
  pub fn integrate_gpu(&mut self,
                       dt:         &f64,
                       fext:       &CLVec3f64,
                       integrator: &Kernel,
                       ctx:        @ComputeContext)
  {
    let cl_pos = Vector::from_vec(ctx, self.positions);
    let cl_vel = Vector::from_vec(ctx, self.velocities);

    integrator.set_arg(0, &cl_vel);
    integrator.set_arg(1, &cl_pos);
    integrator.set_arg(2, fext);
    integrator.set_arg(3, dt);
    integrator.set_arg(4, &(self.positions.len() as u32));

    enqueue_nd_range_kernel(
      &ctx.q,
      integrator,
      1,
      0,
      self.positions.len() as int,
      2);

    self.velocities = cl_vel.to_vec();
    self.positions  = cl_pos.to_vec();
  }
}

impl SoftBodyGpu
{
  pub fn solve_gpu(&mut self, _: f64)
  {
    fail!("Not yet implemented.");
  }
}
