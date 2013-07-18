use std::uint;
use std::vec;
use std::num::Zero;
use OpenCL::hl::*;
use OpenCL::vector::Vector;
use nalgebra::traits::norm::Norm;
use nalgebra::traits::scalar_op::ScalarMul;
use nalgebra::traits::dot::Dot;

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

  // point masses
  positions:  ~[CLVec3f64],
  velocities: ~[CLVec3f64],
  masses:     ~[f64],

  // constants for the solver
  id1s:      ~[i32],
  id2s:      ~[i32],
  real_id1s: ~[i32],
  real_id2s: ~[i32],
  pmasses:   ~[f64],
  impulses:  ~[f64],

  // to compute constraints
  rests:    ~[f64],
  stiffs:   ~[f64]
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

    // init constraints parameters
    let mut id1s:      ~[i32] = ~[];
    let mut id2s:      ~[i32] = ~[];
    let mut real_id1s: ~[i32] = ~[];
    let mut real_id2s: ~[i32] = ~[];
    let mut pmasses:   ~[f64] = ~[];
    let mut rests:     ~[f64] = ~[];
    let mut stiffs:    ~[f64] = ~[];

    for edges.iter().zip(stiffness.iter()).advance |(&(v1, v2), s)|
    {
      rests.push((vbuf[v1] - vbuf[v2]).norm());
      id1s.push(if invmasses[v1] == 0.0 { -1 } else { v1 as i32 });
      id2s.push(if invmasses[v2] == 0.0 { -1 } else { v2 as i32 });
      real_id1s.push(v1 as i32);
      real_id2s.push(v2 as i32);
      stiffs.push(s.clone());
      pmasses.push(invmasses[v1] + invmasses[v2]);
    }

    SoftBodyGpu {
      ext_forces: Zero::zero(),
      positions:  vbuf,
      velocities: vec::from_elem(invmasses.len(), Zero::zero()),
      masses:     invmasses,
      id1s:       id1s,
      id2s:       id2s,
      real_id1s:  real_id1s,
      real_id2s:  real_id2s,
      pmasses:    pmasses,
      rests:      rests,
      stiffs:     stiffs,
      impulses:   vec::from_elem(edges.len(), Zero::zero())
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
    self.ext_forces = fext.clone();

    // FIXME: dont re-create the buffers at each frame!
    let cl_pos = Vector::from_vec(ctx, self.positions);
    let cl_vel = Vector::from_vec(ctx, self.velocities);
    let cl_mas = Vector::from_vec(ctx, self.masses);

    integrator.set_arg(0, &cl_vel);
    integrator.set_arg(1, &cl_pos);
    integrator.set_arg(2, &cl_mas);
    integrator.set_arg(3, fext);
    integrator.set_arg(4, dt);
    integrator.set_arg(5, &(self.positions.len() as i32));

    let work_group_size = 64;
    let num_work_items  =
      work_group_size * ((self.positions.len() + (work_group_size - 1)) / work_group_size);

    enqueue_nd_range_kernel(
      &ctx.q,
      integrator,
      1,
      0,
      num_work_items  as int,
      work_group_size as int);

    self.velocities = cl_vel.to_vec();
    self.positions  = cl_pos.to_vec();
  }

  pub fn solve_gpu(&mut self, dt: &f64, solver: &Kernel, ctx: @ComputeContext)
  {
    let cl_mas = Vector::from_vec(ctx, self.masses);
    let cl_pma = Vector::from_vec(ctx, self.pmasses);
    let cl_id1 = Vector::from_vec(ctx, self.id1s);
    let cl_id2 = Vector::from_vec(ctx, self.id2s);
    let cl_imp = Vector::from_vec(ctx, self.impulses);

    let MJLambdas: ~[CLVec3f64] = vec::from_elem(self.masses.len(), Zero::zero());
    let cl_mjl = Vector::from_vec(ctx, MJLambdas);


    // XXX: initialize the constraints parameters: normals, impulses, objectives and MJlambdas
    // in parallel
    let mut low = ~[];
    let mut hig = ~[];
    let mut nor = ~[];
    let mut obj = ~[];

    for uint::iterate(0u, self.id1s.len()) |i|
    {
      let id1 = self.real_id1s[i];
      let id2 = self.real_id2s[i];

      let mut normal = self.positions[id1] - self.positions[id2];
      let     length = normal.normalize();
      let     limit  = self.stiffs[i] * (length - self.rests[i]).abs();


      let mut dvel = 0.0;

      if !self.masses[id2].is_zero()
      { dvel = dvel - (self.velocities[id2] + self.ext_forces.scalar_mul(dt)).dot(&normal) }

      if !self.masses[id1].is_zero()
      { dvel = dvel + (self.velocities[id1] + self.ext_forces.scalar_mul(dt)).dot(&normal) }

      low.push(-limit);
      hig.push(limit);

      nor.push(normal);

      obj.push(dvel);
    }

    let cl_low = Vector::from_vec(ctx, low);
    let cl_hig = Vector::from_vec(ctx, hig);
    let cl_nor = Vector::from_vec(ctx, nor);
    let cl_obj = Vector::from_vec(ctx, obj);

    // FIXME: there are a lot of arguments which dont change at each iteration.
    // So there might be no need to re-upload them!
    // 0 : niter
    // 1 : num
    // 2 : id1s
    // 3 : id2s
    // 4 : normals
    // 5 : inv_masses
    // 6 : impulses
    // 7 : lobounds
    // 8 : hibounds
    // 9 : objectives
    // 10: pmasses
    // 11: MJLambdas
    solver.set_arg(0, &50i32);
    solver.set_arg(1, &(self.pmasses.len() as i32));
    solver.set_arg(2, &cl_id1);
    solver.set_arg(3, &cl_id2);
    solver.set_arg(4, &cl_nor);
    solver.set_arg(5, &cl_mas);
    solver.set_arg(6, &cl_imp);
    solver.set_arg(7, &cl_low);
    solver.set_arg(8, &cl_hig);
    solver.set_arg(9, &cl_obj);
    solver.set_arg(10, &cl_pma);
    solver.set_arg(11, &cl_mjl);

    enqueue_nd_range_kernel(&ctx.q, solver, 1, 0, 1, 1);

    let dvs = cl_mjl.to_vec();

    for self.velocities.mut_iter().zip(dvs.iter()).advance |(v, dv)|
    { *v = *v + *dv }

    self.impulses = cl_imp.to_vec();
  }
}
