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
  real_id1s:   ~[i32],
  cl_real_id1: Vector<i32>,
  real_id2s:   ~[i32],
  cl_real_id2: Vector<i32>,

  cl_id1:    Vector<i32>,
  cl_id2:    Vector<i32>,

  pmasses:   ~[f64],

  impulses:  ~[f64],
  cl_imp:    Vector<f64>,

  low:       ~[f64],
  hig:       ~[f64],

  // to compute constraints
  rests:    ~[f64],
  cl_rest:  Vector<f64>,
  stiffs:   ~[f64],
  cl_stiff: Vector<f64>,

  // cl buffers
  normals:  ~[CLVec3f64],
  cl_nor:   Vector<CLVec3f64>,

  objectives: ~[f64],
  cl_obj:     Vector<f64>,

  cl_mas:     Vector<f64>,
  cl_pma:     Vector<f64>,
  cl_low:     Vector<f64>,
  cl_hig:     Vector<f64>,

  cl_pos:     Vector<CLVec3f64>,
  cl_vel:     Vector<CLVec3f64>,

}

impl SoftBodyGpu
{
  pub fn from_mesh(vbuf:      ~[CLVec3f64],
                   edges:     ~[(uint, uint)],
                   invmasses: ~[f64],
                   stiffness: ~[f64],
                   solver:    &Kernel,
                   ctx:       @ComputeContext) -> SoftBodyGpu
  {
    assert!(vbuf.len() == invmasses.len(),
            "Vertex buffer and mass informations must have the same size.");

    assert!(stiffness.len() == edges.len(),
            "Edge buffer and stiffness informations must have the same size.");

    // init constraints parameters
    let mut cl_id1s:   ~[i32] = ~[];
    let mut cl_id2s:   ~[i32] = ~[];
    let mut id1s:      ~[i32] = ~[];
    let mut id2s:      ~[i32] = ~[];
    let mut pmasses:   ~[f64] = ~[];
    let mut rests:     ~[f64] = ~[];
    let mut stiffs:    ~[f64] = ~[];

    for edges.iter().zip(stiffness.iter()).advance |(&(v1, v2), s)|
    {
      rests.push((vbuf[v1] - vbuf[v2]).norm());
      cl_id1s.push(if invmasses[v1] == 0.0 { -1 } else { v1 as i32 });
      cl_id2s.push(if invmasses[v2] == 0.0 { -1 } else { v2 as i32 });
      id1s.push(v1 as i32);
      id2s.push(v2 as i32);
      stiffs.push(s.clone());
      pmasses.push(invmasses[v1] + invmasses[v2]);
    }

    let imps       = vec::from_elem(edges.len(), Zero::zero());
    let normals    = vec::from_elem(edges.len(), Zero::zero());
    let objectives = vec::from_elem(edges.len(), Zero::zero());
    let low        = vec::from_elem(rests.len(), -Bounded::max_value::<f64>());
    let hig        = vec::from_elem(rests.len(), Bounded::max_value::<f64>());
    let vels       = vec::from_elem(invmasses.len(), Zero::zero());

    let res = SoftBodyGpu {
      ext_forces: Zero::zero(),
      cl_pos:     Vector::from_vec(ctx, vbuf),
      positions:  vbuf,
      cl_vel:     Vector::from_vec(ctx, vels),
      velocities: vels,
      cl_mas:     Vector::from_vec(ctx, invmasses),
      masses:     invmasses,
      cl_real_id1: Vector::from_vec(ctx, id1s),
      real_id1s:   id1s,
      cl_real_id2: Vector::from_vec(ctx, id2s),
      real_id2s:   id2s,
      cl_id1:     Vector::from_vec(ctx, cl_id1s),
      cl_id2:     Vector::from_vec(ctx, cl_id2s),
      cl_pma:     Vector::from_vec(ctx, pmasses),
      pmasses:    pmasses,
      cl_stiff:   Vector::from_vec(ctx, stiffs),
      stiffs:     stiffs,
      cl_imp:     Vector::from_vec(ctx, imps),
      impulses:   imps,
      cl_nor:     Vector::from_vec(ctx, normals),
      normals:    normals,
      cl_obj:     Vector::from_vec(ctx, objectives),
      objectives: objectives,
      cl_low:     Vector::from_vec(ctx, low),
      low:        low,
      cl_hig:     Vector::from_vec(ctx, hig),
      hig:        hig,
      cl_rest:    Vector::from_vec(ctx, rests),
      rests:      rests,
    };

    solver.set_arg(0, &50i32);
    solver.set_arg(1, &(res.pmasses.len() as i32));
    solver.set_arg(2, &res.cl_id1);
    solver.set_arg(3, &res.cl_id2);
    solver.set_arg(4, &res.cl_nor);
    solver.set_arg(5, &res.cl_mas);
    solver.set_arg(6, &res.cl_imp);
    solver.set_arg(7, &res.cl_low);
    solver.set_arg(8, &res.cl_hig);
    solver.set_arg(9, &res.cl_obj);
    solver.set_arg(10, &res.cl_pma);

    res
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
    self.cl_vel.rewrite(self.velocities);
    self.cl_pos.rewrite(self.positions);

    integrator.set_arg(0, &self.cl_vel);
    integrator.set_arg(1, &self.cl_pos);
    integrator.set_arg(2, &self.cl_mas);
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

    self.cl_vel.to_existing_vec(self.velocities);
    self.cl_pos.to_existing_vec(self.positions);
  }

  pub fn solve_gpu(&mut self, dt: &f64, solver: &Kernel, initializer: &Kernel, ctx: @ComputeContext)
  {
    let MJLambdas: ~[CLVec3f64] = vec::from_elem(self.masses.len(), Zero::zero());

    let cl_mjl = Vector::from_vec(ctx, MJLambdas);


    // 0 : dt
    // 1 : num_elements
    // 2 : id1s
    // 3 : id2s
    // 4 : velocities
    // 5 : positions
    // 6 : normals
    // 7 : fext
    // 8 : invmasses
    // 9 : objectives
    // 10: rests
    // 11: stiffs
    initializer.set_arg(0,  dt);
    initializer.set_arg(1,  &(self.pmasses.len() as i32));
    initializer.set_arg(2,  &self.cl_real_id1);
    initializer.set_arg(3,  &self.cl_real_id2);
    initializer.set_arg(4,  &self.cl_vel);
    initializer.set_arg(5,  &self.cl_pos);
    initializer.set_arg(6,  &self.cl_nor);
    initializer.set_arg(7,  &self.ext_forces);
    initializer.set_arg(8,  &self.cl_mas);
    initializer.set_arg(9,  &self.cl_obj);
    initializer.set_arg(10, &self.cl_rest);
    initializer.set_arg(11, &self.cl_stiff);

    let work_group_size = 64;
    let num_work_items  =
      work_group_size * ((self.pmasses.len() + (work_group_size - 1)) / work_group_size);

    enqueue_nd_range_kernel(
      &ctx.q,
      initializer,
      1,
      0,
      num_work_items  as int,
      work_group_size as int);

    self.cl_imp.rewrite(self.impulses);

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
    /*
     * solver.set_arg(0, &50i32);
     * solver.set_arg(1, &(self.pmasses.len() as i32));
     * solver.set_arg(2, &self.cl_id1);
     * solver.set_arg(3, &self.cl_id2);
     * solver.set_arg(4, &self.cl_nor);
     * solver.set_arg(5, &self.cl_mas);
     * solver.set_arg(6, &self.cl_imp);
     * solver.set_arg(7, &self.cl_low);
     * solver.set_arg(8, &self.cl_hig);
     * solver.set_arg(9, &self.cl_obj);
     * solver.set_arg(10, &self.cl_pma);
     */
    solver.set_arg(11, &cl_mjl);

    enqueue_nd_range_kernel(&ctx.q, solver, 1, 0, 1, 1);

    let dvs = cl_mjl.to_vec();

    for self.velocities.mut_iter().zip(dvs.iter()).advance |(v, dv)|
    { *v = *v + *dv }

    self.cl_imp.to_existing_vec(self.impulses);
  }
}
