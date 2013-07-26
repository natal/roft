use std::num::Zero;
use nalgebra::traits::scalar_op::ScalarMul;
use nalgebra::traits::dot::Dot;
use nalgebra::traits::norm::Norm;
use rs2cl::kernel::Kernel;
use rs2cl::nalgebra2cl::CLVec3f64;
use rs2cl::pragma;
use rs2cl::expr;
use rs2cl::cl_logic::ClOrd;

pub fn integration_kernel() -> ~str
{
  let k = @mut Kernel::new(~"integrate");

  k.enable_extension(pragma::cl_khr_fp64);

  let velocities   = k.param::<~[CLVec3f64]>(expr::Global);
  let positions    = k.param::<~[CLVec3f64]>(expr::Global);
  let invmasses    = k.param::<~[f64]>(expr::Global);
  let fext         = k.param::<CLVec3f64>(expr::Const);
  let dt           = k.param::<f64>(expr::Const);
  let num_elements = k.param::<i32>(expr::Const);

  let id = k.var::<i32>();

  id.assign(k.get_global_id(0));

  do k.if_(id.cl_lt(&num_elements))
  {
    do k.if_(invmasses[id].cl_gt(&expr::literal(0.0)))
    {
      velocities[id].assign(velocities[id] + fext.scalar_mul(&dt));
      positions[id].assign(positions[id] + velocities[id].scalar_mul(&dt));
    }
  }

  k.to_str()
}

pub fn init_constraints_kernel() -> ~str
{
  let k = @mut Kernel::new(~"init_constraints");

  k.enable_extension(pragma::cl_khr_fp64);

  let dt           = k.param::<f64>(expr::Const);
  let num_elements = k.param::<i32>(expr::Const);
  let id1s         = k.named_param::<~[i32]>(~"id1s", expr::Global);
  let id2s         = k.named_param::<~[i32]>(~"id2s", expr::Global);
  let velocities   = k.param::<~[CLVec3f64]>(expr::Global);
  let positions    = k.param::<~[CLVec3f64]>(expr::Global);
  let normals      = k.named_param::<~[CLVec3f64]>(~"normals", expr::Global);
  let fext         = k.param::<CLVec3f64>(expr::Const);
  let invmasses    = k.named_param::<~[f64]>(~"invmasses", expr::Global);
  let objectives   = k.named_param::<~[f64]>(~"objectives", expr::Global);
  let rests        = k.named_param::<~[f64]>(~"rests", expr::Global);
  let stiffs       = k.named_param::<~[f64]>(~"stiffs", expr::Global);

  let id = k.var::<i32>();

  id.assign(k.get_global_id(0));

  do k.if_(id.cl_lt(&num_elements))
  {
    let id1 = id1s[id];
    let id2 = id2s[id];

    let normal = k.var::<CLVec3f64>();

    normal.assign(positions[id1] - positions[id2]);

    let length = k.var::<f64>();

    length.assign(normal.norm());
    normal.assign(normal.normalized());

    let dvel = k.var::<f64>();

    dvel.assign(dt * ((length - rests[id]) * stiffs[id]));

    do k.if_(invmasses[id2].cl_gt(&expr::literal(0.0)))
    { dvel.assign(dvel - (velocities[id2] + fext.scalar_mul(&dt)).dot(&normal)); }

    do k.if_(invmasses[id1].cl_gt(&expr::literal(0.0)))
    { dvel.assign(dvel + (velocities[id1] + fext.scalar_mul(&dt)).dot(&normal)); }

    normals[id].assign(normal);
    objectives[id].assign(dvel);
  }

  k.to_str()
}

pub fn lin_pgs_solver_kernel() -> ~str
{
  let k = @mut Kernel::new(~"lin_pgs_solve");

  k.enable_extension(pragma::cl_khr_fp64);

  /*
   * Params
   */
  let _           = k.named_param::<i32>(~"num", expr::Const);
  let id1s        = k.named_param::<~[i32]>(~"id1s", expr::Global);
  let id2s        = k.named_param::<~[i32]>(~"id2s", expr::Global);
  let normals     = k.named_param::<~[CLVec3f64]>(~"normals", expr::Global);
  let inv_masses  = k.named_param::<~[f64]>(~"inv_masses", expr::Global);
  let impulses    = k.named_param::<~[f64]>(~"impulses", expr::Global);
  let lobounds    = k.named_param::<~[f64]>(~"lobounds", expr::Global);
  let hibounds    = k.named_param::<~[f64]>(~"hibounds", expr::Global);
  let objectives  = k.named_param::<~[f64]>(~"objectives", expr::Global);
  let pmasses     = k.named_param::<~[f64]>(~"pmasses", expr::Global);
  let MJLambdas   = k.named_param::<~[CLVec3f64]>(~"MJLambdas", expr::Global);
  let colors      = k.named_param::<~[i32]>(~"colors", expr::Global);
  let _           = k.named_param::<~[i32]>(~"batches", expr::Global);
  let batch_sizes = k.named_param::<~[i32]>(~"batch_sizes", expr::Global);
  let curr_color  = k.named_param::<i32>(~"curr_color", expr::Const);

  let id = k.var::<i32>();

  id.assign(k.get_global_id(0));

  do k.iterate(expr::literal(0), batch_sizes[id]) |_i|
  {
    let i          = k.var::<i32>();
    let d_lambda_i = k.named_var::<f64>(~"d_lambda_i");
    let id1        = k.named_var::<i32>(~"id1");
    let id2        = k.named_var::<i32>(~"id2");

    i.assign(_i + colors[curr_color] + id); // batches[id]); // FIXME: works only because there is one elem per batch
    id1.assign(id1s[i]);
    id2.assign(id2s[i]);

    /*
     * The solver itself
     */
    d_lambda_i.assign(objectives[i]);

    do k.if_(id1.cl_ge(&Zero::zero()))
    { d_lambda_i.assign(d_lambda_i + normals[i].dot(&MJLambdas[id1])); }

    do k.if_(id2.cl_ge(&Zero::zero()))
    { d_lambda_i.assign(d_lambda_i - normals[i].dot(&MJLambdas[id2])); }

    d_lambda_i.assign(d_lambda_i / pmasses[i]);

    let lambda_i_0 = k.var::<f64>();

    lambda_i_0.assign(impulses[i]);

    impulses[i].assign((lambda_i_0 + d_lambda_i).clamp(&lobounds[i], &hibounds[i]));

    d_lambda_i.assign(impulses[i] - lambda_i_0);

    do k.if_(id1.cl_ge(&Zero::zero()))
    { MJLambdas[id1].assign(MJLambdas[id1] - normals[i].scalar_mul(&(inv_masses[id1] * d_lambda_i))); }

    do k.if_(id2.cl_ge(&Zero::zero()))
    { MJLambdas[id2].assign(MJLambdas[id2] + normals[i].scalar_mul(&(inv_masses[id2] * d_lambda_i))); }
  }

  k.to_str()
}
