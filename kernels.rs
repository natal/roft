use nalgebra::traits::scalar_op::ScalarMul;
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
  let fext         = k.param::<CLVec3f64>(expr::Const);
  let dt           = k.param::<f64>(expr::Const);
  let num_elements = k.param::<u32>(expr::Const);

  let id         = k.var::<u32>();

  id.assign(k.get_global_id(0));

  do k.if_(id.cl_lt(&num_elements))
  {
    velocities[id].assign(velocities[id] + fext.scalar_mul(&dt));
    positions[id].assign(positions[id] + velocities[id].scalar_mul(&dt));
  }

  k.to_str()
}
