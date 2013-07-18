use std::rand::random;
use nalgebra::vec::Vec3;
use kiss3d::window;
use OpenCL::hl::*;
use rs2cl::nalgebra2cl::CLVec3f64;
use soft_body_gpu::SoftBodyGpu;
use builder;
use kernels;

#[main]
fn main()
{
  do window::Window::spawn(~"Soft body demo.") |w|
  {
    /*
     * Initialize OpenCL
     */
    let ctx  = create_compute_context();

    let src  = kernels::integration_kernel() + kernels::lin_pgs_solver_kernel();
    let prog = ctx.create_program_from_source(src);

    println(src);

    prog.build(ctx.device);

    let integrator = prog.create_kernel("integrate");
    let solver     = prog.create_kernel("lin_pgs_solve");

    /*
     * Initialize simulation parameters.
     */
    let hsub = 25;
    let quad = w.add_quad(10.0, 10.0, hsub, 25).set_color(random(), random(), random());

    let (mvs, mvi, invmasses, stiffness) = builder::soft_body_parameters(quad, hsub);

    let cl_mvs = mvs.consume_iter().transform(|v| CLVec3f64::new(v)).collect();
    let soft_body = @mut SoftBodyGpu::from_mesh(cl_mvs, mvi, invmasses, stiffness);

    let timestep: f64 = 0.016;

    /*
     * Initialize opencl
     */

    do w.set_loop_callback |_|
    {
      let gravity = CLVec3f64::new(Vec3::new([ 0.0f64, 0.00, 9.81f64 ]));
      soft_body.integrate_gpu(&timestep, &gravity, &integrator, ctx);

      soft_body.solve_gpu(&timestep, &solver, ctx);

      do quad.modify_vertices |vs|
      {
        for vs.mut_iter().zip(soft_body.positions.iter()).advance |(v, p)|
        {
          *v = Vec3::new([ p.val.at[0] as f32,
                           p.val.at[1] as f32,
                           p.val.at[2] as f32 ]);
        }

        true
      }
    }
  }
}
