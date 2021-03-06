use std::rand::random;
use extra::time;
use nalgebra::vec::Vec3;
use kiss3d::window;
use kiss3d::camera;
use soft_body::SoftBody;
use builder;

#[main]
fn main()
{
  do window::Window::spawn("Soft body demo.") |w|
  {
    let hsub = 75;
    let quad = w.add_quad(100.0, 100.0, hsub, 75).set_color(random(), random(), random());

    let (vertices, ids1, ids2, _, _, _, _, invmasses, stiffness) =
      builder::soft_body_parameters(quad, hsub, false);
    let soft_body = @mut SoftBody::from_mesh(vertices, ids1, ids2, invmasses, stiffness);

    let timestep  = 0.016;

    do w.set_loop_callback
    {
      let before = time::precise_time_s();

      soft_body.integrate(&timestep, &Vec3::new(0.0f64, 0.0, -9.81));

      soft_body.solve(timestep.clone());

      // FIXME: solve constraints
      do quad.modify_vertices |vs|
      {
        for (v, p) in vs.mut_iter().zip(soft_body.points.iter())
        {
          *v = Vec3::new(p.position.x as f32,
                         p.position.y as f32,
                         p.position.z as f32);
        }

        true
      }

      println((1.0 / (time::precise_time_s() - before)).to_str() + " fps");
    }

    do w.camera().change_mode |m|
    {
      match *m
      {
        camera::ArcBall(ref mut ab) => {
          ab.dist      = 500.0;
          ab.dist_step = 100.0;
        },

        _ => { }
      }
    }

    w.set_light(window::StickToCamera);
  }
}
