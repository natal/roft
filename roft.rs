use std::vec;
use nalgebra::vec::Vec3;
use kiss3d::window;
use kiss3d::object::{VerticesNormalsTriangles, Object};
use soft_body::SoftBody;
use graph::{Mesh, Graph};

#[main]
fn main()
{
  do window::Window::spawn(~"Soft body demo.") |w|
  {
    let hsub = 25;
    let quad = w.add_quad(10.0, 10.0, hsub, 25).set_texture(~"texture.jpg");

    let soft_body = @mut quad_to_soft_body(quad, hsub);

    // w.set_wireframe_mode(true);

    do w.set_loop_callback |_|
    {
      soft_body.integrate(&0.016, &Vec3::new([ 0.0f64, 0.00, 9.81f64 ]));

      soft_body.solve(0.016);

      // FIXME: solve constraints
      do quad.modify_vertices |vs|
      {
        for vs.mut_iter().zip(soft_body.points.iter()).advance |(v, p)|
        {
          *v = Vec3::new([ p.position.at[0] as f32,
                           p.position.at[1] as f32,
                           p.position.at[2] as f32 ]);
        }

        true
      }
    }
  }
}

fn quad_to_soft_body(quad: @mut Object, w: uint) -> SoftBody<f64, Vec3<f64>>
{
  match quad.geometry()
  {
    &VerticesNormalsTriangles(ref vs, _, ref ts) =>
    {
      let     mesh  = Mesh::new(vs.clone(), ts.clone());
      let mut graph = Graph::new(mesh);

      graph.augment();
      graph.build_edge_graph();

      let (mvs, mvi) = graph.export();

      let mut invmasses = vec::from_elem(mvs.len(), 1.0f64);
      invmasses[0] = 0.0;
      invmasses[w] = 0.0;
      invmasses[mvs.len() - 1]     = 0.0;
      // invmasses[mvs.len() / 2 + w / 2]     = 0.0;
      invmasses[mvs.len() - w - 1] = 0.0;
      // for uint::iterate(mvs.len() - w - 1, mvs.len()) |i|
      // { invmasses[i] = 0.0 }

      let stiffness = vec::from_elem(mvi.len(), 1000f64);

      SoftBody::from_mesh(mvs, mvi, invmasses, stiffness)
    },
    _ => fail!("Unable to build the soft body without geometric informations.")
  }
}
