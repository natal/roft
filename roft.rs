use std::uint;
use std::vec;
use std::rand::random;
use nalgebra::vec::Vec3;
use kiss3d::window;
use soft_body::SoftBody;
use builder;
use kiss3d::object::{VerticesNormalsTriangles, Object};
use graph::{Mesh, Graph};

#[main]
fn main()
{
  do window::Window::spawn(~"Soft body demo.") |w|
  {
    let hsub = 50;
    let quad = w.add_quad(20.0, 20.0, hsub, 50).set_color(random(), random(), random());

    let (mvs, mvi, invmasses, stiffness) = builder::soft_body_parameters(quad, hsub);
    let soft_body = @mut SoftBody::from_mesh(mvs, mvi, invmasses, stiffness);

    let timestep  = 0.016;

    do w.set_loop_callback |_|
    {
      soft_body.integrate(&timestep, &Vec3::new([ 0.0f64, 0.00f64, -9.81 ]));

      soft_body.solve(timestep.clone());

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
      println("Building blob graph");

      // second parameter is the minimum number of internal connections to
      //consider 2 blobs connected
      // graph.build_blob_graph(6, 1);
      // graph.color_blob_graph();
      // graph.write_blob_graph();
      // graph.write_line_graph();

      // Pour recuperer les blobs:
      // graph.blobs[i].color() -> recupere la couleur du blob i
      // graph.blobs[i].content.sub_nodes[j].content ->
      //       recupere l'edge j qui est dans le blob i
      // graph.blobs[i].content.sub_nodes[j].content.color() ->
      //       recupere la couleur de l'edge j dans le blob i
      // graph.blobs[i].content.sub_nodes[j].content.node_1 ->
      //       recupere le noeud 1 de l'edge j dans le blob i
      // graph.blobs[i].content.sub_nodes[j].content.node_2 ->
      //       recupere le noeud 2 de l'edge j dans le blob i
      // graph.blobs[i].content.sub_nodes[j].content.node_1.pos ->
      //       recupere la position du noeud 1 de l'edge j dans le blob i
      //


      let (mvs, mvi) = graph.export();

      let mut invmasses = vec::from_elem(mvs.len(), 1.0f64);
      // invmasses[0] = 0.0;
      // invmasses[w] = 0.0;
      for uint::iterate(mvs.len() - w - 1, mvs.len()) |i|
      { invmasses[i] = 0.0 }

      let stiffness = vec::from_elem(mvi.len(), 1000f64);

      SoftBody::from_mesh(mvs, mvi, invmasses, stiffness)
    },
    _ => fail!("Unable to build the soft body without geometric informations.")
  }
}
