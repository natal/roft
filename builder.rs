use std::vec;
use nalgebra::vec::Vec3;
use kiss3d::object::{VerticesNormalsTriangles, Object};
use graph::{Mesh, Graph};

pub fn soft_body_parameters(quad: @mut Object, w: uint) -> (~[Vec3<f64>], ~[(uint, uint)], ~[f64], ~[f64])
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
      graph.build_blob_graph(8, 2);
      graph.color_blob_graph();
      graph.write_blob_graph();
      graph.write_line_graph();


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

      (mvs, mvi, invmasses, stiffness)
    },
    _ => fail!("Unable to build the soft body without geometric informations.")
  }
}
