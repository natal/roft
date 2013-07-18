use std::vec;
use nalgebra::vec::Vec3;
use kiss3d::object::{VerticesNormalsTriangles, Object};
use graph::{Mesh, Graph};

pub fn cg2ids(graph: &Graph) -> (~[Vec3<f64>],
                                 ~[i32],
                                 ~[i32],
                                 ~[i32],
                                 ~[i32],
                                 ~[i32])
{
  let cgs = graph.export_batches();

  let mut vertices: ~[Vec3<f64>];
  let mut ids1: ~[i32] = ~[];
  let mut ids2: ~[i32] = ~[];
  let mut colors: ~[i32] = ~[];
  let mut batches: ~[i32] = ~[];
  let mut batch_sizes: ~[i32] = ~[];

  let mut indices: ~[(uint, uint)];

  (vertices, indices) = graph.export();



  let mut i: i32 = 0;
  for cgs.iter().advance |cg|
  {
    colors.push(i);
    let mut j: i32 = 0;

    for cg.batches.iter().advance |batch|
    {
      batches.push(j);
      batch_sizes.push(batch.edges.len() as i32);
      for batch.edges.iter().advance |e|
      {
        ids1.push(e.node_1.id as i32);
        ids2.push(e.node_2.id as i32);
        i = i + 1;
        j = j + 1;
      }
    }
  }
  (vertices, ids1, ids2, colors, batches, batch_sizes)
}

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
      graph.build_blob_graph(0, 0);
      graph.color_blob_graph();

      let color_groups = graph.export_batches();
      let (v1, v2, v3, v4, v5 ,v6) = cg2ids(&graph);

      println("Size color groups : " + color_groups.len().to_str());




      graph.write_blob_graph(~"./blob.dot");
      graph.write_line_graph(~"./line.dot");


      let (mvs, mvi) = graph.export();

      let mut invmasses = vec::from_elem(mvs.len(), 1.0f64);
      invmasses[0] = 0.0;
      invmasses[w] = 0.0;
      invmasses[mvs.len() - 1]     = 0.0;
      // invmasses[mvs.len() / 2 + w / 2]     = 0.0;
      invmasses[mvs.len() - w - 1] = 0.0;
      // for uint::iterate(mvs.len() - w - 1, mvs.len()) |i|
      // { invmasses[i] = 0.0 }

      let stiffness = vec::from_elem(mvi.len(), 50.0f64);

      (mvs, mvi, invmasses, stiffness)
    },
    _ => fail!("Unable to build the soft body without geometric informations.")
  }
}
