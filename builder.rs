use std::vec;
use nalgebra::vec::Vec3;
use kiss3d::object::{VerticesNormalsTriangles, Object};
use graph::{Mesh, Graph};

pub fn cg2ids(graph: &mut Graph) -> (~[Vec3<f64>],
                                 ~[i32],
                                 ~[i32],
                                 ~[i32],
                                 ~[i32],
                                 ~[i32],
                                 ~[i32])
{
  let cgs = graph.export_batches();

  let mut ids1: ~[i32] = ~[];
  let mut ids2: ~[i32] = ~[];
  let mut colors: ~[i32] = ~[];
  let mut colors_sizes: ~[i32] = ~[];
  let mut batches: ~[i32] = ~[];
  let mut batch_sizes: ~[i32] = ~[];

  let (vertices, _, _) = graph.export();



  let mut i: i32 = 0;
  for cgs.iter().advance |cg|
  {
    colors.push(i);
    colors_sizes.push(cg.batches.len() as i32);
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
  (vertices, ids1, ids2, colors, colors_sizes, batches, batch_sizes)
}

pub fn cg2ids_no_blob(graph: &mut Graph) -> (~[Vec3<f64>],
                                             ~[i32],
                                             ~[i32],
                                             ~[i32],
                                             ~[i32],
                                             ~[i32],
                                             ~[i32])
{
  let cgs = graph.export_edges();

  let mut ids1: ~[i32] = ~[];
  let mut ids2: ~[i32] = ~[];
  let mut colors: ~[i32] = ~[];
  let mut colors_sizes: ~[i32] = ~[];
  let mut batches: ~[i32] = ~[];
  let mut batch_sizes: ~[i32] = ~[];

  let (vertices, _, _) = graph.export();



  let mut i: i32 = 0;
  for cgs.iter().advance |cg|
  {
    colors.push(i);
    colors_sizes.push(cg.len() as i32);
    let mut j: i32 = 0;

    for cg.iter().advance |e|
    {
      batches.push(j);
      batch_sizes.push(1);
      ids1.push(e.node_1.id as i32);
      ids2.push(e.node_2.id as i32);
      i = i + 1;
      j = j + 1;
    }
  }
  (vertices, ids1, ids2, colors, colors_sizes, batches, batch_sizes)
}

pub fn soft_body_parameters(quad: @mut Object, w: uint, color_graph: bool) -> (~[Vec3<f64>], ~[i32], ~[i32], ~[i32], ~[i32], ~[i32], ~[i32], ~[f64], ~[f64])
{
  match quad.geometry()
  {
    &VerticesNormalsTriangles(ref vs, _, ref ts) =>
    {
      let     mesh  = Mesh::new(vs.clone(), ts.clone());
      let mut graph = Graph::new(mesh);

      graph.augment();
      graph.build_edge_graph();


      let (vertices, ids1, ids2, colors, colors_sizes, batches, batch_sizes) =
      if color_graph
      {
        println("Preprocessing, please wait...");
        //graph.build_blob_graph(0, 0);
        graph.color_edge_graph();
        cg2ids_no_blob(&mut graph)
      }

      else
      {
        let (mvs, ids1_cpu, ids2_cpu) = graph.export();
        (mvs, ids1_cpu, ids2_cpu, ~[], ~[], ~[], ~[])
      };

      let mut invmasses = vec::from_elem(vertices.len(), 1.0f64);
      invmasses[0] = 0.0;
      invmasses[w] = 0.0;
      invmasses[vertices.len() - 1]     = 0.0;
      // invmasses[mvs.len() / 2 + w / 2]     = 0.0;
      invmasses[vertices.len() - w - 1] = 0.0;
      // for uint::iterate(mvs.len() - w - 1, mvs.len()) |i|
      // { invmasses[i] = 0.0 }

      let stiffness = vec::from_elem(ids1.len(), 50.0f64);

      (vertices, ids1, ids2, colors, colors_sizes, batches, batch_sizes, invmasses, stiffness)
    },
    _ => fail!("Unable to build the soft body without geometric informations.")
  }
}
