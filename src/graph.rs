extern mod extra;
extern mod nalgebra;
extern mod kiss3d;

use std::uint;
use std::vec;
use nalgebra::vec::Vec3;
use node::Node;
use edge::Edge;
use vertex::Vertex;

type Vec3f = Vec3<f32>;


pub struct Blob<T>
{
  sub_nodes: ~[@mut Node<T>]
}

impl<T> Blob<T>
{
  pub fn new() -> Blob<T>
  {
    Blob
    {
      sub_nodes: ~[]
    }
  }

  pub fn merge(&mut self, b: &Blob<T>)
  {
    for b.sub_nodes.iter().advance |sb|
    {
      self.sub_nodes.push(*sb);
    }
  }

  pub fn is_singleton(&self) -> bool
  {
    self.sub_nodes.len() == 1
  }

  pub fn nb_adj_elements(&self, b2: &Blob<T>) -> uint
  {
    let mut count = 0u;
    for self.sub_nodes.iter().advance |e1|
    {
      for b2.sub_nodes.iter().advance |e2|
      {
        if e1.is_adj_to(*e2)
        { count = count + 1 }
      }
    }
    count
  }

  pub fn disconnect_adj_elements(&self, b2: &Blob<T>)
  {
    for self.sub_nodes.iter().advance |e1|
    {
      for b2.sub_nodes.iter().advance |e2|
      {
        if e1.is_adj_to(*e2)
        { Node::disconnect(*e1, *e2) }
      }
    }
  }
}

pub struct Mesh
{
  vbuff: ~[Vec3f],
  ibuff: ~[(u32, u32, u32)]
}

impl Mesh
{
  pub fn new(vb: ~[Vec3f], ib: ~[(u32, u32, u32)]) -> Mesh
  {
    Mesh
    {
      vbuff: vb,
      ibuff: ib
    }
  }
}

#[deriving(Clone)]
pub struct Batch
{
  edges: ~[Edge]
}

impl Batch
{
  pub fn new(edges: ~[@mut Node<Edge>]) -> Batch
  {
    let mut new_batch = Batch
    {
      edges: ~[]
    };

    for edges.iter().advance |e|
    { new_batch.edges.push(e.content) }

    new_batch
  }
}

#[deriving(Clone)]
pub struct ColorGroup
{
  batches: ~[Batch]
}

impl ColorGroup
{
  pub fn new() -> ColorGroup
  {
    ColorGroup
    {
      batches: ~[]
    }
  }
}

pub struct Graph
{
  nodes:              ~[@mut Node<Vertex>],
  edges:              ~[@mut Node<Edge>],
  blobs:              ~[@mut Node<Blob<Edge>>],
  priv blob_chrom_nb: int,
  priv edge_chrom_nb: int
}

impl Graph
{
  pub fn new(mesh: Mesh) -> Graph
  {
    let mut nodes: ~[@mut Node<Vertex>] = ~[];
    for mesh.vbuff.iter().enumerate().advance |(vid, v)|
    { nodes.push(@mut Node::new(vid, Vertex::new(*v), *v)) }


    for mesh.ibuff.iter().advance |&(id1, id2, id3)|
    {
      Node::connect(nodes[id1], nodes[id2]);
      Node::connect(nodes[id1], nodes[id3]);
      Node::connect(nodes[id2], nodes[id3]);
    }

    Graph
    {
      nodes: nodes,
      edges: ~[],
      blobs: ~[],
      blob_chrom_nb: 0,
      edge_chrom_nb: 0
    }
  }

  pub fn augment(&mut self)
  {
    let mut to_connect : ~[~[@mut Node<Vertex>]] = ~[];
    for uint::iterate(0u, self.nodes.len()) |i|
    {
      self.unmark();
      let n1 = self.nodes[i];
      to_connect.push(n1.distant_nodes(|n2| (n2.dist() == 2 && n1.share_k_adjs(n2, 2)), 2));
    }

    for self.nodes.iter().enumerate().advance |(i, n)|
    {
      for to_connect[i].iter().advance |n2|
      { Node::connect(*n, *n2) }
    }
  }

  pub fn unmark(&mut self)
  {
     for self.nodes.iter().advance |n|
     { n.unmark() }
     for self.edges.iter().advance |e|
     { e.unmark() }
     for self.blobs.iter().advance |b|
     { b.unmark() }
  }

  pub fn intern_unmark(&mut self)
  {
     for self.nodes.iter().advance |n|
     { n.intern_unmark() }
     for self.edges.iter().advance |e|
     { e.intern_unmark() }
     for self.blobs.iter().advance |b|
     { b.intern_unmark() }
  }

  pub fn reset_dists(&mut self)
  {
     for self.nodes.iter().advance |n1|
     { n1.reset_dist() }
     for self.edges.iter().advance |e1|
     { e1.unmark() }
  }


  pub fn build_blob_graph(&mut self, dist: uint, min_connections: uint)
  {
    let mut c = 0u;
    self.unmark();
    for self.edges.iter().advance |e|
    {
      if !e.is_marked()
      {
        let mut blob = Blob::new();
        blob.sub_nodes.push(*e);

        e.mark();


        let dist_nodes = e.distant_nodes(|n| (n.dist() <= dist), dist);
        for dist_nodes.iter().advance |de|
        {
          if !de.is_marked()
          {
            blob.sub_nodes.push(*de);
            de.mark();
          }
        }
        self.blobs.push(@mut Node::new(c, blob, e.pos));
        c = c + 1;
      }
    }

    for self.blobs.iter().advance |b1|
    {
      for self.blobs.iter().advance |b2|
      {
        if b1.content.nb_adj_elements(&b2.content) > min_connections
        { Node::connect(*b1, *b2) }
        else
        { b1.content.disconnect_adj_elements(&b2.content) }
      }
    }
  }

  pub fn export_batches(&self) -> ~[ColorGroup]
  {
    let mut color_groups = vec::from_elem(self.blob_chrom_nb as uint, ColorGroup::new());

    for self.blobs.iter().advance |blob|
    {
      if blob.color() < 0
      { fail!("blob graph has not been colored correctly") }
      color_groups[blob.color() as uint].batches.push(Batch::new(blob.content.sub_nodes.clone()));
    }

    color_groups
  }

  pub fn export_edges(&self) -> ~[~[Edge]]
  {
    let mut color_groups = vec::from_elem(self.edge_chrom_nb as uint, ~[]);

    for self.edges.iter().advance |e|
    {
      if e.color() < 0
      { fail!("edge graph has not been colored correctly") }
      color_groups[e.color() as uint].push(e.content);
    }

    color_groups
  }


  pub fn export(&mut self) -> (~[Vec3<f64>], ~[i32], ~[i32])
  {
    let mut vertices: ~[Vec3<f64>] = ~[];
    let mut ids1: ~[i32] = ~[];
    let mut ids2: ~[i32] = ~[];

     self.unmark();
     for self.nodes.iter().advance |n|
     {
       vertices.push(Vec3::new(n.pos.x as f64,
                               n.pos.y as f64,
                               n.pos.z as f64));
     }

     for self.edges.iter().advance |e|
     {
       ids1.push(e.content.node_1.index() as i32);
       ids2.push(e.content.node_2.index() as i32);
     }

     (vertices, ids1, ids2)
  }

  pub fn build_edge_graph(&mut self)
  {
    self.unmark();
    for self.nodes.iter().advance |n|
    { Vertex::split(*n, &mut self.edges) }

    for self.nodes.iter().advance |n|
    {
      n.adj.clear();
      n.content.connect_edges();
    }
  }



  // Warning changes order of edges in edge array
  // DSATUR algorithm


  pub fn color_edge_graph(&mut self)
  {
    self.edge_chrom_nb = Node::color_graph(self.edges);
  }

  pub fn color_blob_graph(&mut self)
  {
    self.blob_chrom_nb = Node::color_graph(self.blobs);
    for self.blobs.iter().advance |b|
    {
      for b.content.sub_nodes.iter().advance |e|
      { e.set_color(b.color()) }
    }
  }

}
