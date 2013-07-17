extern mod extra;
extern mod nalgebra;
extern mod kiss3d;

use std::io;
use std::uint;
use nalgebra::vec::Vec3;
use kiss3d::window;
use kiss3d::object::VerticesNormalsTriangles;
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

pub struct Graph
{
  nodes: ~[@mut Node<Vertex>],
  edges: ~[@mut Node<Edge>],
  blobs: ~[@mut Node<Blob<Edge>>]
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
      blobs: ~[]
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


  pub fn write_simple(&mut self)
  {
     let path = Path("./simple.dot");
     let file = io::file_writer(&path, [io::Create]).get();
     file.write_str("graph simple {\n");
     self.unmark();
     for self.nodes.iter().advance |n1|
     {
       file.write_str(n1.to_str() + " [pos=\"" + n1.pos.at[0].to_str() + "," +
                      n1.pos.at[1].to_str() + "!\"]\n");
     }
     for self.edges.iter().advance |e|
     {
       file.write_str(e.content.node_1.to_str() + " -- " +
                      e.content.node_2.to_str() + "\n");
     }

     for self.nodes.iter().advance |e1|
     {
       for e1.adj.iter().advance |e2|
       {
         if (!e2.is_marked())
         {
           file.write_str(e1.to_str() + " -- " +
                          e2.to_str() + "\n");
         }
       }
       e1.mark();
     }
     file.write_str("}");
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


  pub fn export(&mut self) -> (~[Vec3<f64>], ~[(uint, uint)])
  {
    let mut vertices: ~[Vec3<f64>] = ~[];
    let mut edges: ~[(uint, uint)] = ~[];

     self.unmark();
     for self.nodes.iter().advance |n|
     {
       vertices.push(Vec3::new([n.pos.at[0] as f64,
                                n.pos.at[1] as f64,
                                n.pos.at[2] as f64]));
     }

     for self.edges.iter().advance |e|
     { edges.push((e.content.node_1.index(), e.content.node_2.index())) }

     (vertices, edges)

  }

  pub fn write_line_graph(&mut self)
  {
     let path = Path("./line.dot");
     let file = io::file_writer(&path, [io::Create]).get();
     let colors = ~["azure", "skyblue", "pink", "crimson", "peru",
                    "orange", "gold", "lawngreen", "cyan", "blueviolet",
                    "lavender", "mediumblue", "limegreen", "chocolate", "plum",
                    "yellowgreen", "royalblue", "hotpink", "darkslategray",
                    "darkorange", "beige", "aliceblue", "tomato", "salmon"];
     file.write_str("graph line {\n");
     self.unmark();
     for self.edges.iter().advance |e|
     {
       file.write_str(e.to_str() + " [pos=\"" + e.pos.at[0].to_str() + "," +
                      e.pos.at[1].to_str() + "!\", color="+
                      colors[(e.color() + 1) as uint] +
                      ", style=filled]\n");
     }

     for self.edges.iter().advance |e1|
     {
       for e1.adj.iter().advance |e2|
       {
         if (!e2.is_marked())
         {
           file.write_str(e1.to_str() + " -- " +
                          e2.to_str() + "\n");
         }
       }
       e1.mark();
     }
     file.write_str("}");
  }

  pub fn write_blob_graph(&mut self)
  {
     let path = Path("./blob.dot");
     let file = io::file_writer(&path, [io::Create]).get();
     let colors = ~["azure", "skyblue", "pink", "crimson", "peru",
                    "orange", "gold", "lawngreen", "cyan", "blueviolet",
                    "lavender", "mediumblue", "limegreen", "chocolate", "plum",
                    "yellowgreen", "royalblue", "hotpink", "darkslategray",
                    "darkorange", "beige", "aliceblue", "tomato", "salmon"];
     file.write_str("graph blob {\n");
     self.unmark();

     for self.blobs.iter().advance |e|
     {
       file.write_str(e.to_str() + " [pos=\"" + e.pos.at[0].to_str() + "," +
                      e.pos.at[1].to_str() + "!\", color="+
                      colors[(e.color() + 1) as uint] +
                      ", style=filled]\n");
     }

     for self.blobs.iter().advance |b1|
     {
       for b1.adj.iter().advance |b2|
       {
         if (!b2.is_marked())
         {
           file.write_str(b1.to_str() + " -- " +
                          b2.to_str() + "\n");
         }
       }
       b1.mark();
     }
     file.write_str("}");
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
    Node::color_graph(self.edges);
  }

  pub fn color_blob_graph(&mut self)
  {
    Node::color_graph(self.blobs);
    for self.blobs.iter().advance |b|
    {
      for b.content.sub_nodes.iter().advance |e|
      { e.set_color(b.color()) }
    }
  }

}



fn main()
{
  do window::Window::spawn(~"Mesh") |w|
  {
    let q = w.add_quad(10.0, 10.0, 5, 5);
    match *q.geometry()
    {
      VerticesNormalsTriangles(ref v, _, ref t) =>
      {
        let mesh = Mesh::new(v.clone(), t.clone());
        let mut graph = Graph::new(mesh);
        println("Augmenting graph...");
        graph.augment();
        println("Creating line graph...");
        graph.build_edge_graph();
        println("Coloring edge graph...");
//        graph.color_edge_graph();
        println("Writting to file...");
        graph.write_simple();
//        graph.write_line_graph();
        let (v, e) = graph.export();
        println(v.to_str());
        println(e.to_str());
        println("Done");
      }
      _ => { }
    }
  }
}
