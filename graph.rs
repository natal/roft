extern mod extra;
extern mod nalgebra;
extern mod kiss3d;

use std::io;
use std::uint;
use std::vec;
use extra::sort::Sort;
use extra::container::Deque;
use extra::ringbuf::RingBuf;
use nalgebra::vec::Vec3;
use nalgebra::traits::scalar_op::ScalarMul;
use kiss3d::window;
use kiss3d::object::VerticesNormalsTriangles;

type Vec3f = Vec3<f32>;

pub struct Node<T>
{
  priv id:      uint,
  priv color:   int,
  priv marked:  bool,
  priv dist:    uint,
  adj:          ~[@mut Node<T>],
  content:      T,
  pos:          Vec3f
}

impl<T> Ord for Node<T>
{
  fn lt(&self, other: &Node<T>) -> bool
  { self.adj.len() > other.adj.len() }

  fn le(&self, other: &Node<T>) -> bool
  { self.adj.len() >= other.adj.len() }

  fn ge(&self, other: &Node<T>) -> bool
  { self.adj.len() <= other.adj.len() }

  fn gt(&self, other: &Node<T>) -> bool
  { self.adj.len() < other.adj.len() }
}

impl<T> Eq for Node<T>
{
  fn eq(&self, other: &Node<T>) -> bool
  { self.adj.len() == other.adj.len() }

  fn ne(&self, other: &Node<T>) -> bool
  { self.adj.len() != other.adj.len() }
}

impl<T> Node<T>
{
  fn new(id: uint, content: T, pos: Vec3f) -> Node<T>
  {
    Node
    {
      id :      id,
      color :   -1,
      marked :  false,
      adj :     ~[],
      content : content,
      pos:      pos,
      dist:     0
    }
  }

  pub fn dsat(&self, nb_colors: uint) -> uint
  {
    let mut buckets = vec::from_elem(nb_colors + 1, false);

    for self.adj.iter().advance |n|
    {
      if (n.color >= 0)
      { buckets[n.color as uint] = true }
    }

    let mut count = 0u;

    for buckets.iter().advance |u|
    {
      if (*u)
      { count = count + 1 }
    }
    count
  }

  pub fn index(&self) -> uint
  { self.id as uint }

  pub fn degree(&self) -> uint
  { self.adj.len() }

  pub fn reset_dist(&mut self)
  { self.dist = 0 }

  pub fn unmark(&mut self)
  { self.marked = false }

  pub fn mark(&mut self)
  { self.marked = true }

  pub fn is_marked(&self) -> bool
  { self.marked }

  pub fn to_str(&self) -> ~str
  { self.id.to_str() }

  pub fn equals(&self, n: &Node<T>) -> bool
  { n.id == self.id }

  pub fn color(&self) -> int
  { self.color }

  pub fn set_color(&mut self, col: int)
  { self.color = col; }

  pub fn color_with_min(&mut self) -> int
  {
    let mut max_col = self.adj.head().color;
    for self.adj.iter().advance |n|
    {
      if (max_col < n.color)
      { max_col = n.color }
    }
    self.color = max_col + 1;
    self.color
  }

  pub fn is_adj_to(&self, node: &Node<T>) -> bool
  {
    for self.adj.iter().advance |n|
    {
      if n.equals(node)
      { return true }
    }
    return false;
  }

  pub fn nb_common_adj(&self, node: &Node<T>) -> uint
  {
    let mut nb_common = 0;
    for self.adj.iter().advance |n1|
    {
      for node.adj.iter().advance |n2|
      {
        if n1.equals(*n2)
        { nb_common = nb_common + 1}
      }
    }
    nb_common
  }

  priv fn share_k_adjs(&self, node: &Node<T>, k: uint) -> bool
  {
    let mut count = 0;
    for self.adj.iter().advance |n1|
    {
      for node.adj.iter().advance |n2|
      {
        if n1.equals(*n2)
        { count = count + 1 }
        if count >= k
        {  return true }
      }
    }
    return false
  }

  // Must be used with graph unmarked
  pub fn distant_nodes(@mut self, pred: &fn (&Node<T>) -> bool, d: uint) -> ~[@mut Node<T>]
  {
    let mut node_queue: ~RingBuf<@mut Node<T>> = ~RingBuf::new();
    let mut res_nodes:  ~[@mut Node<T>] = ~[];

    self.marked = true;
    self.dist   = 0;

    node_queue.push_back(self);

    while !node_queue.is_empty()
    {
      let n = node_queue.pop_front().unwrap();
      for uint::iterate(0u, n.adj.len()) |i|
      {
        let n2 : @mut Node<T> = n.adj[i];
        if !n2.marked
        {
          n2.dist = n.dist + 1;
          if n2.dist < d
          { node_queue.push_back(n2) }
          else if n2.dist == d && pred(n2)
          { res_nodes.push(n2) }
        }
        n2.marked = true;
      }
    }
    res_nodes
  }

  pub fn connect(n1: @mut Node<T>, n2: @mut Node<T>)
  {
    if !n2.is_adj_to(n1)
    {
      n1.adj.push(n2);
      n2.adj.push(n1);
    }
  }
}


pub struct Edge
{
  node_1: @mut Node<Vertex>,
  node_2: @mut Node<Vertex>
}

impl Edge
{
  pub fn new(n1: @mut Node<Vertex>, n2: @mut Node<Vertex>) -> Edge
  {
    Edge
    {
      node_1: n1,
      node_2: n2
    }
  }
}

pub struct EdgeSet
{
  edges: ~[@mut Node<Edge>]
}

impl EdgeSet
{
  pub fn new(ne: @mut Node<Edge>) -> EdgeSet
  {
    EdgeSet
    {
      edges: ~[ne]
    }
  }

  pub fn is_singleton(&self) -> bool
  {
    self.edges.len() == 1
  }

}

pub struct Vertex
{
  edges: ~[@mut Node<Edge>],
  pos:   Vec3f
}

impl Vertex
{
  pub fn new(pos: Vec3f) -> Vertex
  {
    Vertex
    {
      edges: ~[],
      pos:   pos
    }
  }

  pub fn connect_edges(&mut self)
  {
    for self.edges.iter().advance |e1|
    {
      for self.edges.iter().advance |e2|
      {
        if (e1.id != e2.id) && !e1.is_adj_to(*e2)
        { Node::connect(*e1, *e2) }
      }
    }
  }

  //To be used with graph unmarked
  pub fn split(node: @mut Node<Vertex>, all_edges: &mut ~[@mut Node<Edge>])
  {
    for uint::iterate(0u, node.adj.len()) |i|
    {
      let a = node.adj[i];
      if !a.is_marked()
      {
        let edge = Edge::new(a, node);
        let pos = (a.content.pos + node.content.pos).scalar_mul(&0.5);
        let edge_node = @mut Node::new(all_edges.len(), edge, pos);

        node.content.edges.push(edge_node);
        a.content.edges.push(edge_node);
        all_edges.push(edge_node);
      }
    }
    node.mark();
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
  edge_sets: ~[@mut Node<EdgeSet>]
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
      edge_sets: ~[]
    }
  }

  pub fn augment(&mut self)
  {
    let mut to_connect : ~[~[@mut Node<Vertex>]] = ~[];
    for uint::iterate(0u, self.nodes.len()) |i|
    {
      self.unmark();
      let n1 = self.nodes[i];
      to_connect.push(n1.distant_nodes(|n2| n1.share_k_adjs(n2, 2), 2));
    }

    for self.nodes.iter().enumerate().advance |(i, n)|
    {
      for to_connect[i].iter().advance |n2|
      { Node::connect(*n, *n2) }
    }
  }

  pub fn unmark(&mut self)
  {
     for self.nodes.iter().advance |n1|
     { n1.unmark() }
     for self.edges.iter().advance |e1|
     { e1.unmark() }
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
                      colors[e.color() as uint] +
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


  // Warning : erases color
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

  // Requires to have built the edge graph
//  pub fn build_edge_set_graph(&mut self, dist: uint)
//  {
//    self.unmark();
//
//    for self.edges.iter().enumerate().advance |(i, e)|
//    { edge_sets.push(Node::new(i, EdgeSet::new(e), e.pos)); }
//
//    for self.edges.iter().advance |e|
//    {
//      let dist_nodes = e.distant_nodes(|e2| es.is_singleton(), dist);
//      for dist_nodes.iter().advance |dn|
//      {
//      }
//    }
//
//  }

  // Warning changes order of edges in edge array
  // DSATUR algorithm

  pub fn color_edge_graph(&mut self)
  {
    assert!(self.edges.len() > 0)

    self.edges.qsort();
    self.edges.head().set_color(0);

    let mut nb_chrom : int = 0;
    let mut uncolored = self.edges.len();


    while uncolored > 1
    {
      assert!(nb_chrom >= 0);

      let mut max_node = self.edges.iter().find_(|n| n.color() < 0).unwrap();
      let mut max_dsat = max_node.dsat(nb_chrom as uint);
      for self.edges.iter().advance |n|
      {
        if (n.color() < 0)
        {
          let cur_dsat = n.dsat(nb_chrom as uint);
          if (max_dsat < cur_dsat) || (max_dsat == cur_dsat) && (max_node.degree() < n.degree())
          {
            max_dsat = cur_dsat;
            max_node = n;
          }
        }
      }
      uncolored = uncolored - 1;
      nb_chrom = max_node.color_with_min().max(&nb_chrom);
    }
    println("nb_chrom : " + nb_chrom.to_str());
    println("mean : " + (self.edges.len() as float / (nb_chrom as float)).to_str());
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
