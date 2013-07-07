extern mod nalgebra;

use std::iterator::Counter;
use std::io;
use std::uint;
use nalgebra::vec::Vec3;
use nalgebra::traits::scalar_op::ScalarMul;

type Vec3f = Vec3<f64>;


pub struct Edge
{
  priv id:      uint,
  priv color:   int,
  adj_edges:    ~[@mut Edge],
  node_1:       @mut Node,
  node_2:       @mut Node,
  pos:          Vec3f
}

impl Edge
{
  pub fn new(identifier: uint, col: int, n1: @mut Node, n2: @mut Node) -> Edge
  {
    Edge
    {
      id:    identifier,
      adj_edges:   ~[],
      node_1: n1,
      node_2: n2,
      color: col,
      pos: (n1.pos + n2.pos).scalar_mul(&0.5)
    }
  }
  pub fn to_str(&self) -> ~str
  {
    "e" + self.id.to_str()
  }

  pub fn equals(&self, e: &Edge) -> bool
  {
    e.id == self.id
  }

  pub fn color(&self) -> int
  {
    self.color
  }

  pub fn set_color(&mut self, col: int)
  {
    self.color = col;
  }

  pub fn is_adj_to(&self, edge: &Edge) -> bool
  {
    for self.adj_edges.iter().advance |n1|
    {
      if n1.equals(edge)
      { return true }
    }
    return false;
  }
}


pub struct Node
{
  priv id:      uint,
  pos:          Vec3f,
  adj_edges:    ~[@mut Edge],
  adj_nodes:    ~[@mut Node],
  priv marked: bool
}

impl Node
{
  pub fn new(identifier: uint, pos: Vec3f) -> Node
  {
    Node
    {
      id:    identifier,
      adj_edges:   ~[],
      adj_nodes:   ~[],
      pos: pos,
      marked: false
    }
  }

  pub fn unmark(&mut self)
  {
    self.marked = false
  }

  pub fn mark(&mut self)
  {
    self.marked = true
  }

  pub fn is_marked(&self) -> bool
  {
    self.marked
  }

  pub fn to_str(&self) -> ~str
  {
    self.id.to_str()
  }

  pub fn equals(&self, e: &Node) -> bool
  {
    e.id == self.id
  }

  pub fn clear_adj_nodes(&mut self)
  {
    self.adj_nodes.clear();
  }

  pub fn is_adj_to(&self, node: &Node) -> bool
  {
    self.adj_nodes.iter().any_(|n1| n1.equals(node))
  }

  pub fn nb_common_adj(&self, node: &Node) -> uint
  {
    let mut nb_common = 0;
    for self.adj_nodes.iter().advance |n1|
    {
      for node.adj_nodes.iter().advance |n2|
      {
        if n1.equals(*n2)
        { nb_common = nb_common + 1}
      }
    }
    nb_common
  }

  priv fn share_2_adjs(&self, node: &Node) -> bool
  {
    let mut count = 0;
    for self.adj_nodes.iter().advance |n1|
    {
      for node.adj_nodes.iter().advance |n2|
      {
        if n1.equals(*n2)
        { count = count + 1 }
        if count >= 2
        { return true }
      }
    }
    return false
  }

  pub fn get_2_distant_adjs(&self) -> ~[@mut Node]
  {
    let mut dist_nodes : ~[@mut Node] = ~[];
    for self.adj_nodes.iter().advance |n1|
    {
      for n1.adj_nodes.iter().advance |n2|
      {
        if (!self.is_adj_to(*n2)) && (self.share_2_adjs(*n2)) && (!self.equals(*n2))
        { dist_nodes.push(*n2) }
      }
    }
    dist_nodes
  }

  pub fn connect_nodes(n1: @mut Node, n2: @mut Node)
  {
    if !n2.is_adj_to(n1)
    {
      n1.adj_nodes.push(n2);
      n2.adj_nodes.push(n1);
    }
  }


  pub fn connect_edges(&mut self)
  {
    for self.adj_edges.iter().advance |e1|
    {
      for self.adj_edges.iter().advance |e2|
      {
        if (e1.id != e2.id) && !e1.is_adj_to(*e2)
        {
          e1.adj_edges.push(*e2);
          e2.adj_edges.push(*e1);
        }
      }
    }
  }

  pub fn split(node: @mut Node, all_edges: &mut ~[@mut Edge])
  {
    for uint::iterate(0u, node.adj_nodes.len()) |i|
    {
      let a = node.adj_nodes[i];
      if !a.marked
      {
        let edge = @mut Edge::new(all_edges.len() + 1, -1, a, node);

        node.adj_edges.push(edge);
        a.adj_edges.push(edge);

        all_edges.push(edge);
      }
    }
    node.marked = true;
  }
}

pub struct Mesh
{
  vbuff: ~[Vec3f],
  ibuff: ~[uint]
}

impl Mesh
{
  pub fn new(vb: ~[Vec3f], ib: ~[uint]) -> Mesh
  {
    assert!(ib.len() % 3 == 0);
    Mesh
    {
      vbuff: vb,
      ibuff: ib
    }
  }
}

pub struct Graph
{
  nodes: ~[@mut Node],
  edges: ~[@mut Edge]
}

impl Graph
{
  pub fn new(mesh: Mesh) -> Graph
  {
    let mut nodes: ~[@mut Node] = ~[];
    for mesh.vbuff.iter().enumerate().advance |(vid, v)|
    {
      nodes.push(@mut Node::new (vid, *v));
    }

    for Counter::new(0u, 3).advance |i|
    {
      if i >= mesh.ibuff.len()
      { break }

      let id1 = mesh.ibuff[i];
      let id2 = mesh.ibuff[i + 1];
      let id3 = mesh.ibuff[i + 2];

      Node::connect_nodes(nodes[id1], nodes[id2]);
      Node::connect_nodes(nodes[id1], nodes[id3]);
      Node::connect_nodes(nodes[id2], nodes[id1]);
      Node::connect_nodes(nodes[id2], nodes[id3]);
      Node::connect_nodes(nodes[id3], nodes[id2]);
      Node::connect_nodes(nodes[id3], nodes[id1]);
    }

    for nodes.iter().advance |n|
    {
      println(n.to_str() + " has " + n.adj_nodes.len().to_str() + " neighbors");
    }

    Graph
    {
      nodes: nodes,
      edges: ~[]
    }
  }

  pub fn augment(&mut self)
  {
    let mut to_connect : ~[~[@mut Node]] = ~[];
    for self.nodes.iter().advance |n|
    {
      to_connect.push(n.get_2_distant_adjs());
    }
    for self.nodes.iter().enumerate().advance |(i, n)|
    {
      for to_connect[i].iter().advance |n2|
      { Node::connect_nodes(*n, *n2) }
    }
  }

  pub fn unmark(&mut self)
  {
     for self.nodes.iter().advance |n1|
     {
       n1.unmark();
     }
  }

  pub fn write_to_file(&mut self)
  {
     let path = Path("./out.dot");
     let file = io::file_writer(&path, [io::Create]).get();
     file.write_str("graph graphname {\n");

     self.unmark();
     for self.nodes.iter().advance |n1|
     {
       file.write_str(n1.to_str() + " [pos=\"" + n1.pos.at[0].to_str() + "," +
                      n1.pos.at[1].to_str() + "!\"]\n");
     }

     for self.edges.iter().advance |e|
     {
       file.write_str(e.to_str() + " [pos=\"" + e.pos.at[0].to_str() + "," +
                      e.pos.at[1].to_str() + "!\"]\n");
     }

     for self.nodes.iter().advance |n1|
     {
       for n1.adj_edges.iter().advance |e|
       {
         file.write_str(n1.to_str() + " -- " +
                        e.to_str() + "\n");
       }
     }

     for self.nodes.iter().advance |n1|
     {
       for n1.adj_nodes.iter().advance |n2|
       {
         if (!n2.is_marked())
         {
           file.write_str(n1.to_str() + " -- " +
                          n2.to_str() + "\n");
         }
       }
       n1.mark();
     }
     file.write_str("}");
  }

  // Warning : erases color
  pub fn build_edge_graph(&mut self)
  {
    self.unmark();
    for self.nodes.iter().advance |n|
    {
      Node::split(*n, &mut self.edges);
    }

    for self.nodes.iter().advance |n|
    {
      n.clear_adj_nodes();
      n.connect_edges();
    }
  }

  pub fn color(&self)
  {
  }
}

fn main()
{
  let vb = ~[Vec3::new([0.0f64,0.0,0.0]), Vec3::new([0.0f64,2.0,0.0]), Vec3::new([0.0f64,4.0,0.0]),Vec3::new([0.0f64,6.0,0.0]),
            Vec3::new([2.0f64,0.0,0.0]), Vec3::new([2.0f64,2.0,0.0]), Vec3::new([2.0f64,4.0,0.0]), Vec3::new([2.0f64,6.0,0.0]),
            Vec3::new([4.0f64,0.0,0.0]), Vec3::new([4.0f64,2.0,0.0]), Vec3::new([4.0f64,4.0,0.0]), Vec3::new([4.0f64,6.0,0.0])];

  let ib = ~[0u, 1, 5,
             1, 2, 5,
             2, 3, 6,
             0, 4, 5,
             2, 5, 6,
             3, 6, 7,
             4, 5, 8,
             8, 9, 5,
             9, 5, 6,
             9, 10, 6,
             6, 10, 11,
             6, 7, 11];

  let mesh = Mesh::new(vb, ib);
  let mut graph = Graph::new(mesh);
  graph.augment();
  graph.build_edge_graph();
  graph.write_to_file();
}
