extern mod extra;
extern mod nalgebra;

use std::iterator::Counter;
use std::io;
use std::uint;
use extra::sort::Sort;
use nalgebra::vec::Vec3;
use nalgebra::traits::scalar_op::ScalarMul;

type Vec3f = Vec3<f64>;


pub struct Edge
{
  priv id:      uint,
  priv color:   int,
  priv marked:  bool,
  adj_edges:    ~[@mut Edge],
  node_1:       @mut Node,
  node_2:       @mut Node,
  pos:          Vec3f
}

impl Ord for Edge
{
  fn lt(&self, other: &Edge) -> bool
  {
    self.adj_edges.len() > other.adj_edges.len()
  }

  fn le(&self, other: &Edge) -> bool
  {
    self.adj_edges.len() >= other.adj_edges.len()
  }

  fn ge(&self, other: &Edge) -> bool
  {
    self.adj_edges.len() <= other.adj_edges.len()
  }

  fn gt(&self, other: &Edge) -> bool
  {
    self.adj_edges.len() < other.adj_edges.len()
  }
}

impl Eq for Edge
{
  fn eq(&self, other: &Edge) -> bool
  {
    self.adj_edges.len() == other.adj_edges.len()
  }

  fn ne(&self, other: &Edge) -> bool
  {
    self.adj_edges.len() != other.adj_edges.len()
  }
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
      marked: false,
      pos: (n1.pos + n2.pos).scalar_mul(&0.5)
    }
  }

  pub fn dsat(&self, nb_colors: uint) -> uint
  {
    let mut buckets : ~[bool] = ~[];
    for uint::iterate(0u, nb_colors + 1) |i|
    { buckets.push(false) }

    for self.adj_edges.iter().advance |e|
    {
      if (e.color >= 0)
      { buckets[e.color as uint] = true }
    }

    let mut count = 0u;

    for buckets.iter().advance |u|
    {
      if (*u)
      { count = count + 1 }
    }
    count
  }

  pub fn degree(&self) -> uint
  {
    self.adj_edges.len()
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

  pub fn min_color(&mut self) -> int
  {
    let mut max_col = self.adj_edges.head().color;
    for self.adj_edges.iter().advance |e|
    {
      if (max_col < e.color)
      { max_col = e.color }
    }

    self.color = max_col + 1;
    self.color
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
     for self.edges.iter().advance |e1|
     {
       e1.unmark();
     }
  }

  pub fn write_to_file(&mut self)
  {
     let path = Path("./out.dot");
     let file = io::file_writer(&path, [io::Create]).get();
     let colors = ~["azure", "skyblue", "pink", "crimson", "peru",
                    "orange", "gold", "lawngreen", "cyan", "blueviolet",
                    "lavender", "mediumblue", "limegreen"];
     file.write_str("graph graphname {\n");

     self.unmark();
 //    for self.nodes.iter().advance |n1|
 //    {
 //      file.write_str(n1.to_str() + " [pos=\"" + n1.pos.at[0].to_str() + "," +
 //                     n1.pos.at[1].to_str() + "!\"]\n");
 //    }

     for self.edges.iter().advance |e|
     {
       file.write_str(e.to_str() + " [pos=\"" + e.pos.at[0].to_str() + "," +
                      e.pos.at[1].to_str() + "!\", color="+
                      colors[e.color() as uint] +
                      ", style=filled]\n");
     }

 //    for self.nodes.iter().advance |n1|
 //    {
 //      for n1.adj_edges.iter().advance |e|
 //      {
 //        file.write_str(n1.to_str() + " -- " +
 //                       e.to_str() + "\n");
 //      }
 //    }

     for self.edges.iter().advance |e1|
     {
       for e1.adj_edges.iter().advance |e2|
       {
         if (!e2.is_marked())
         {
           file.write_str(e1.to_str() + " -- " +
                          e2.to_str() + "\n");
         }
       }
       e1.mark();
     }

 //    for self.nodes.iter().advance |n1|
 //    {
 //      for n1.adj_nodes.iter().advance |n2|
 //      {
 //        if (!n2.is_marked())
 //        {
 //          file.write_str(n1.to_str() + " -- " +
 //                         n2.to_str() + "\n");
 //        }
 //      }
 //      n1.mark();
 //    }
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
      for self.edges.iter().enumerate().advance |(i, n)|
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
      nb_chrom = max_node.min_color().max(&nb_chrom);
    }
    println("nb_chrom : " + nb_chrom.to_str());
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

//  let vb = ~[Vec3::new([0.0f64,0.0,0.0]), Vec3::new([0.0f64,2.0,0.0]), Vec3::new([0.0f64,4.0,0.0]),
//            Vec3::new([2.0f64,0.0,0.0]), Vec3::new([2.0f64,2.0,0.0]), Vec3::new([2.0f64,4.0,0.0]) ,
//            Vec3::new([4.0f64,0.0,0.0]), Vec3::new([4.0f64,2.0,0.0]), Vec3::new([4.0f64,4.0,0.0]) ];
//
//  let ib = ~[0u, 1, 4,
//             1, 2, 4,
//             4, 2, 5,
//             0, 3, 4,
//             3, 4, 6,
//             6, 7, 4,
//             4, 5, 7,
//             7, 5, 8];

  let mesh = Mesh::new(vb, ib);
  let mut graph = Graph::new(mesh);
  graph.augment();
  graph.build_edge_graph();
  graph.color_edge_graph();
  graph.write_to_file();
}
