extern mod nalgebra;

use nalgebra::vec::Vec3;
use nalgebra::traits::scalar_op::ScalarMul;
use edge::Edge;
use node::Node;

type Vec3f = Vec3<f32>;

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
    for e1 in self.edges.iter()
    {
      for e2 in self.edges.iter()
      {
        if (!e1.equals_(*e2)) && !e1.is_adj_to(*e2)
        { Node::connect(*e1, *e2) }
      }
    }
  }

  //To be used with graph unmarked
  pub fn split(node: @mut Node<Vertex>, all_edges: &mut ~[@mut Node<Edge>])
  {
    for i in range(0u, node.adj.len())
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
