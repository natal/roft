use node::Node;
use vertex::Vertex;


#[deriving(Clone)]
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
