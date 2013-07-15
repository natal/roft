extern mod extra;
extern mod nalgebra;
extern mod kiss3d;

use std::uint;
use std::vec;
use nalgebra::vec::Vec3;
use extra::container::Deque;
use extra::ringbuf::RingBuf;

type Vec3f = Vec3<f32>;

pub struct Node<T>
{
  id:      uint,
  priv color:   int,
  priv marked:  bool,
  priv intern_mark: bool,
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
  pub fn new(id: uint, content: T, pos: Vec3f) -> Node<T>
  {
    Node
    {
      id :      id,
      color :   -1,
      marked :  false,
      adj :     ~[],
      content : content,
      pos:      pos,
      dist:     0,
      intern_mark: false
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

  pub fn dist(&self) -> uint
  {
    self.dist
  }

  pub fn index(&self) -> uint
  { self.id as uint }

  pub fn degree(&self) -> uint
  { self.adj.len() }

  pub fn reset_dist(&mut self)
  { self.dist = 0 }

  pub fn unmark(&mut self)
  { self.marked = false }

  pub fn intern_unmark(&mut self)
  { self.intern_mark = false }

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

  pub fn share_k_adjs(&self, node: &Node<T>, k: uint) -> bool
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

    self.intern_mark = true;
    self.dist   = 0;

    node_queue.push_back(self);

    while !node_queue.is_empty()
    {
      let n = node_queue.pop_front().unwrap();
      for uint::iterate(0u, n.adj.len()) |i|
      {
        let n2 : @mut Node<T> = n.adj[i];
        if !n2.intern_mark
        {
          n2.dist = n.dist + 1;
          if n2.dist < d
          { node_queue.push_back(n2) }
          if pred(n2)
          { res_nodes.push(n2) }
        }
        n2.intern_mark = true;
      }
    }
    res_nodes
  }

  pub fn connect(n1: @mut Node<T>, n2: @mut Node<T>)
  {
    if !n1.equals(n2) && !n2.is_adj_to(n1)
    {
      n1.adj.push(n2);
      n2.adj.push(n1);
    }
  }

  pub fn id(&self) -> uint
  {
    self.id
  }

  pub fn eat(n1: @mut Node<T>, n2: @mut Node<T>)
  {
    let mut to_connect: ~[@mut Node<T>] = ~[];
    for n1.adj.iter().advance |n3|
    { to_connect.push(*n3) }

    for to_connect.iter().advance |n3|
    { Node::connect(n1, *n3) }

    n2.adj.clear();
  }
}

