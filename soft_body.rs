use std::num::{Zero, abs};
use nalgebra::vec::Vec0;
use nalgebra::traits::division_ring::DivisionRing;
use nalgebra::traits::norm::Norm;
use nalgebra::traits::dot::Dot;
use nalgebra::traits::vector_space::VectorSpace;
use nphysics::constraint::velocity_constraint::VelocityConstraint;
use nphysics::constraint::projected_gauss_seidel_solver::projected_gauss_seidel_solve;

pub struct PointMass<N, V>
{
  invmass:  N,
  velocity: V,
  position: V,
}

pub struct ConstraintsGeometry<N>
{
  stiffness:   N,
  rest_length: N,
  impulse:     N,
  rb1:         uint,
  rb2:         uint
}

pub struct SoftBody<N, V>
{
  points:      ~[PointMass<N, V>],
  constraints: ~[ConstraintsGeometry<N>]
}

impl<N: DivisionRing + Eq + Ord + Clone,
     V: VectorSpace<N> + Norm<N> + Dot<N> + Clone>
    SoftBody<N, V>
{
  pub fn from_mesh(vbuf:      ~[V],
                   edges:     ~[(uint, uint)],
                   invmasses: ~[N],
                   stiffness: ~[N]) -> SoftBody<N, V>
  {
    assert!(vbuf.len() == invmasses.len(),
            "Vertex buffer and mass informations must have the same size.");

    assert!(stiffness.len() == edges.len(),
            "Edge buffer and stiffness informations must have the same size.");

    // create points mass
    let mut points = ~[];

    for vbuf.iter().zip(invmasses.iter()).advance |(v, m)|
    {
      points.push(PointMass {
        invmass:  m.clone(),
        velocity: Zero::zero(),
        position: v.clone(),
      });
    }

    // create constraints
    let mut constraints = ~[];

    for edges.iter().zip(stiffness.iter()).advance |(&(v1, v2), s)|
    {
      constraints.push(ConstraintsGeometry {
        stiffness:   s.clone(),
        rest_length: (vbuf[v1] - vbuf[v2]).norm(),
        impulse:     Zero::zero(),
        rb1:         v1,
        rb2:         v2
      });
    }

    SoftBody {
      points:      points,
      constraints: constraints
    }
  }

  pub fn integrate(&mut self, dt: &N, fext: &V)
  {
    for self.points.mut_iter().advance |p|
    {
      if !p.invmass.is_zero()
      {
        p.velocity = p.velocity + fext.scalar_mul(dt);
        p.position = p.position + p.velocity.scalar_mul(dt);
      }
    }
  }

  pub fn collect_constraints(&self, out: &mut ~[VelocityConstraint<V, Vec0<N>, N>])
  {
    for self.constraints.iter().advance |c|
    {
      let mut normal = self.points[c.rb1].position - self.points[c.rb2].position;
      let     length = normal.normalize();

      let m1 = self.points[c.rb1].invmass.clone();
      let m2 = self.points[c.rb2].invmass.clone();

      let limit = abs(c.stiffness * (length - c.rest_length));

      if true // length != c.rest_length
      {
        out.push(
          VelocityConstraint {
            weighted_normal1:   normal.scalar_mul(&m1),
            weighted_normal2:   normal.scalar_mul(&m2),

            rot_axis1:          Zero::zero(),
            weighted_rot_axis1: Zero::zero(),

            rot_axis2:          Zero::zero(),
            weighted_rot_axis2: Zero::zero(),

            projected_mass:     m1 + m2,

            impulse:            c.impulse.clone(), 
            unit_impulse:       Zero::zero(),
            lobound:            -limit,
            hibound:            limit,
            objective:          -(self.points[c.rb2].velocity -
                                  self.points[c.rb1].velocity).dot(&normal),
            id1:                if m1.is_zero() { -1 } else { c.rb1 as int },
            id2:                if m2.is_zero() { -1 } else { c.rb2 as int },

            normal:             normal
          }
        )
      }
    }
  }
}

impl<V: VectorSpace<N> + Dot<N> + Norm<N> + Copy + Clone + ToStr,
     N:  DivisionRing + Orderable + Ord + ToStr + Eq + Clone + Copy>
     SoftBody<N, V>
{
  pub fn solve(&mut self)
  {
    let mut constraints = ~[];

    self.collect_constraints(&mut constraints);
  
    let res = projected_gauss_seidel_solve(constraints,
                                           self.points.len(),
                                           50,
                                           false);

    for self.points.mut_iter().enumerate().advance |(i, p)|
    { p.velocity = p.velocity + res[i].lv }

    for constraints.iter().enumerate().advance |(i, c)|
    { self.constraints[i].impulse = c.impulse.clone() }
  }
}
