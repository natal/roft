use std::uint;
use std::num::Zero;
use nalgebra::vec::Vec0;
use nalgebra::traits::division_ring::DivisionRing;
use nalgebra::traits::norm::Norm;
use nalgebra::traits::dot::Dot;
use nalgebra::traits::vector_space::VectorSpace;
use nphysics::constraint::velocity_constraint::VelocityConstraint;
use nphysics::constraint::projected_gauss_seidel_solver::projected_gauss_seidel_solve;

pub struct PointMass<N, V>
{
  invmass:    N,
  velocity:   V,
  position:   V,
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
  ext_forces:  V,
  points:      ~[PointMass<N, V>],
  constraints: ~[ConstraintsGeometry<N>]
}

impl<N: DivisionRing + NumCast + Signed + Bounded + Eq + Ord + Clone,
     V: VectorSpace<N> + Norm<N> + Dot<N> + Clone>
    SoftBody<N, V>
{
  pub fn from_mesh(vbuf:      ~[V],
                   ids1:      ~[i32],
                   ids2:      ~[i32],
                   invmasses: ~[N],
                   stiffness: ~[N]) -> SoftBody<N, V>
  {
    assert!(vbuf.len() == invmasses.len(),
            "Vertex buffer and mass informations must have the same size.");

    // create points mass
    let mut points = ~[];

    for vbuf.iter().zip(invmasses.iter()).advance |(v, m)|
    {
      points.push(PointMass {
        invmass:    m.clone(),
        velocity:   Zero::zero(),
        position:   v.clone(),
      });
    }

    // create constraints
    let mut constraints = ~[];

    for uint::iterate(0u, ids2.len()) |i|
    {
      let v1 = ids1[i];
      let v2 = ids2[i];
      let s  = stiffness[i].clone();

      constraints.push(ConstraintsGeometry {
        stiffness:   s.clone(),
        rest_length: (vbuf[v1] - vbuf[v2]).norm(),
        impulse:     Zero::zero(),
        rb1:         v1 as uint,
        rb2:         v2 as uint
      });
    }

    SoftBody {
      points:      points,
      constraints: constraints,
      ext_forces:  Zero::zero()
    }
  }

  pub fn integrate(&mut self, dt: &N, fext: &V)
  {
    self.ext_forces = fext.clone();

    for self.points.mut_iter().advance |p|
    {
      if !p.invmass.is_zero()
      {
        p.velocity = p.velocity + fext.scalar_mul(dt);;
        p.position = p.position + p.velocity.scalar_mul(dt);
      }
    }
  }

  pub fn collect_constraints(&self,
                             dt:          N,
                             out:         &mut ~[VelocityConstraint<V, Vec0<N>, N>],
                             first_order: bool)
  {
    for self.constraints.iter().advance |c|
    {
      let mut normal = self.points[c.rb1].position - self.points[c.rb2].position;
      let     length = normal.normalize();

      let m1 = self.points[c.rb1].invmass.clone();
      let m2 = self.points[c.rb2].invmass.clone();

      let mut dvel = Zero::zero::<N>();

      if !first_order
      {
        dvel = dvel + dt * ((length - c.rest_length) * c.stiffness);

        if !m2.is_zero()
        { dvel = dvel - (self.points[c.rb2].velocity + self.ext_forces.scalar_mul(&dt)).dot(&normal) }

        if !m1.is_zero()
        { dvel = dvel + (self.points[c.rb1].velocity + self.ext_forces.scalar_mul(&dt)).dot(&normal) }
      }
      else
      {
        dvel = (length - c.rest_length) * NumCast::from::<N, float>(0.4) / dt;
      }

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

            impulse:            if false { Zero::zero() } else { c.impulse.clone() },
            unit_impulse:       Zero::zero(),
            lobound:            -Bounded::max_value::<N>(), // limit,
            hibound:            Bounded::max_value(), // limit,
            objective:          dvel,
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
     N:  DivisionRing + Orderable + NumCast + Signed + Bounded + Ord + ToStr + Eq + Clone + Copy>
     SoftBody<N, V>
{
  pub fn solve(&mut self, dt: N)
  {
    let mut constraints = ~[];

    // second order resolution
    self.collect_constraints(dt.clone(), &mut constraints, false);
  
    let res = projected_gauss_seidel_solve(constraints,
                                           self.points.len(),
                                           50,
                                           false);

    for self.points.mut_iter().enumerate().advance |(i, p)|
    { p.velocity = p.velocity + res[i].lv }

    for constraints.iter().enumerate().advance |(i, c)|
    { self.constraints[i].impulse = c.impulse.clone() }

    // first order resolution
    // self.collect_constraints(dt.clone(), &mut constraints, true);
  
    // let res = projected_gauss_seidel_solve(constraints,
    //                                        self.points.len(),
    //                                        50,
    //                                        true);

    // for self.points.mut_iter().enumerate().advance |(i, p)|
    // { p.position = p.position + res[i].lv.scalar_mul(&dt) }
  }
}
