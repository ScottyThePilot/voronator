//! Implements the Delaunay triangulation algorithm.
//!
//! This module was ported from the original [Delaunator](https://github.com/mapbox/delaunator), by Mapbox.
//! If a triangulation is possible a given set of points in the 2D space, it returns a [`Triangulation`] structure.
//! This structure contains three main components: [`triangles`], [`halfedges`] and [`hull`]:
//! ```no_run
//! # use voronator::Point;
//! # use voronator::delaunator::Triangulation;
//! let coords = vec![Point::new(0.0, 0.0), Point::new(1.0, 0.0), Point::new(1.0, 1.0), Point::new(0.0, 1.0)];
//! let delaunay = Triangulation::new(&coords).unwrap();
//! ```
//! - `triangles`: A `Vec<usize>` that contains the indices for each vertex of a triangle in the original array.
//!   All triangles are directed counter-clockwise. To get the coordinates of all triangles, use:
//!   ```no_run
//!   # use voronator::Point;
//!   # use voronator::delaunator::Triangulation;
//!   # let coords: Vec<Point> = Vec::new();
//!   # let delaunay = Triangulation::new(&coords).unwrap();
//!   let mut t = 0;
//!   loop {
//!     println!("[{:?}, {:?}, {:?}]",
//!       coords[delaunay.triangles[t]],
//!       coords[delaunay.triangles[t + 1]],
//!       coords[delaunay.triangles[t + 2]]
//!     );
//!     t += 3;
//!   }
//!   ```
//! - `halfedges`: `Vec<usize>` array of triangle half-edge indices that allows you to traverse the triangulation.
//!   i-th half-edge in the array corresponds to vertex `triangles[i]` the half-edge is coming from. `halfedges[i]`
//!   is the index of a twin half-edge in an adjacent triangle (or `INVALID_INDEX` for outer half-edges on the convex
//!   hull). The flat array-based data structures might be counterintuitive, but they're one of the key reasons this
//!   library is fast.
//! - `hull`: A `Vec<usize>` array of indices that reference points on the convex hull of the input data,
//!   counter-clockwise.
//!
//! The last two components, `inedges` and `outedges`, are for voronator internal use only.
//!
//! # Example
//!
//! ```
//! extern crate voronator;
//!
//! use voronator::delaunator::Triangulation;
//! use voronator::Point;
//!
//! fn main() {
//!   let points = vec![
//!     Point::new(0.0, 0.0),
//!     Point::new(1.0, 0.0),
//!     Point::new(1.0, 1.0),
//!     Point::new(0.0, 1.0),
//!   ];
//!
//!   let t = Triangulation::new(&points).expect("No triangulation exists for this input.");
//!
//!   for i in 0..t.triangles_count() {
//!     let i0 = t.triangles[3 * i];
//!     let i1 = t.triangles[3 * i + 1];
//!     let i2 = t.triangles[3 * i + 2];
//!
//!     let p = [points[i0], points[i1], points[i2]];
//!
//!     println!("triangle {}: {:?}", i, p);
//!   }
//! }
//! ```
//!
//! [`Triangulation`]: ./struct.Triangulation.html
//! [`triangles`]: ./struct.Triangulation.html#structfield.triangles
//! [`halfedges`]: ./struct.Triangulation.html#structfield.halfedges
//! [`hull`]: ./struct.Triangulation.html#structfield.hull

use maybe_parallel_iterator::IntoMaybeParallelRefIterator;

use crate::Point;

/// Defines a comparison epsilon used for floating-point comparisons
pub const EPSILON: f64 = f64::EPSILON * 2.0;

/// Defines an invalid index in the Triangulation vectors
pub const INVALID_INDEX: usize = usize::max_value();

fn equals_with_span(p: Point, q: Point, span: f64) -> bool {
  let dist = Point::distance_squared(p, q) / span;
  dist < 1e-20 // dunno about this
}

/// Returs the next halfedge for a given halfedge
///
/// # Arguments
///
/// * `i` - The current halfedge index
pub fn next_halfedge(i: usize) -> usize {
  if i % 3 == 2 {
    i - 2
  } else {
    i + 1
  }
}

/// Returs the previous halfedge for a given halfedge
///
/// # Arguments
///
/// * `i` - The current halfedge index
pub fn prev_halfedge(i: usize) -> usize {
  if i % 3 == 0 {
    i + 2
  } else {
    i - 1
  }
}

/// Returns a vec containing indices for the 3 edges of a triangle t
///
/// # Arguments
///
/// * `t` - The triangle index
pub fn edges_of_triangle(t: usize) -> [usize; 3] {
  [3 * t, 3 * t + 1, 3 * t + 2]
}

/// Returns the triangle associated with the given edge
///
/// # Arguments
///
/// * `e` - The edge index
pub fn triangle_of_edge(e: usize) -> usize {
  ((e as f64) / 3.0).floor() as usize
}

/// Returns a vec containing the indices of the corners of the given triangle
///
/// # Arguments
///
/// * `t` - The triangle index
/// * `delaunay` - A reference to a fully constructed Triangulation
pub fn points_of_triangle(t: usize, delaunay: &Triangulation) -> [usize; 3] {
  edges_of_triangle(t).map(|e| delaunay.triangles[e])
}

/// Returns a vec containing the indices for the adjacent triangles of the given triangle
///
/// # Arguments
///
/// * `t` - The triangle index
/// * `delaunay` - A reference to a fully constructed Triangulation
pub fn triangles_adjacent_to_triangle(t: usize, delaunay: &Triangulation) -> Vec<usize> {
  edges_of_triangle(t)
    .into_iter()
    .filter_map(|e| {
      let opposite = delaunay.halfedges[e];
      if opposite != INVALID_INDEX {
        Some(triangle_of_edge(opposite))
      } else {
        None
      }
    })
    .collect()
}

/// Returns a vec containing all edges around a point
///
/// # Arguments
///
/// * `start` - The start point index
/// * `delaunay` - A reference to a fully constructed Triangulation
pub fn edges_around_point(start: usize, delaunay: &Triangulation) -> Vec<usize> {
  let mut result: Vec<usize> = Vec::new();

  // If the starting index is invalid we can't continue
  if start == INVALID_INDEX {
    return result;
  }

  let mut incoming = start;
  loop {
    result.push(incoming);
    let outgoing = next_halfedge(incoming);
    incoming = delaunay.halfedges[outgoing];
    if incoming == INVALID_INDEX || incoming == start {
      break;
    }
  }

  result
}

/// Represents a Delaunay triangulation for a given set of points. See example in [`delaunator`] for usage details.
///
/// [`delaunator`]: ./index.html#example

#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug)]
pub struct Triangulation {
  /// Contains the indices for each vertex of a triangle in the original array. All triangles are directed counter-clockwise.
  pub triangles: Vec<usize>,
  /// A `Vec<usize>` of triangle half-edge indices that allows you to traverse the triangulation. i-th half-edge in the array corresponds to vertex `triangles[i]` the half-edge is coming from. `halfedges[i]` is the index of a twin half-edge in an adjacent triangle (or `INVALID_INDEX` for outer half-edges on the convex hull).
  pub halfedges: Vec<usize>,
  /// A `Vec<usize>` array of indices that reference points on the convex hull of the input data, counter-clockwise.
  pub hull: Vec<usize>,
  /// A `Vec<usize>` that contains indices for halfedges of points in the hull that points inwards to the diagram. Only for [`voronator`] internal use.
  ///
  /// [`voronator`]: ../index.html
  pub inedges: Vec<usize>,
  /// A `Vec<usize>` that contains indices for halfedges of points in the hull that points outwards to the diagram. Only for [`voronator`] internal use.
  ///
  /// [`voronator`]: ../index.html
  pub outedges: Vec<usize>
}

impl Triangulation {
  /// Calculates the Delaunay triangulation, if it exists, for a given set of 2D points.
  pub fn new(points: &[Point]) -> Option<Self> {
    triangulate(points)
  }

  fn empty(n: usize) -> Self {
    let max_triangles = 2 * n - 5;
    Self {
      triangles: Vec::with_capacity(max_triangles * 3),
      halfedges: Vec::with_capacity(max_triangles * 3),
      hull: Vec::new(),
      inedges: vec![INVALID_INDEX; n],
      outedges: vec![INVALID_INDEX; n]
    }
  }

  /// Returns the number of triangles calculated in the triangulation. Same as `triangles.len() / 3`.
  pub fn triangles_count(&self) -> usize {
    self.triangles.len() / 3
  }

  fn legalize(&mut self, p: usize, points: &[Point], hull: &mut Hull) -> usize {
    /* if the pair of triangles doesn't satisfy the Delaunay condition
     * (p1 is inside the circumcircle of [p0, pl, pr]), flip them,
     * then do the same check/flip recursively for the new pair of triangles
     *
     *           pl                    pl
     *          /||\                  /  \
     *       al/ || \bl            al/    \a
     *        /  ||  \              /      \
     *       /  a||b  \    flip    /___ar___\
     *     p0\   ||   /p1   =>   p0\---bl---/p1
     *        \  ||  /              \      /
     *       ar\ || /br             b\    /br
     *          \||/                  \  /
     *           pr                    pr
     */
    let mut i: usize = 0;
    let mut ar;
    let mut a = p;

    let mut edge_stack: Vec<usize> = Vec::new();

    loop {
      let b = self.halfedges[a];
      ar = prev_halfedge(a);

      if b == INVALID_INDEX {
        if i > 0 {
          i -= 1;
          a = edge_stack[i];
          continue;
        } else {
          break;
        }
      }

      let al = next_halfedge(a);
      let bl = prev_halfedge(b);

      let p0 = self.triangles[ar];
      let pr = self.triangles[a];
      let pl = self.triangles[al];
      let p1 = self.triangles[bl];

      let illegal = crate::math::in_circle(points[p1], points[p0], points[pr], points[pl]);
      if illegal {
        self.triangles[a] = p1;
        self.triangles[b] = p0;

        let hbl = self.halfedges[bl];

        // Edge swapped on the other side of the hull (rare).
        // Fix the halfedge reference
        if hbl == INVALID_INDEX {
          let mut e = hull.start;
          loop {
            if hull.tri[e] == bl {
              hull.tri[e] = a;
              break;
            }

            e = hull.prev[e];

            if e == hull.start {
              break;
            }
          }
        }

        self.link(a, hbl);
        self.link(b, self.halfedges[ar]);
        self.link(ar, bl);

        let br = next_halfedge(b);

        if i < edge_stack.len() {
          edge_stack[i] = br;
        } else {
          edge_stack.push(br);
        }

        i += 1;
      } else if i > 0 {
        i -= 1;
        a = edge_stack[i];
        continue;
      } else {
        break;
      }
    }

    ar
  }

  fn link(&mut self, a: usize, b: usize) {
    let s: usize = self.halfedges.len();

    if a == s {
      self.halfedges.push(b);
    } else if a < s {
      self.halfedges[a] = b;
    } else {
      // todo: fix hard error, make it recoverable or graceful
      panic!("Cannot link edge")
    }

    if b != INVALID_INDEX {
      let s2: usize = self.halfedges.len();
      if b == s2 {
        self.halfedges.push(a);
      } else if b < s2 {
        self.halfedges[b] = a;
      } else {
        // todo: fix hard error, make it recoverable or graceful
        panic!("Cannot link edge")
      }
    }
  }

  fn add_triangle(&mut self, i: [usize; 3], x: [usize; 3]) -> usize {
    let t: usize = self.triangles.len();
    self.triangles.extend(i);
    self.link(t, x[0]);
    self.link(t + 1, x[1]);
    self.link(t + 2, x[2]);

    t
  }
}

//@see https://stackoverflow.com/questions/33333363/built-in-mod-vs-custom-mod-function-improve-the-performance-of-modulus-op/33333636#33333636
fn fast_mod(i: usize, c: usize) -> usize {
  if i >= c {
    i % c
  } else {
    i
  }
}

// monotonically increases with real angle,
// but doesn't need expensive trigonometry
fn pseudo_angle(d: Point) -> f64 {
  let p = d.x / (d.x.abs() + d.y.abs());
  if d.y > 0.0 {
    (3.0 - p) / 4.0
  } else {
    (1.0 + p) / 4.0
  }
}

#[derive(Debug, Clone, PartialEq)]
struct Hull {
  prev: Vec<usize>,
  next: Vec<usize>,
  tri: Vec<usize>,
  hash: Vec<usize>,
  start: usize,
  center: Point
}

impl Hull {
  fn new(n: usize, center: Point, i: [usize; 3], points: &[Point]) -> Self {
    // initialize a hash table for storing edges of the advancing convex hull
    let hash_len = (n as f64).sqrt().ceil() as usize;

    let mut hull = Self {
      prev: vec![0; n],
      next: vec![0; n],
      tri: vec![0; n],
      hash: vec![INVALID_INDEX; hash_len],
      start: i[0],
      center
    };

    hull.next[i[0]] = i[1];
    hull.prev[i[2]] = i[1];
    hull.next[i[1]] = i[2];
    hull.prev[i[0]] = i[2];
    hull.next[i[2]] = i[0];
    hull.prev[i[1]] = i[0];

    hull.tri[i[0]] = 0;
    hull.tri[i[1]] = 1;
    hull.tri[i[2]] = 2;

    // todo here

    hull.hash_edge(points[i[0]], i[0]);
    hull.hash_edge(points[i[1]], i[1]);
    hull.hash_edge(points[i[2]], i[2]);

    hull
  }

  fn hash_key(&self, p: Point) -> usize {
    let d = p - self.center;
    let angle: f64 = pseudo_angle(d);
    let len = self.hash.len();

    fast_mod((angle * (len as f64)).floor() as usize, len)
  }

  fn hash_edge(&mut self, p: Point, i: usize) {
    let key = self.hash_key(p);
    self.hash[key] = i;
  }

  fn find_visible_edge(&self, p: Point, span: f64, points: &[Point]) -> (usize, bool) {
    // find a visible edge on the convex hull using edge hash
    let mut start = 0;
    let key = self.hash_key(p);
    for j in 0..self.hash.len() {
      let index = fast_mod(key + j, self.hash.len());
      start = self.hash[index];
      if start != INVALID_INDEX && start != self.next[start] {
        break;
      }
    }

    // Make sure what we found is on the hull.
    // todo: return something that represents failure to fail gracefully instead
    if self.prev[start] == start || self.prev[start] == INVALID_INDEX {
      panic!("not in the hull");
    }

    start = self.prev[start];
    let mut e = start;
    let mut q: usize;

    // Advance until we find a place in the hull where our current point
    // can be added.
    loop {
      q = self.next[e];

      if equals_with_span(p, points[e], span) || equals_with_span(p, points[q], span) {
        e = INVALID_INDEX;
        break;
      }

      if crate::math::counter_clockwise(p, points[e], points[q]) {
        break;
      }

      e = q;
      if e == start {
        e = INVALID_INDEX;
        break;
      }
    }

    (e, e == start)
  }
}

fn find_seed_triangle(center: Point, points: &[Point]) -> Option<[usize; 3]> {
  let i0 = crate::math::find_closest_point(points, center).unwrap_or(INVALID_INDEX);
  let p0 = points[i0];

  let mut i1 = crate::math::find_closest_point(points, p0).unwrap_or(INVALID_INDEX);
  let p1 = points[i1];

  // find the third point which forms the smallest circumcircle
  // with the first two
  let mut min_radius = f64::MAX;
  let mut i2 = INVALID_INDEX;
  for (i, &p) in points.iter().enumerate() {
    if i != i0 && i != i1 {
      let r = crate::math::circumradius(p0, p1, p);

      if r < min_radius {
        i2 = i;
        min_radius = r;
      }
    }
  }

  if min_radius == f64::MAX {
    None
  } else {
    let p2 = points[i2];

    if crate::math::counter_clockwise(p0, p1, p2) {
      std::mem::swap(&mut i1, &mut i2)
    }

    Some([i0, i1, i2])
  }
}

fn triangulate(points: &[Point]) -> Option<Triangulation> {
  if points.len() < 3 {
    return None;
  }

  let (center_bbox, span) = crate::math::calculate_bbox_center(points);
  let seed_triangle @ [i0, i1, i2] = find_seed_triangle(center_bbox, points)?;

  let center = crate::math::circumcenter(points[i0], points[i1], points[i2]).unwrap();

  // Calculate the distances from the center once to avoid having to
  // calculate for each compare.
  let mut dists: Vec<(usize, f64)> = points
    .maybe_par_iter()
    .enumerate()
    .map(|(i, &p)| (i, Point::distance_squared(p, center)))
    .collect();

  // sort the points by distance from the seed triangle circumcenter
  dists.sort_unstable_by(|(_, a), (_, b)| f64::total_cmp(a, b));

  let mut hull = Hull::new(points.len(), center, seed_triangle, points);

  let mut triangulation = Triangulation::empty(points.len());
  triangulation.add_triangle(seed_triangle, [INVALID_INDEX; 3]);

  let mut pp = Point::new(f64::NAN, f64::NAN);

  // go through points based on distance from the center.
  for (k, (i, _)) in dists.iter().copied().enumerate() {
    let p = points[i];

    // skip near-duplicate points
    if k > 0 && p.abs_diff_eq(pp, EPSILON) {
      continue;
    }

    // skip seed triangle points
    if seed_triangle.contains(&i) {
      continue;
    }

    pp = p;

    let (mut e, backwards) = hull.find_visible_edge(p, span, points);
    if e == INVALID_INDEX {
      continue;
    }

    // add the first triangle from the point
    let mut t = triangulation.add_triangle([e, i, hull.next[e]], [
      INVALID_INDEX,
      INVALID_INDEX,
      hull.tri[e]
    ]);

    hull.tri[i] = triangulation.legalize(t + 2, points, &mut hull);
    hull.tri[e] = t;

    // walk forward through the hull, adding more triangles and
    // flipping recursively
    let mut next = hull.next[e];
    loop {
      let q = hull.next[next];
      if !crate::math::counter_clockwise(p, points[next], points[q]) {
        break;
      }
      t = triangulation.add_triangle([next, i, q], [hull.tri[i], INVALID_INDEX, hull.tri[next]]);

      hull.tri[i] = triangulation.legalize(t + 2, points, &mut hull);
      hull.next[next] = next;
      next = q;
    }

    // walk backward from the other side, adding more triangles
    // and flipping
    if backwards {
      loop {
        let q = hull.prev[e];
        if !crate::math::counter_clockwise(p, points[q], points[e]) {
          break;
        }
        t = triangulation.add_triangle([q, i, e], [INVALID_INDEX, hull.tri[e], hull.tri[q]]);
        triangulation.legalize(t + 2, points, &mut hull);
        hull.tri[q] = t;
        hull.next[e] = e;
        e = q;
      }
    }

    // update the hull indices
    hull.prev[i] = e;
    hull.next[e] = i;
    hull.prev[next] = i;
    hull.next[i] = next;
    hull.start = e;

    hull.hash_edge(p, i);
    hull.hash_edge(points[e], e);
  }

  for e in 0..triangulation.triangles.len() {
    let endpoint = triangulation.triangles[next_halfedge(e)];
    if triangulation.halfedges[e] == INVALID_INDEX
      || triangulation.inedges[endpoint] == INVALID_INDEX
    {
      triangulation.inedges[endpoint] = e;
    }
  }

  let mut vert0: usize;
  let mut vert1 = hull.start;
  loop {
    vert0 = vert1;
    vert1 = hull.next[vert1];
    triangulation.inedges[vert1] = hull.tri[vert0];
    triangulation.outedges[vert0] = hull.tri[vert1];
    if vert1 == hull.start {
      break;
    }
  }

  let mut e = hull.start;
  loop {
    triangulation.hull.push(e);
    e = hull.next[e];
    if e == hull.start {
      break;
    }
  }

  Some(triangulation)
}
