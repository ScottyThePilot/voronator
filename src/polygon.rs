//! Provides functions for handling polygons.
//!
//! Polygons are stored as a `Vec<Point>`.
//!
//! # Example
//!
//! ```no_run
//! extern crate voronator;
//!
//! use voronator::Point;
//! use voronator::polygon::Polygon;
//!
//! fn main() {
//!   let points = vec![Point::new(0.0, 0.0), Point::new(1.0, 0.0), Point::new(1.0, 1.0), Point::new(0.0, 1.0)];
//!   let polygon = Polygon::new(points);
//! }

use crate::Point;

/// Represents a polygon.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Polygon {
  /// The points that make up this polygon.
  pub points: Vec<Point>
}

impl Polygon {
  /// Create a polygon consisting of the points supplied.
  pub fn new(points: Vec<Point>) -> Self {
    Polygon { points }
  }

  /// Performs Sutherland-Hodgman clipping on this polygon.
  pub fn clip(self, other: &Self) -> Self {
    sutherland_hodgman(self, other)
  }
}

/// Sutherland-Hodgman clipping modified from https://rosettacode.org/wiki/Sutherland-Hodgman_polygon_clipping#C.2B.2B
fn sutherland_hodgman(subject: Polygon, clip: &Polygon) -> Polygon {
  let clip = clip.points.as_slice();
  let mut output_polygon = subject.points;

  for j in 0..clip.len() {
    let input_polygon = std::mem::take(&mut output_polygon);

    // Get clipping polygon edge
    let cp1 = clip[j];
    let cp2 = clip[(j + 1) % clip.len()];

    for i in 0..input_polygon.len() {
      // Get subject polygon edge
      let s = input_polygon[i];
      let e = input_polygon[(i + 1) % input_polygon.len()];

      let v1 = crate::math::inside(s, cp1, cp2);
      let v2 = crate::math::inside(e, cp1, cp2);

      // Case 1: Both vertices are inside:
      // Only the second vertex is added to the output list
      if v1 && v2 {
        output_polygon.push(e);

      // Case 2: First vertex is outside while second one is inside:
      // Both the point of intersection of the edge with the clip boundary
      // and the second vertex are added to the output list
      } else if !v1 && v2 {
        output_polygon.push(crate::math::intersection(cp1, cp2, s, e));
        output_polygon.push(e);

      // Case 3: First vertex is inside while second one is outside:
      // Only the point of intersection of the edge with the clip boundary
      // is added to the output list
      } else if v1 && !v2 {
        output_polygon.push(crate::math::intersection(cp1, cp2, s, e));
      }

      // Case 4: Both vertices are outside
      // No vertices are added to the output list
    }
  }

  Polygon::new(output_polygon)
}
