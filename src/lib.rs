//! Constructs a Voronoi diagram given a set of points.
//!
//! This module was adapted from [d3-delaunay](https://github.com/d3/d3-delaunay) and from
//! [Red Blob Games](https://www.redblobgames.com/x/2022-voronoi-maps-tutorial/) Voronoi maps tutorial.
//! It implements the Delaunay triangulation dual extraction, which is the
//! Voronoi diagram. It also implements a centroidal tesselation based on the
//! Voronoi diagram, but using centroids instead of circumcenters for the
//! vertices of the cell polygons.
//!
//! Apart from the triangle center they are using, the Voronoi and Centroidal
//! diagrams differ in how they handle the hull cells. The Voronoi diagram
//! implements a clipping algorithm that clips the diagram into a bounding box,
//! thus extracting neat polygons around the hull. The Centroid diagram, in the
//! other hand, doesn't. The outer cells can be missing or be distorted,
//! as triangles calculated by the Delaunay triangulation can be too thin in the
//! hull, causing centroid calculation to be "bad".
//!
//! If you have a robust solution for this particular problem, please let me
//! know either by creating an issue or through a pull-request, and I will make
//! sure to add your solution with the proper credits.
//!
//! # Example
//!
//! ## Voronoi Diagram
//! ```rust
//! extern crate voronator;
//! extern crate rand;
//!
//! use rand::distributions::Uniform;
//! use rand::prelude::*;
//! use voronator::{Diagram, Point};
//!
//! fn main() {
//!   let mut rng = rand::thread_rng();
//!   let range1 = Uniform::new(0.0, 100.0);
//!   let range2 = Uniform::new(0.0, 100.0);
//!   let mut points: Vec<Point> = (0..10)
//!     .map(|_| Point::new(rng.sample(&range1), rng.sample(&range2)))
//!     .collect();
//!
//!   let min = Point::new(0.0, 0.0);
//!   let max = Point::new(100.0, 100.0);
//!   let diagram = Diagram::new_voronoi(points, min, max).unwrap();
//!
//!   for cell in diagram.cells.iter() {
//!     let p = cell.points.iter().copied().collect::<Vec<Point>>();
//!
//!     println!("{:?}", p);
//!   }
//! }
//! ```
//!
//! ## Centroidal Tesselation Diagram
//! ```rust
//! extern crate voronator;
//! extern crate rand;
//!
//! use rand::distributions::Uniform;
//! use rand::prelude::*;
//! use voronator::{Diagram, Point};
//!
//! fn main() {
//!   let mut rng = rand::thread_rng();
//!   let range1 = Uniform::new(0.0, 100.0);
//!   let range2 = Uniform::new(0.0, 100.0);
//!   let mut points: Vec<Point> = (0..10)
//!     .map(|_| Point::new(rng.sample(&range1), rng.sample(&range2)))
//!     .collect();
//!
//!   let diagram = Diagram::new_centroidal(points).unwrap();
//!
//!   for cell in diagram.cells {
//!     let p = cell.points.iter().copied().collect::<Vec<Point>>();
//!
//!     println!("{:?}", p);
//!   }
//! }
//! ```

#![warn(
  future_incompatible,
  missing_copy_implementations,
  missing_debug_implementations,
  missing_docs
)]

pub extern crate glam;

pub mod delaunator;
mod math;
pub mod polygon;

#[doc(no_inline)]
pub use glam::DVec2 as Point;
use maybe_parallel_iterator::IntoMaybeParallelIterator;

use crate::delaunator::{Triangulation, INVALID_INDEX};
use crate::polygon::Polygon;

/// Can represent a Voronoi diagram or centroidal tesselation diagram.
#[derive(Debug, Clone, PartialEq)]
pub struct Diagram {
  /// A list of input site points.
  pub sites: Vec<Point>,
  /// A [`Triangulation`] that contains the Delaunay triangulation information.
  pub delaunay: Triangulation,
  /// A list containing the circumcenters of each triangle in the Delaunay
  /// triangulation.
  pub centers: Vec<Point>,
  /// A list of [`Polygon`]s that make up this Voronoi diagram.
  pub cells: Vec<Polygon>,
  /// Contains a list of cell neighbors for each site/cell.
  pub neighbors: Vec<Vec<usize>>
}

impl Diagram {
  /// Creates a centroidal tesselation, if it exists, for a given set of points.
  pub fn new_centroidal(sites: Vec<Point>) -> Option<Self> {
    let delaunay = Triangulation::new(&sites)?;
    let centers = calculate_centroids(&sites, &delaunay);
    let cells = calculate_polygons(&sites, &centers, &delaunay, None);
    let neighbors = calculate_neighbors(&sites, &delaunay);
    Some(Diagram {
      sites,
      delaunay,
      centers,
      cells,
      neighbors
    })
  }

  /// Creates a Voronoi diagram, if it exists, for a given set of points.
  pub fn new_voronoi(sites: Vec<Point>, min: Point, max: Point) -> Option<Self> {
    // Create a polygon defining the region to clip to (rectangle from min point to max point)
    let clip_polygon = Polygon::new(vec![
      Point::new(min.x, min.y),
      Point::new(max.x, min.y),
      Point::new(max.x, max.y),
      Point::new(min.x, max.y),
    ]);

    Diagram::with_bounding_polygon_voronoi(sites, &clip_polygon)
  }

  /// Creates a Voronoi diagram, if it exists, for a given set of points bounded by the supplied polygon.
  pub fn with_bounding_polygon_voronoi(
    mut sites: Vec<Point>,
    clip_polygon: &Polygon
  ) -> Option<Self> {
    sites.extend(helper_points(clip_polygon));

    let delaunay = Triangulation::new(&sites)?;
    let centers = calculate_circumcenters(&sites, &delaunay);
    let mut cells = calculate_polygons(&sites, &centers, &delaunay, Some(clip_polygon));
    let mut neighbors = calculate_neighbors(&sites, &delaunay);

    sites.truncate(sites.len() - NUM_HELPER_POINTS);
    cells.truncate(cells.len() - NUM_HELPER_POINTS);
    neighbors.truncate(neighbors.len() - NUM_HELPER_POINTS);

    Some(Diagram {
      sites,
      delaunay,
      centers,
      cells,
      neighbors
    })
  }
}

fn calculate_polygons(
  points: &[Point],
  centers: &[Point],
  delaunay: &Triangulation,
  clip_polygon: Option<&Polygon>
) -> Vec<Polygon> {
  (0..points.len())
    .into_maybe_par_iter()
    .map(|t| {
      let incoming = delaunay.inedges[t];
      let edges = crate::delaunator::edges_around_point(incoming, delaunay);
      let triangles: Vec<usize> = edges
        .into_iter()
        .map(crate::delaunator::triangle_of_edge)
        .collect();
      let polygon = Polygon::new(triangles.into_iter().map(|t| centers[t]).collect());
      if let Some(clip_polygon) = clip_polygon {
        polygon.clip(clip_polygon)
      } else {
        polygon
      }
    })
    .collect()
}

fn calculate_centroids(points: &[Point], delaunay: &Triangulation) -> Vec<Point> {
  (0..delaunay.triangles_count())
    .into_maybe_par_iter()
    .map(|t| {
      let mut sum = Point::ZERO;
      for i in 0..3 {
        let s = 3 * t + i; // triangle coord index
        sum += points[delaunay.triangles[s]];
      }

      sum / 3.0
    })
    .collect()
}

fn calculate_circumcenters(points: &[Point], delaunay: &Triangulation) -> Vec<Point> {
  (0..delaunay.triangles_count())
    .into_maybe_par_iter()
    .map(|t| {
      let v = crate::delaunator::points_of_triangle(t, delaunay).map(|p| points[p]);
      crate::math::circumcenter(v[0], v[1], v[2]).unwrap_or(Point::ZERO)
    })
    .collect()
}

fn calculate_neighbors(points: &[Point], delaunay: &Triangulation) -> Vec<Vec<usize>> {
  (0..points.len())
    .into_maybe_par_iter()
    .map(|t| {
      let mut neighbours: Vec<usize> = Vec::new();

      let e0 = delaunay.inedges[t];
      if e0 != INVALID_INDEX {
        let mut e = e0;
        loop {
          neighbours.push(delaunay.triangles[e]);
          e = crate::delaunator::next_halfedge(e);
          if delaunay.triangles[e] != t {
            break;
          }
          e = delaunay.halfedges[e];
          if e == INVALID_INDEX {
            neighbours.push(delaunay.triangles[delaunay.outedges[t]]);
            break;
          }
          if e == e0 {
            break;
          }
        }
      }

      neighbours
    })
    .collect()
}

const NUM_HELPER_POINTS: usize = 4;

fn helper_points(polygon: &Polygon) -> [Point; NUM_HELPER_POINTS] {
  let (min, max) = crate::math::calculate_bbox(&polygon.points).unwrap_or((Point::MAX, Point::MIN));
  let size = max - min;
  let center = (max + min) / 2.0;

  [
    Point::new(min.x - size.x * 2.0, center.y),
    Point::new(max.x + size.x * 2.0, center.y),
    Point::new(center.x, min.y - size.y * 2.0),
    Point::new(center.x, max.y + size.y * 2.0)
  ]
}
