extern crate voronator;

use voronator::delaunator::Triangulation;
use voronator::Point;

fn main() {
  let points = vec![
    Point::new(0.0, 0.0),
    Point::new(1.0, 0.0),
    Point::new(1.0, 1.0),
    Point::new(0.0, 1.0),
  ];

  let triangulation = Triangulation::new(&points).unwrap();

  println!(
    "Triangulated {} points, generating {} triangles.",
    points.len(),
    triangulation.triangles.len() / 3
  );
}
