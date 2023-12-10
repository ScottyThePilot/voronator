extern crate serde_json;
extern crate voronator;

use std::f64;
use std::fs::File;

use voronator::delaunator::*;
use voronator::{Diagram, Point};

#[test]
fn basic() {
  let points = load_fixture("tests/fixtures/ukraine.json");

  validate(&points);
}

#[test]
fn square() {
  let points = vec![
    Point::new(0.0, 0.0),
    Point::new(1.0, 0.0),
    Point::new(1.0, 1.0),
    Point::new(0.0, 1.0),
  ];

  validate(&points);
}

#[test]
fn issue_11() {
  let points: Vec<Point> = vec![
    Point::new(516.0, 661.0),
    Point::new(369.0, 793.0),
    Point::new(426.0, 539.0),
    Point::new(273.0, 525.0),
    Point::new(204.0, 694.0),
    Point::new(747.0, 750.0),
    Point::new(454.0, 390.0),
  ];

  validate(&points);
}

#[test]
fn issue_13() {
  let points = load_fixture("tests/fixtures/issue13.json");

  validate(&points);
}

#[test]
fn duplicated_points() {
  use std::collections::HashSet;

  let points = [
    Point::new(2520.0, 856.0),
    Point::new(794.0, 66.0),
    Point::new(974.0, 446.0)
  ];

  let voronoi = Diagram::new_voronoi(
    points.to_vec(),
    Point::new(0.0, 0.0),
    Point::new(2560.0, 2560.0)
  )
  .unwrap();

  println!("# cells: {}", voronoi.cells.len());
  for polygon in voronoi.cells.iter() {
    let cell_vertices = polygon.points.as_slice();

    let expected = cell_vertices.len();
    let actual = cell_vertices
      .iter()
      .map(|n| format!("{:?}", n))
      .collect::<HashSet<String>>()
      .len();
    println!("{} != {}", expected, actual);
    println!(
      "{}",
      cell_vertices
        .iter()
        .map(|n| format!("{:?}", n))
        .collect::<String>()
    );

    assert!(expected == actual)
  }
}

#[test]
fn issue_24() {
  let points = vec![
    Point::new(382.0, 302.0),
    Point::new(382.0, 328.0),
    Point::new(382.0, 205.0),
    Point::new(623.0, 175.0),
    Point::new(382.0, 188.0),
    Point::new(382.0, 284.0),
    Point::new(623.0, 87.0),
    Point::new(623.0, 341.0),
    Point::new(141.0, 227.0),
  ];

  validate(&points);
}

#[test]
fn issue_43() {
  let points = load_fixture("tests/fixtures/issue43.json");

  validate(&points);
}

#[test]
fn issue_44() {
  let points = load_fixture("tests/fixtures/issue44.json");

  validate(&points);
}

#[test]
fn robustness() {
  let robustness1 = load_fixture("tests/fixtures/robustness1.json");

  validate(&robustness1);
  validate(&(scale_points(&robustness1, 1e-9)));
  validate(&(scale_points(&robustness1, 1e-2)));
  validate(&(scale_points(&robustness1, 100.0)));
  validate(&(scale_points(&robustness1, 1e9)));

  let robustness2 = load_fixture("tests/fixtures/robustness2.json");
  validate(&robustness2[0..100]);
  validate(&robustness2);
}

#[test]
fn few_points() {
  let points = load_fixture("tests/fixtures/ukraine.json");

  let d = Triangulation::new(&points[0..0]);
  assert!(d.is_none());

  let d = Triangulation::new(&points[0..1]);
  assert!(d.is_none());

  let d = Triangulation::new(&points[0..2]);
  assert!(d.is_none());
}

#[test]
fn collinear() {
  let points = vec![
    Point::new(0.0, 0.0),
    Point::new(1.0, 0.0),
    Point::new(3.0, 0.0),
    Point::new(2.0, 0.0),
  ];

  let d = Triangulation::new(&points);
  assert!(d.is_none());
}

fn scale_points(points: &[Point], scale: f64) -> Vec<Point> {
  points.iter().map(|&p| p * scale).collect()
}

fn load_fixture(path: &str) -> Vec<Point> {
  let file = File::open(path).unwrap();
  let u: Vec<(f64, f64)> = serde_json::from_reader(file).unwrap();
  u.into_iter().map(Point::from).collect()
}

fn triangle_area(points: &[Point], delaunay: &Triangulation) -> f64 {
  let mut vals: Vec<f64> = Vec::new();
  let mut t = 0;
  while t < delaunay.triangles.len() {
    let a = &points[delaunay.triangles[t]];
    let b = &points[delaunay.triangles[t + 1]];
    let c = &points[delaunay.triangles[t + 2]];
    let val = ((b.y - a.y) * (c.x - b.x) - (b.x - a.x) * (c.y - b.y)).abs();
    vals.push(val);
    t += 3;
  }

  better_sum(&vals)
}

fn hull_area(points: &[Point], delaunay: &Triangulation) -> f64 {
  let mut hull_areas = Vec::new();
  let mut i = 0;
  let mut j = delaunay.hull.len() - 1;
  while i < delaunay.hull.len() {
    let p0 = &points[delaunay.hull[j]];
    let p = &points[delaunay.hull[i]];
    hull_areas.push((p.x - p0.x) * (p.y + p0.y));
    j = i;
    i += 1;
  }

  better_sum(&hull_areas)
}

fn better_sum(x: &[f64]) -> f64 {
  let mut sum = x[0];
  let mut err: f64 = 0.0;
  for i in 1..x.len() {
    let k = x[i];
    let m = sum + k;
    err += if sum.abs() >= k.abs() {
      sum - m + k
    } else {
      k - m + sum
    };
    sum = m;
  }

  sum + err
}

fn validate(points: &[Point]) {
  let triangulation = Triangulation::new(&points).expect("No triangulation exists for this input");

  // Validate halfedges
  for (i, &h) in triangulation.halfedges.iter().enumerate() {
    if h != INVALID_INDEX && triangulation.halfedges[h] != i {
      panic!("Invalid halfedge connection");
    }
  }

  // Validate triangulation
  let hull_area = hull_area(points, &triangulation);
  let triangles_area = triangle_area(points, &triangulation);

  let err = ((hull_area - triangles_area) / hull_area).abs();
  if err > EPSILON {
    panic!("Triangulation is broken: {} error", err);
  }
}
