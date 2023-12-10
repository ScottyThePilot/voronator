extern crate voronator;

use std::f64;
use std::time::{Duration, Instant};

use rand::distributions::Uniform;
use rand::Rng;
use rand_distr::StandardNormal;
use voronator::delaunator::Triangulation;
use voronator::Point;

fn report(n: usize, elapsed: Duration, t: &Triangulation) {
  println!(
    "{} points ({} tris, {} hull points): {}.{}ms.",
    n,
    t.triangles_count(),
    t.hull.len(),
    elapsed.as_millis(),
    elapsed.subsec_micros()
  );
}

fn uniform(count: &[usize]) {
  let mut rng = rand::thread_rng();
  let range = Uniform::new(0.0, 1000.0);

  println!("Uniform distribution:");

  for &c in count {
    let points = (0..c)
      .map(|_| Point::new(rng.sample(&range), rng.sample(&range)))
      .collect::<Vec<Point>>();

    let now = Instant::now();
    let t = Triangulation::new(&points).expect("No triangulation exists for this input.");
    let elapsed = now.elapsed();

    report(points.len(), elapsed, &t);
  }
}

fn gaussian(count: &[usize]) {
  let mut rng = rand::thread_rng();

  println!("Gaussian distribution:");

  for &c in count {
    let points = (0..c)
      .map(|_| Point::new(rng.sample(StandardNormal), rng.sample(StandardNormal)) * 1000.0)
      .collect::<Vec<Point>>();

    let now = Instant::now();
    let t = Triangulation::new(&points).expect("No triangulation exists for this input.");
    let elapsed = now.elapsed();

    report(points.len(), elapsed, &t);
  }
}

fn grid(count: &[usize]) {
  println!("Grid distribution:");

  for &c in count {
    let size = (c as f64).sqrt().floor() as usize;
    let mut points: Vec<Point> = Vec::new();

    for i in 0..size {
      for j in 0..size {
        points.push(Point::new(i as f64, j as f64));
      }
    }

    let now = Instant::now();
    let t = Triangulation::new(&points).expect("No triangulation exists for this input.");
    let elapsed = now.elapsed();

    report(points.len(), elapsed, &t);
  }
}

fn degenerate(count: &[usize]) {
  println!("Degenerate distribution:");

  for &c in count {
    let mut points: Vec<Point> = vec![Point { x: 0.0, y: 0.0 }];
    for i in 0..c {
      let angle = 2.0 * f64::consts::PI * (i as f64) / (c as f64);
      points.push(Point::new(1e10 * angle.sin(), 1e10 * angle.cos()));
    }

    let now = Instant::now();
    let t = Triangulation::new(&points).expect("No triangulation exists for this input.");
    let elapsed = now.elapsed();

    report(points.len(), elapsed, &t);
  }
}

fn main() {
  let count = [20000, 100000, 200000, 500000, 1000000];

  gaussian(&count);
  uniform(&count);
  grid(&count);
  degenerate(&count);
}
