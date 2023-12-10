extern crate plotters;
extern crate rand;
extern crate voronator;

use plotters::prelude::*;
use rand::prelude::*;
use voronator::{Diagram, Point};

const IMG_WIDTH: u32 = 500;
const IMG_HEIGHT: u32 = 500;

fn get_points(n: i32, jitter: f64) -> Vec<Point> {
  let mut rng = rand::thread_rng();
  let mut points: Vec<Point> = Vec::new();
  for i in 0..n + 1 {
    for j in 0..n + 1 {
      points.push(Point {
        x: (i as f64) + jitter * (rng.gen::<f64>() - rng.gen::<f64>()),
        y: (j as f64) + jitter * (rng.gen::<f64>() - rng.gen::<f64>())
      });
    }
  }

  points
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
  let size = 10;
  let points: Vec<Point> = get_points(size, 0.6)
    .into_iter()
    .map(|p| Point {
      x: ((IMG_WIDTH as f64) / 20.0 + p.x * (IMG_WIDTH as f64)) / (size as f64),
      y: ((IMG_HEIGHT as f64) / 20.0 + p.y * (IMG_HEIGHT as f64)) / (size as f64)
    })
    .collect();

  let now = std::time::Instant::now();
  let diagram = Diagram::new_voronoi(
    points.to_vec(),
    Point::new(0.0, 0.0),
    Point::new(IMG_WIDTH as f64, IMG_HEIGHT as f64)
  )
  .expect("No triangulation exists for this input.");

  println!(
    "time it took to generating a diagram for {} points: {}ms",
    points.len(),
    now.elapsed().as_millis()
  );

  let root = BitMapBackend::new("plot.png", (IMG_WIDTH, IMG_HEIGHT)).into_drawing_area();
  root.fill(&WHITE)?;

  let root = root.apply_coord_spec(RangedCoord::<RangedCoordf32, RangedCoordf32>::new(
    0f32..IMG_WIDTH as f32,
    0f32..IMG_HEIGHT as f32,
    (0..IMG_WIDTH as i32, 0..IMG_HEIGHT as i32)
  ));

  println!("triangles: {}", diagram.delaunay.triangles_count());
  println!("cells: {}", diagram.cells.len());

  for cell in diagram.cells.iter() {
    let p: Vec<(f32, f32)> = cell
      .points
      .iter()
      .map(|x| (x.x as f32, x.y as f32))
      .collect();

    for _ in 0..p.len() {
      let plot = PathElement::new(p.clone(), ShapeStyle {
        color: BLACK.to_rgba(),
        filled: true,
        stroke_width: 1
      });
      root.draw(&plot)?;
    }
  }

  Ok(())
}
