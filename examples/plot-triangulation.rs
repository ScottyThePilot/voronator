extern crate plotters;
extern crate rand;
extern crate voronator;

use plotters::prelude::*;
use rand::prelude::*;
use rand_distr::Uniform;
use voronator::delaunator::Triangulation;
use voronator::Point;

fn main() -> Result<(), Box<dyn std::error::Error>> {
  let mut rng = rand::thread_rng();
  let range = Uniform::new(0.0, 1000.0);
  let points = (0..1000)
    .map(|_| Point::new(rng.sample(&range), rng.sample(&range)))
    .collect::<Vec<Point>>();
  let t = Triangulation::new(&points).expect("No triangulation exists for this input.");

  let root = BitMapBackend::new("plot.png", (960, 400)).into_drawing_area();
  root.fill(&WHITE)?;

  let root = root.apply_coord_spec(RangedCoord::<RangedCoordf32, RangedCoordf32>::new(
    0f32..1000f32,
    0f32..1000f32,
    (0..1000, 0..1000)
  ));

  for i in 0..t.triangles_count() {
    let i0 = t.triangles[3 * i];
    let i1 = t.triangles[3 * i + 1];
    let i2 = t.triangles[3 * i + 2];

    let p = [
      (points[i0].x as f32, points[i0].y as f32),
      (points[i1].x as f32, points[i1].y as f32),
      (points[i2].x as f32, points[i2].y as f32)
    ];

    let color = RGBColor(rng.gen(), rng.gen(), rng.gen());
    let poly = Polygon::new(p, Into::<ShapeStyle>::into(&color));
    root.draw(&poly)?;
  }

  for &point in &points {
    let p = (point.x as f32, point.y as f32);
    root.draw(&Circle::new(p, 3, ShapeStyle::from(&BLACK).filled()))?;
  }

  Ok(())
}
