# voronator
Port of the [`d3-delaunay`](https://github.com/d3/d3-delaunay) and [`delaunator`](https://github.com/mapbox/delaunator) libraries in Rust.

This package implements the Voronoi diagram construction as a dual of the Delaunay triangulation for a set of points. It also implements the construction of a centroidal tesselation of a Delaunay triangulation, inspired by [Red Blob Games](https://www.redblobgames.com/x/2022-voronoi-maps-tutorial/).

## Examples

```rust
extern crate voronator;
extern crate rand;

use rand::distributions::Uniform;
use rand::prelude::*;
use voronator::{Diagram, Point};

fn main() {
  let mut rng = rand::thread_rng();
  let range1 = Uniform::new(0.0, 100.0);
  let range2 = Uniform::new(0.0, 100.0);
  let mut points: Vec<Point> = (0..10)
    .map(|_| Point::new(rng.sample(&range1), rng.sample(&range2)))
    .collect();

  let min = Point::new(0.0, 0.0);
  let max = Point::new(100.0, 100.0);
  let diagram = Diagram::new_voronoi(points, min, max).unwrap();

  for cell in diagram.cells.iter() {
    let p = cell.points.iter().copied().collect::<Vec<Point>>();

    println!("{:?}", p);
  }
}
```
Possible output:

![Possible output](example.png?raw=true "Possible output")
