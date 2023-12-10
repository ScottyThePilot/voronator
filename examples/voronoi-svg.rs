extern crate voronator;

use svg::node::element::Path;
use svg::{node, Document, Node};
use voronator::{Diagram, Point};

fn main() {
  let points = vec![
    Point::new(2520.0, 856.0),
    Point::new(794.0, 66.0),
    Point::new(974.0, 446.0),
  ];

  let voronoi = Diagram::new_voronoi(
    points.clone(),
    Point::new(0.0, 0.0),
    Point::new(2560.0, 2560.0)
  )
  .unwrap();

  let mut document = Document::new().set("viewBox", (0, 0, 2560, 2560));
  let colours = ["blue", "green", "red"];
  let mut i = 0;

  for polygon in voronoi.cells.iter() {
    let cell_vertices = polygon.points.as_slice();

    let mut is_start = true;
    let mut d: Vec<node::element::path::Command> = Vec::new();
    for point in cell_vertices.into_iter() {
      if is_start {
        d.push(node::element::path::Command::Move(
          node::element::path::Position::Absolute,
          (point.x, point.y).into()
        ));
        is_start = false;
      } else {
        d.push(node::element::path::Command::Line(
          node::element::path::Position::Absolute,
          (point.x, point.y).into()
        ));
      }
    }
    d.push(node::element::path::Command::Close);
    let data = node::element::path::Data::from(d);

    let path = Path::new()
      .set("fill", colours[i])
      .set("stroke", "black")
      .set("stroke-width", 1)
      .set("d", data);

    document.append(path);

    i += 1;
  }

  for point in points {
    document.append(
      node::element::Rectangle::new()
        .set("x", point.x - 15.0)
        .set("y", point.y - 15.0)
        .set("width", 30)
        .set("height", 30)
    );
  }

  svg::save("example.svg", &document).unwrap();
}
