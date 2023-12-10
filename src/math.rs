use crate::Point;

pub fn in_circle(p: Point, a: Point, b: Point, c: Point) -> bool {
  let d = a - p;
  let e = b - p;
  let f = c - p;

  let ap = d.length_squared();
  let bp = e.length_squared();
  let cp = f.length_squared();

  let res =
    d.x * (e.y * cp - bp * f.y) - d.y * (e.x * cp - bp * f.x) + ap * (e.x * f.y - e.y * f.x);

  res < 0.0
}

pub fn circumdelta(a: Point, b: Point, c: Point) -> Option<Point> {
  let d = b - a;
  let e = c - a;

  let bl = d.length_squared();
  let cl = e.length_squared();
  let det = Point::perp_dot(d, e);

  let x = (e.y * bl - d.y * cl) * (0.5 / det);
  let y = (d.x * cl - e.x * bl) * (0.5 / det);

  if (bl != 0.0) && (cl != 0.0) && (det != 0.0) {
    Some(Point { x, y })
  } else {
    None
  }
}

pub fn circumradius(a: Point, b: Point, c: Point) -> f64 {
  match circumdelta(a, b, c) {
    Some(delta) => delta.length_squared(),
    None => f64::MAX
  }
}

pub fn circumcenter(a: Point, b: Point, c: Point) -> Option<Point> {
  circumdelta(a, b, c).map(|point| (a + point))
}

pub fn counter_clockwise(p0: Point, p1: Point, p2: Point) -> bool {
  let v0 = p1 - p0;
  let v1 = p2 - p0;
  let det = Point::perp_dot(v0, v1);
  let dist = v0.length_squared() + v1.length_squared();

  if det == 0.0 {
    return false;
  }

  let reldet = (dist / det).abs();

  if reldet > 1e14 {
    return false;
  }

  det > 0.0
}

pub fn calculate_bbox(points: &[Point]) -> Option<(Point, Point)> {
  points.iter().copied().fold(None, |acc, point| match acc {
    Some((min, max)) => Some((point.min(min), point.max(max))),
    None => Some((point, point))
  })
}

pub fn calculate_bbox_center(points: &[Point]) -> (Point, f64) {
  let (min, max) = calculate_bbox(points).unwrap_or((Point::NEG_INFINITY, Point::INFINITY));
  let width = max.x - min.x;
  let height = max.y - min.y;
  let span = width * width + height * height;
  (((min + max) / 2.0), span)
}

pub fn find_closest_point(points: &[Point], p: Point) -> Option<usize> {
  let mut min_dist = f64::MAX;
  let mut min_index = None;

  for (i, &q) in points.iter().enumerate() {
    if Some(i) != min_index {
      let d = Point::distance_squared(p, q);

      if d < min_dist && d > 0.0 {
        min_index = Some(i);
        min_dist = d;
      }
    }
  }

  min_index
}

pub fn inside(p: Point, p1: Point, p2: Point) -> bool {
  (p2.y - p1.y) * p.x + (p1.x - p2.x) * p.y + Point::perp_dot(p2, p1) < 0.0
}

pub fn intersection(cp1: Point, cp2: Point, s: Point, e: Point) -> Point {
  let dc = cp1 - cp2;
  let dp = s - e;

  let n1 = Point::perp_dot(cp1, cp2);
  let n2 = Point::perp_dot(s, e);
  let n3 = Point::perp_dot(dc, dp).recip();

  Point::new((n1 * dp.x - n2 * dc.x) * n3, (n1 * dp.y - n2 * dc.y) * n3)
}
