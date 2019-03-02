extern crate nalgebra as na;
extern crate num_complex;
#[allow(unused_variables)]
use na::{Matrix2, Vector2,VectorN};
use num_complex::Complex;
use na::linalg::FourierTransform;

fn main() {
    let f = FourierTransform::new(Vector2::new(1.,1.5));

    println!("{}",f.as_real());
}
