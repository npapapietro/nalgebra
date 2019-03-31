extern crate nalgebra as na;
extern crate num_complex;
#[allow(unused_variables)]
use na::{Vector4};
// use num_complex::Complex;
use na::linalg::FourierTransform;

fn main() {
    let f = FourierTransform::new(::convert(Vector4::new(-1.0,1.5,2.0,3.0)));

    println!("{}",f.fft());
}
