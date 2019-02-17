extern crate nalgebra as na;
extern crate num_complex;
extern crate alga;
use na::{convert, DMatrix, Real, Scalar};
use num_complex::Complex;



fn fourier_transform_matrix<N: Real>(size: usize) -> DMatrix<Complex<N>> {
    let mut mat = DMatrix::<Complex<N>>::zeros(size, size);

    for i in 0..size {
        for j in 0..i {
            let arg =  N::pi() * convert(2.0f64) * convert(i as f64) * convert(j as f64) / convert(size as f64);
            mat[(i, j)] = Complex::new(N::cos(arg), -N::sin(arg));
        }
    }
    mat = (mat + mat.transpose()) / convert(2.0f64);

    return mat;
}

// fn dft<N1: Real, D: Dim, S: Storage<N,D>, N2: Complex, S2 : Storage<N2,D>>(vector: Vector<N1,D,S>) ->  Vector<N2,D,S>{
//     let N = vector.len();
//     for i in 0..N {

//     }
// }

fn main() {
    let M = fourier_transform_matrix::<f32>(3);
    println!("{:?}", M);
}
