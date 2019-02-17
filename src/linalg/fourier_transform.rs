use num_complex::Complex;
use alga::general::{ClosedAdd};
use std::ops::{Add, AddAssign};
use base::Scalar;

impl<T, Right> ClosedAdd<Right> for Complex<T>
where T : Real + Scalar + Add<Right, Output = T> + AddAssign<Right>
{
}

pub fn fourier_transform_matrix<N: Real>(size: usize) -> DMatrix<Complex<N>> {
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