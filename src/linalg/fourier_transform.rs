use alga::general::{ClosedAdd, ClosedMul, Real, SubsetOf};
use base::allocator::Allocator;
use base::storage::Storage;
use base::{DefaultAllocator, Dim, DimName, MatrixN, VectorN};
use num_complex::Complex;

///
#[derive(Clone, Debug)]
pub struct FourierTransform<N, D>
where
    N: Real,
    D: Dim,
    DefaultAllocator: Allocator<N, D, D>
        + Allocator<N, D>
        + Allocator<Complex<N>, D>
        + Allocator<Complex<N>, D, D>,
{
    ///
    pub v: VectorN<N, D>,
    ///
    pub z: VectorN<Complex<N>, D>,
    m: MatrixN<Complex<N>, D>,
}

impl<N, D> FourierTransform<N, D>
where
    N: Real,
    D: DimName,
    usize: SubsetOf<N>,
    N: SubsetOf<Complex<N>>,
    Complex<N>: ClosedMul + ClosedAdd,
    MatrixN<Complex<N>, D>: ClosedMul,
    DefaultAllocator: Allocator<N, D, D>
        + Allocator<N, D>
        + Allocator<Complex<N>, D>
        + Allocator<Complex<N>, D, D>,
{
    ///
    pub fn new(v: VectorN<N, D>) -> Self {
        let (nrows, _) = v.data.shape();
        let m = MatrixN::<Complex<N>, D>::from_fn(|i, j| {
            Complex::new(
                N::cos(N::two_pi() * ::convert(i) * ::convert(j) / ::convert(nrows.value())),
                -N::sin(N::two_pi() * ::convert(i) * ::convert(j) / ::convert(nrows.value())),
            )
        });
        let z: VectorN<Complex<N>, D> = ::convert(v.clone());
        Self { v: v, m: m, z: z }
    }

    ///
    pub fn as_real(&self) -> VectorN<N, D> {
        VectorN::<N, D>::from_fn(|i, _| self.z.dot(&self.m.row(i).transpose()).re)
    }

    ///
    pub fn as_complex(&self) -> VectorN<Complex<N>, D> {
        VectorN::<Complex<N>, D>::from_fn(|i, _| self.z.dot(&self.m.row(i).transpose()))
    }

    ///
    pub fn as_imag(&self) -> VectorN<N, D> {
        VectorN::<N, D>::from_fn(|i, _| self.z.dot(&self.m.row(i).transpose()).im)
    }
}
