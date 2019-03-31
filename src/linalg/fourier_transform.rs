use alga::general::{ClosedAdd, ClosedMul, Real, SubsetOf};
use base::allocator::Allocator;
use base::storage::Storage;
use base::{DefaultAllocator, Dim, DimName, MatrixN, VectorN};
use num_complex::Complex;
use base::dimension::{U1,U2,DimDiv,DimSub,DimMin,DimQuot};
use std::ops::Rem;

pub trait DimRem{
    
}

impl Rem for Dim {
    pub fn rem(self)
}


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
    pub v: VectorN<Complex<N>, D>,
    ///
    pub r: VectorN<N, D>,
    is_real: bool,
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
    /// Initializes with a complex vector
    pub fn new(vector: VectorN<Complex<N>, D>) -> Self {
        Self {
            r: VectorN::zeros_generic(vector.data.shape().0, vector.data.shape().1),
            v: vector,
            is_real: false,
        }
    }
    /// Initializes with a real vector
    pub fn from_real(real: VectorN<N, D>) -> Self {
        Self {
            v: VectorN::zeros_generic(real.data.shape().0, real.data.shape().1),
            r: real,
            is_real: true,
        }
    }

    ///
    /// Returns the real component of Fourier Transformed vector
    ///
    #[inline]
    pub fn dft(&self) -> VectorN<Complex<N>, D> {
        let (nrows, ncols) = self.v.data.shape();
        let mut vec = VectorN::zeros_generic(nrows, ncols);

        for i in 0..nrows.value() {
            for j in 0..nrows.value() {
                let arg =
                    N::frac_2_pi() * ::convert(i) * ::convert(j) / ::convert(nrows.value());
                vec[i] += Complex::new(N::cos(arg), -N::sin(arg));
            }
            if self.is_real {
                vec[i] = vec[i] * self.r[i];
            } else {
                vec[i] *= self.v[i];
            }
        }
        vec
    }

    ///
    /// Returns the real component of Fourier Transformed vector
    ///
    #[inline]
    pub fn fft<C1,C2>(&self) -> VectorN<Complex<N>, D> 
    where
        C1: DimSub<D>,
        D: DimDiv<U2, Output=DimQuot<D,U2>>,
        C2: DimMin<DimQuot<D,U2>>        
    {
        let (nrows, ncols) = self.v.data.shape();
        let mut vec = VectorN::zeros_generic(nrows, ncols);
        if nrows.value() == U1.value() {
            if self.is_real {
                vec[0] = ::convert(self.r[0]);
            } else {
                vec[0] = self.v[0];
            }
        } else {
            let c1 = self.v.data.shape().0.div(U2);
            println!("{:?}",c1);
            // let _n = VectorN::zeros_generic(C1)

        }

        vec
    }
}
