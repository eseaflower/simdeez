use super::*;
use core::arch::x86_64::*;
use core::mem;

mod avx512;
mod overloads;
pub use self::avx512::*;
pub use self::overloads::*;

#[derive(Copy, Debug, Clone)]
pub struct I16x32(pub __m512i);
impl SimdBase<I16x32, i16> for I16x32 {}
impl SimdSmallInt<I16x32, i16> for I16x32 {}

#[derive(Copy, Debug, Clone)]
pub struct I32x16(pub __m512i);
impl SimdBase<I32x16, i32> for I32x16 {}
impl SimdSmallInt<I32x16, i32> for I32x16 {}

#[derive(Copy, Debug, Clone)]
pub struct I64x8(pub __m512i);
impl SimdBase<I64x8, i64> for I64x8 {}

#[derive(Copy, Debug, Clone)]
pub struct F32x16(pub __m512);
impl SimdBase<F32x16, f32> for F32x16 {}
impl SimdFloat<F32x16, f32> for F32x16 {}

#[derive(Copy, Debug, Clone)]
pub struct F64x8(pub __m512d);
impl SimdBase<F64x8, f64> for F64x8 {}
impl SimdFloat<F64x8, f64> for F64x8 {}
