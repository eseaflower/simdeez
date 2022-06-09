use super::*;

pub struct Avx512;
impl Simd for Avx512 {
    type Vi16 = I16x32;
    type Vi32 = I32x16;
    type Vf32 = F32x16;
    type Vf64 = F64x8;
    type Vi64 = I64x8;

    const VF32_WIDTH: usize = 16;
    const VF64_WIDTH: usize = 8;
    const VI16_WIDTH: usize = 32;
    const VI32_WIDTH: usize = 16;
    const VI64_WIDTH: usize = 8;

    #[inline(always)]
    unsafe fn abs_ps(a: Self::Vf32) -> Self::Vf32 {
        let b = _mm512_set1_ps(-0.0f32);
        F32x16(_mm512_castsi512_ps(_mm512_andnot_si512(
            _mm512_castps_si512(b),
            _mm512_castps_si512(a.0),
        )))
        // F32x16(_mm512_andnot_ps(b, a.0))
    }
    #[inline(always)]
    unsafe fn abs_pd(a: Self::Vf64) -> Self::Vf64 {
        let b = _mm512_set1_pd(-0.0f64);
        F64x8(_mm512_castsi512_pd(_mm512_andnot_si512(
            _mm512_castpd_si512(b),
            _mm512_castpd_si512(a.0),
        )))
        // F64x8(_mm512_andnot_pd(b, a.0))
    }
    #[inline(always)]
    unsafe fn mullo_epi16(a: Self::Vi16, b: Self::Vi16) -> Self::Vi16 {
        I16x32(_mm512_mullo_epi16(a.0, b.0))
    }
    #[inline(always)]
    unsafe fn andnot_ps(a: Self::Vf32, b: Self::Vf32) -> Self::Vf32 {
        F32x16(_mm512_castsi512_ps(_mm512_andnot_si512(
            _mm512_castps_si512(a.0),
            _mm512_castps_si512(b.0),
        )))
        // F32x16(_mm512_andnot_ps(a.0, b.0))
    }
    #[inline(always)]
    unsafe fn andnot_pd(a: Self::Vf64, b: Self::Vf64) -> Self::Vf64 {
        F64x8(_mm512_castsi512_pd(_mm512_andnot_si512(
            _mm512_castpd_si512(a.0),
            _mm512_castpd_si512(b.0),
        )))
        // F64x8(_mm512_andnot_pd(a.0, b.0))
    }
    #[inline(always)]
    unsafe fn andnot_epi32(a: Self::Vi32, b: Self::Vi32) -> Self::Vi32 {
        I32x16(_mm512_andnot_si512(a.0, b.0))
    }
    #[inline(always)]
    unsafe fn andnot_epi64(a: Self::Vi64, b: Self::Vi64) -> Self::Vi64 {
        I64x8(_mm512_andnot_si512(a.0, b.0))
    }
    #[inline(always)]
    unsafe fn blendv_epi32(a: Self::Vi32, b: Self::Vi32, mask: Self::Vi32) -> Self::Vi32 {
        let kmask = _mm512_cmpneq_epi32_mask(mask.0, _mm512_setzero_epi32());
        I32x16(_mm512_mask_blend_epi32(kmask, a.0, b.0))
        // I32x16(_mm512_castps_si512(_mm512_blendv_ps(
        //     _mm512_castsi512_ps(a.0),
        //     _mm512_castsi512_ps(b.0),
        //     _mm512_castsi512_ps(mask.0),
        // )))
    }
    #[inline(always)]
    unsafe fn blendv_epi64(a: Self::Vi64, b: Self::Vi64, mask: Self::Vi64) -> Self::Vi64 {
        let kmask = _mm512_cmpneq_epi64_mask(mask.0, _mm512_setzero_epi32());
        I64x8(_mm512_mask_blend_epi64(kmask, a.0, b.0))
        // I64x8(_mm512_castpd_si512(_mm512_blendv_pd(
        //     _mm512_castsi512_pd(a.0),
        //     _mm512_castsi512_pd(b.0),
        //     _mm512_castsi512_pd(mask.0),
        // )))
    }
    #[inline(always)]
    unsafe fn blendv_ps(a: Self::Vf32, b: Self::Vf32, mask: Self::Vf32) -> Self::Vf32 {
        let kmask = _mm512_cmpneq_ps_mask(mask.0, _mm512_setzero_ps());
        F32x16(_mm512_mask_blend_ps(kmask, a.0, b.0))
        // F32x16(_mm512_blendv_ps(a.0, b.0, mask.0))
    }
    #[inline(always)]
    unsafe fn blendv_pd(a: Self::Vf64, b: Self::Vf64, mask: Self::Vf64) -> Self::Vf64 {
        let kmask = _mm512_cmpneq_pd_mask(mask.0, _mm512_setzero_pd());
        F64x8(_mm512_mask_blend_pd(kmask, a.0, b.0))
        // F64x8(_mm512_blendv_pd(a.0, b.0, mask.0))
    }
    #[inline(always)]
    unsafe fn castps_epi32(a: Self::Vf32) -> Self::Vi32 {
        I32x16(_mm512_castps_si512(a.0))
    }
    #[inline(always)]
    unsafe fn castpd_epi64(a: Self::Vf64) -> Self::Vi64 {
        I64x8(_mm512_castpd_si512(a.0))
    }
    #[inline(always)]
    unsafe fn castepi32_ps(a: Self::Vi32) -> Self::Vf32 {
        F32x16(_mm512_castsi512_ps(a.0))
    }
    #[inline(always)]
    unsafe fn castepi64_pd(a: Self::Vi64) -> Self::Vf64 {
        F64x8(_mm512_castsi512_pd(a.0))
    }
    #[inline(always)]
    unsafe fn castps_pd(a: Self::Vf32) -> Self::Vf64 {
        F64x8(_mm512_castps_pd(a.0))
    }
    #[inline(always)]
    unsafe fn castpd_ps(a: Self::Vf64) -> Self::Vf32 {
        F32x16(_mm512_castpd_ps(a.0))
    }
    #[inline(always)]
    unsafe fn cmpeq_epi32(a: Self::Vi32, b: Self::Vi32) -> Self::Vi32 {
        let kmask = _mm512_cmpeq_epi32_mask(a.0, b.0);
        I32x16(_mm512_maskz_mov_epi32(kmask, _mm512_set1_epi32(-1)))
        // I32x16(_mm512_cmpeq_epi32(a.0, b.0))
    }
    #[inline(always)]
    unsafe fn cmpneq_epi32(a: Self::Vi32, b: Self::Vi32) -> Self::Vi32 {
        let kmask = _mm512_cmpneq_epi32_mask(a.0, b.0);
        I32x16(_mm512_maskz_mov_epi32(kmask, _mm512_set1_epi32(-1)))
        // Self::not_epi32(I32x16(_mm512_cmpeq_epi32(a.0, b.0)))
    }
    #[inline(always)]
    unsafe fn cmpge_epi32(a: Self::Vi32, b: Self::Vi32) -> Self::Vi32 {
        let kmask = _mm512_cmpge_epi32_mask(a.0, b.0);
        I32x16(_mm512_maskz_mov_epi32(kmask, _mm512_set1_epi32(-1)))
        // Self::not_epi32(I32x16(_mm512_cmpgt_epi32(b.0, a.0)))
    }
    #[inline(always)]
    unsafe fn cmpgt_epi32(a: Self::Vi32, b: Self::Vi32) -> Self::Vi32 {
        let kmask = _mm512_cmpgt_epi32_mask(a.0, b.0);
        I32x16(_mm512_maskz_mov_epi32(kmask, _mm512_set1_epi32(-1)))
        // I32x16(_mm512_cmpgt_epi32(a.0, b.0))
    }
    #[inline(always)]
    unsafe fn cmple_epi32(a: Self::Vi32, b: Self::Vi32) -> Self::Vi32 {
        let kmask = _mm512_cmple_epi32_mask(a.0, b.0);
        I32x16(_mm512_maskz_mov_epi32(kmask, _mm512_set1_epi32(-1)))
        // Self::not_epi32(I32x16(_mm512_cmpgt_epi32(a.0, b.0)))
    }
    #[inline(always)]
    unsafe fn cmplt_epi32(a: Self::Vi32, b: Self::Vi32) -> Self::Vi32 {
        let kmask = _mm512_cmplt_epi32_mask(a.0, b.0);
        I32x16(_mm512_maskz_mov_epi32(kmask, _mm512_set1_epi32(-1)))
        // I32x16(_mm512_cmpgt_epi32(b.0, a.0))
    }
    #[inline(always)]
    unsafe fn cmpeq_epi64(a: Self::Vi64, b: Self::Vi64) -> Self::Vi64 {
        let kmask = _mm512_cmpeq_epi64_mask(a.0, b.0);
        I64x8(_mm512_maskz_mov_epi64(kmask, _mm512_set1_epi64(-1)))
        // I64x8(_mm512_cmpeq_epi64(a.0, b.0))
    }
    #[inline(always)]
    unsafe fn cmpneq_epi64(a: Self::Vi64, b: Self::Vi64) -> Self::Vi64 {
        let kmask = _mm512_cmpneq_epi64_mask(a.0, b.0);
        I64x8(_mm512_maskz_mov_epi64(kmask, _mm512_set1_epi64(-1)))
        // Self::not_epi64(I64x8(_mm512_cmpeq_epi64(a.0, b.0)))
    }
    #[inline(always)]
    unsafe fn cmpge_epi64(a: Self::Vi64, b: Self::Vi64) -> Self::Vi64 {
        let kmask = _mm512_cmpge_epi64_mask(a.0, b.0);
        I64x8(_mm512_maskz_mov_epi64(kmask, _mm512_set1_epi64(-1)))
        // Self::not_epi64(I64x8(_mm512_cmpgt_epi64(b.0, a.0)))
    }
    #[inline(always)]
    unsafe fn cmpgt_epi64(a: Self::Vi64, b: Self::Vi64) -> Self::Vi64 {
        let kmask = _mm512_cmpgt_epi64_mask(a.0, b.0);
        I64x8(_mm512_maskz_mov_epi64(kmask, _mm512_set1_epi64(-1)))
        // I64x8(_mm512_cmpgt_epi64(a.0, b.0))
    }
    #[inline(always)]
    unsafe fn cmple_epi64(a: Self::Vi64, b: Self::Vi64) -> Self::Vi64 {
        let kmask = _mm512_cmple_epi64_mask(a.0, b.0);
        I64x8(_mm512_maskz_mov_epi64(kmask, _mm512_set1_epi64(-1)))
        // Self::not_epi64(I64x8(_mm512_cmpgt_epi64(a.0, b.0)))
    }
    #[inline(always)]
    unsafe fn cmplt_epi64(a: Self::Vi64, b: Self::Vi64) -> Self::Vi64 {
        let kmask = _mm512_cmplt_epi64_mask(a.0, b.0);
        I64x8(_mm512_maskz_mov_epi64(kmask, _mm512_set1_epi64(-1)))
        // I64x8(_mm512_cmpgt_epi64(b.0, a.0))
    }
    #[inline(always)]
    unsafe fn cmpeq_ps(a: Self::Vf32, b: Self::Vf32) -> Self::Vf32 {
        let kmask = _mm512_cmpeq_ps_mask(a.0, b.0);
        F32x16(_mm512_maskz_mov_ps(
            kmask,
            _mm512_castsi512_ps(_mm512_set1_epi32(-1)),
        ))
        // F32x16(_mm512_cmp_ps(a.0, b.0, _CMP_EQ_OQ))
    }
    #[inline(always)]
    unsafe fn cmpneq_ps(a: Self::Vf32, b: Self::Vf32) -> Self::Vf32 {
        let kmask = _mm512_cmpneq_ps_mask(a.0, b.0);
        F32x16(_mm512_maskz_mov_ps(
            kmask,
            _mm512_castsi512_ps(_mm512_set1_epi32(-1)),
        ))
        // F32x16(_mm512_cmp_ps(a.0, b.0, _CMP_NEQ_OQ))
    }
    #[inline(always)]
    unsafe fn cmpge_ps(a: Self::Vf32, b: Self::Vf32) -> Self::Vf32 {
        let kmask = _mm512_cmp_ps_mask::<_CMP_GE_OQ>(a.0, b.0);
        F32x16(_mm512_maskz_mov_ps(
            kmask,
            _mm512_castsi512_ps(_mm512_set1_epi32(-1)),
        ))
        // F32x16(_mm512_cmp_ps(a.0, b.0, _CMP_GE_OQ))
    }
    #[inline(always)]
    unsafe fn cmpgt_ps(a: Self::Vf32, b: Self::Vf32) -> Self::Vf32 {
        let kmask = _mm512_cmp_ps_mask::<_CMP_GT_OQ>(a.0, b.0);
        F32x16(_mm512_maskz_mov_ps(
            kmask,
            _mm512_castsi512_ps(_mm512_set1_epi32(-1)),
        ))
        // F32x16(_mm512_cmp_ps(a.0, b.0, _CMP_GT_OQ))
    }
    #[inline(always)]
    unsafe fn cmple_ps(a: Self::Vf32, b: Self::Vf32) -> Self::Vf32 {
        let kmask = _mm512_cmp_ps_mask::<_CMP_LE_OQ>(a.0, b.0);
        F32x16(_mm512_maskz_mov_ps(
            kmask,
            _mm512_castsi512_ps(_mm512_set1_epi32(-1)),
        ))
        // F32x16(_mm512_cmp_ps(a.0, b.0, _CMP_LE_OQ))
    }
    #[inline(always)]
    unsafe fn cmplt_ps(a: Self::Vf32, b: Self::Vf32) -> Self::Vf32 {
        let kmask = _mm512_cmp_ps_mask::<_CMP_LT_OQ>(a.0, b.0);
        F32x16(_mm512_maskz_mov_ps(
            kmask,
            _mm512_castsi512_ps(_mm512_set1_epi32(-1)),
        ))
        // F32x16(_mm512_cmp_ps(a.0, b.0, _CMP_LT_OQ))
    }
    #[inline(always)]
    unsafe fn cmpeq_pd(a: Self::Vf64, b: Self::Vf64) -> Self::Vf64 {
        let kmask = _mm512_cmpeq_pd_mask(a.0, b.0);
        F64x8(_mm512_maskz_mov_pd(
            kmask,
            _mm512_castsi512_pd(_mm512_set1_epi64(-1)),
        ))
        // F64x8(_mm512_cmp_pd(a.0, b.0, _CMP_EQ_OQ))
    }
    #[inline(always)]
    unsafe fn cmpneq_pd(a: Self::Vf64, b: Self::Vf64) -> Self::Vf64 {
        let kmask = _mm512_cmpneq_pd_mask(a.0, b.0);
        F64x8(_mm512_maskz_mov_pd(
            kmask,
            _mm512_castsi512_pd(_mm512_set1_epi64(-1)),
        ))
        // F64x8(_mm512_cmp_pd(a.0, b.0, _CMP_NEQ_OQ))
    }
    #[inline(always)]
    unsafe fn cmpge_pd(a: Self::Vf64, b: Self::Vf64) -> Self::Vf64 {
        let kmask = _mm512_cmp_pd_mask::<_CMP_GE_OQ>(a.0, b.0);
        F64x8(_mm512_maskz_mov_pd(
            kmask,
            _mm512_castsi512_pd(_mm512_set1_epi64(-1)),
        ))
        // F64x8(_mm512_cmp_pd(a.0, b.0, _CMP_GE_OQ))
    }
    #[inline(always)]
    unsafe fn cmpgt_pd(a: Self::Vf64, b: Self::Vf64) -> Self::Vf64 {
        let kmask = _mm512_cmp_pd_mask::<_CMP_GT_OQ>(a.0, b.0);
        F64x8(_mm512_maskz_mov_pd(
            kmask,
            _mm512_castsi512_pd(_mm512_set1_epi64(-1)),
        ))
        // F64x8(_mm512_cmp_pd(a.0, b.0, _CMP_GT_OQ))
    }
    #[inline(always)]
    unsafe fn cmple_pd(a: Self::Vf64, b: Self::Vf64) -> Self::Vf64 {
        let kmask = _mm512_cmp_pd_mask::<_CMP_LE_OQ>(a.0, b.0);
        F64x8(_mm512_maskz_mov_pd(
            kmask,
            _mm512_castsi512_pd(_mm512_set1_epi64(-1)),
        ))
        // F64x8(_mm512_cmp_pd(a.0, b.0, _CMP_LE_OQ))
    }
    #[inline(always)]
    unsafe fn cmplt_pd(a: Self::Vf64, b: Self::Vf64) -> Self::Vf64 {
        let kmask = _mm512_cmp_pd_mask::<_CMP_LT_OQ>(a.0, b.0);
        F64x8(_mm512_maskz_mov_pd(
            kmask,
            _mm512_castsi512_pd(_mm512_set1_epi64(-1)),
        ))
        // F64x8(_mm512_cmp_pd(a.0, b.0, _CMP_LT_OQ))
    }
    #[inline(always)]
    unsafe fn cvtepi32_ps(a: Self::Vi32) -> Self::Vf32 {
        F32x16(_mm512_cvtepi32_ps(a.0))
    }
    #[inline(always)]
    unsafe fn cvtepi64_pd(a: Self::Vi64) -> Self::Vf64 {
        let x = _mm512_add_epi64(
            a.0,
            _mm512_castpd_si512(_mm512_set1_pd(core::mem::transmute::<i64, f64>(
                0x0018000000000000,
            ))),
        );
        F64x8(_mm512_sub_pd(
            _mm512_castsi512_pd(x),
            _mm512_set1_pd(core::mem::transmute::<i64, f64>(0x0018000000000000)),
        ))
    }
    #[inline(always)]
    unsafe fn cvtps_epi32(a: Self::Vf32) -> Self::Vi32 {
        I32x16(_mm512_cvtps_epi32(a.0))
    }
    #[inline(always)]
    unsafe fn cvtpd_epi64(a: Self::Vf64) -> Self::Vi64 {
        let x = _mm512_add_pd(
            a.0,
            _mm512_set1_pd(core::mem::transmute::<i64, f64>(0x0018000000000000)),
        );
        I64x8(_mm512_sub_epi64(
            _mm512_castpd_si512(x),
            _mm512_castpd_si512(_mm512_set1_pd(core::mem::transmute::<i64, f64>(
                0x0018000000000000,
            ))),
        ))
    }
    #[inline(always)]
    unsafe fn ceil_ps(a: Self::Vf32) -> Self::Vf32 {
        F32x16(_mm512_roundscale_ps::<_MM_FROUND_CEIL>(a.0))
        // F32x16(_mm512_ceil_ps(a.0))
    }
    #[inline(always)]
    unsafe fn fast_ceil_ps(a: Self::Vf32) -> Self::Vf32 {
        F32x16(_mm512_roundscale_ps::<_MM_FROUND_CEIL>(a.0))
        // F32x16(_mm512_ceil_ps(a.0))
    }
    #[inline(always)]
    unsafe fn ceil_pd(a: Self::Vf64) -> Self::Vf64 {
        F64x8(_mm512_roundscale_pd::<_MM_FROUND_CEIL>(a.0))
        // F64x8(_mm512_ceil_pd(a.0))
    }
    #[inline(always)]
    unsafe fn floor_ps(a: Self::Vf32) -> Self::Vf32 {
        F32x16(_mm512_roundscale_ps::<_MM_FROUND_FLOOR>(a.0))
        // F32x16(_mm512_floor_ps(a.0))
    }
    #[inline(always)]
    unsafe fn floor_pd(a: Self::Vf64) -> Self::Vf64 {
        F64x8(_mm512_roundscale_pd::<_MM_FROUND_FLOOR>(a.0))
        // F64x8(_mm512_floor_pd(a.0))
    }
    #[inline(always)]
    unsafe fn fast_floor_ps(a: Self::Vf32) -> Self::Vf32 {
        F32x16(_mm512_roundscale_ps::<_MM_FROUND_FLOOR>(a.0))
        // F32x16(_mm512_floor_ps(a.0))
    }
    #[inline(always)]
    unsafe fn fast_floor_pd(a: Self::Vf64) -> Self::Vf64 {
        F64x8(_mm512_roundscale_pd::<_MM_FROUND_FLOOR>(a.0))
        // F64x8(_mm512_floor_pd(a.0))
    }
    #[inline(always)]
    unsafe fn fmadd_ps(a: Self::Vf32, b: Self::Vf32, c: Self::Vf32) -> Self::Vf32 {
        F32x16(_mm512_fmadd_ps(a.0, b.0, c.0))
    }
    #[inline(always)]
    unsafe fn fnmadd_ps(a: Self::Vf32, b: Self::Vf32, c: Self::Vf32) -> Self::Vf32 {
        F32x16(_mm512_fnmadd_ps(a.0, b.0, c.0))
    }
    #[inline(always)]
    unsafe fn fmadd_pd(a: Self::Vf64, b: Self::Vf64, c: Self::Vf64) -> Self::Vf64 {
        F64x8(_mm512_fmadd_pd(a.0, b.0, c.0))
    }
    #[inline(always)]
    unsafe fn fnmadd_pd(a: Self::Vf64, b: Self::Vf64, c: Self::Vf64) -> Self::Vf64 {
        F64x8(_mm512_fnmadd_pd(a.0, b.0, c.0))
    }
    #[inline(always)]
    unsafe fn fmsub_ps(a: Self::Vf32, b: Self::Vf32, c: Self::Vf32) -> Self::Vf32 {
        F32x16(_mm512_fmsub_ps(a.0, b.0, c.0))
    }
    #[inline(always)]
    unsafe fn fnmsub_ps(a: Self::Vf32, b: Self::Vf32, c: Self::Vf32) -> Self::Vf32 {
        F32x16(_mm512_fnmsub_ps(a.0, b.0, c.0))
    }
    #[inline(always)]
    unsafe fn fmsub_pd(a: Self::Vf64, b: Self::Vf64, c: Self::Vf64) -> Self::Vf64 {
        F64x8(_mm512_fmsub_pd(a.0, b.0, c.0))
    }
    #[inline(always)]
    unsafe fn fnmsub_pd(a: Self::Vf64, b: Self::Vf64, c: Self::Vf64) -> Self::Vf64 {
        F64x8(_mm512_fnmsub_pd(a.0, b.0, c.0))
    }
    #[inline(always)]
    unsafe fn horizontal_add_ps(a: Self::Vf32) -> f32 {
        _mm512_reduce_add_ps(a.0)
        // let mut vlow = _mm512_castps512_ps128(a.0);
        // let vhigh = _mm512_extractf128_ps(a.0, 1);
        // vlow = _mm_add_ps(vlow, vhigh);
        // let mut shuf = _mm_movehdup_ps(vlow);
        // let mut sums = _mm_add_ps(vlow, shuf);
        // shuf = _mm_movehl_ps(shuf, sums);
        // sums = _mm_add_ss(sums, shuf);
        // _mm_cvtss_f32(sums)
    }
    #[inline(always)]
    unsafe fn horizontal_add_pd(a: Self::Vf64) -> f64 {
        _mm512_reduce_add_pd(a.0)
        // let mut vlow = _mm512_castpd512_pd128(a.0);
        // let vhigh = _mm512_extractf128_pd(a.0, 1);
        // vlow = _mm_add_pd(vlow, vhigh);
        // let high64 = _mm_unpackhi_pd(vlow, vlow);
        // _mm_cvtsd_f64(_mm_add_sd(vlow, high64))
    }
    #[inline(always)]
    unsafe fn i32gather_epi32(arr: &[i32], index: Self::Vi32) -> Self::Vi32 {
        I32x16(_mm512_i32gather_epi32(
            index.0,
            arr.as_ptr() as *const u8,
            4,
        ))
        // I32x16(_mm512_i32gather_epi32(&arr[0] as *const i32, index.0, 4))
    }
    #[inline(always)]
    unsafe fn i64gather_epi64(arr: &[i64], index: Self::Vi64) -> Self::Vi64 {
        I64x8(_mm512_i64gather_epi64(
            index.0,
            arr.as_ptr() as *const u8,
            4,
        ))
        // I64x8(_mm512_i64gather_epi64(&arr[0] as *const i64, index.0, 8))
    }
    #[inline(always)]
    unsafe fn i32gather_ps(arr: &[f32], index: Self::Vi32) -> Self::Vf32 {
        F32x16(_mm512_i32gather_ps(index.0, arr.as_ptr() as *const u8, 4))
        // F32x16(_mm512_i32gather_ps(&arr[0] as *const f32, index.0, 4))
    }
    #[inline(always)]
    unsafe fn load_ps(a: &f32) -> Self::Vf32 {
        F32x16(_mm512_load_ps(a as *const f32))
    }
    #[inline(always)]
    unsafe fn load_pd(a: &f64) -> Self::Vf64 {
        F64x8(_mm512_load_pd(a as *const f64))
    }
    #[inline(always)]
    unsafe fn load_epi16(a: &i16) -> Self::Vi16 {
        I16x32(_mm512_load_si512(a as *const _ as *const i32))
        // let m = mem::transmute::<&i16, &__m512i>(a);
        // I16x32(_mm512_load_si512(m))
    }
    #[inline(always)]
    unsafe fn load_epi32(a: &i32) -> Self::Vi32 {
        I32x16(_mm512_load_si512(a))
        // let m = mem::transmute::<&i32, &__m512i>(a);
        // I32x16(_mm512_load_si512(m))
    }
    #[inline(always)]
    unsafe fn load_epi64(a: &i64) -> Self::Vi64 {
        I64x8(_mm512_load_si512(a as *const _ as *const i32))
        // let m = mem::transmute::<&i64, &__m512i>(a);
        // I64x8(_mm512_load_si512(m))
    }
    #[inline(always)]
    unsafe fn loadu_ps(a: &f32) -> Self::Vf32 {
        F32x16(_mm512_loadu_ps(a as *const f32))
    }
    #[inline(always)]
    unsafe fn loadu_pd(a: &f64) -> Self::Vf64 {
        F64x8(_mm512_loadu_pd(a as *const f64))
    }
    #[inline(always)]
    unsafe fn loadu_epi32(a: &i32) -> Self::Vi32 {
        I32x16(_mm512_loadu_si512(a))
        // let m = mem::transmute::<&i32, &__m512i>(a);
        // I32x16(_mm512_loadu_si512(m))
    }
    #[inline(always)]
    unsafe fn loadu_epi64(a: &i64) -> Self::Vi64 {
        I64x8(_mm512_loadu_si512(a as *const _ as *const i32))
        // let m = mem::transmute::<&i64, &__m512i>(a);
        // I64x8(_mm512_loadu_si512(m))
    }
    #[inline(always)]
    unsafe fn maskload_epi32(mem_addr: &i32, mask: Self::Vi32) -> Self::Vi32 {
        let kmask = _mm512_cmpneq_epi32_mask(mask.0, _mm512_setzero_epi32());
        I32x16(_mm512_maskz_load_epi32(kmask, mem_addr))
        // I32x16(_mm512_maskload_epi32(mem_addr as *const i32, mask.0))
    }
    #[inline(always)]
    unsafe fn maskload_epi64(mem_addr: &i64, mask: Self::Vi64) -> Self::Vi64 {
        let kmask = _mm512_cmpneq_epi64_mask(mask.0, _mm512_setzero_si512());
        I64x8(_mm512_maskz_load_epi64(kmask, mem_addr))
        // I64x8(_mm512_maskload_epi64(mem_addr as *const i64, mask.0))
    }
    #[inline(always)]
    unsafe fn maskload_ps(mem_addr: &f32, mask: Self::Vi32) -> Self::Vf32 {
        let kmask = _mm512_cmpneq_epi32_mask(mask.0, _mm512_setzero_epi32());
        F32x16(_mm512_maskz_load_ps(kmask, mem_addr))
        // F32x16(_mm512_maskload_ps(mem_addr as *const f32, mask.0))
    }
    #[inline(always)]
    unsafe fn maskload_pd(mem_addr: &f64, mask: Self::Vi64) -> Self::Vf64 {
        let kmask = _mm512_cmpneq_epi64_mask(mask.0, _mm512_setzero_si512());
        F64x8(_mm512_maskz_load_pd(kmask, mem_addr))
        // F64x8(_mm512_maskload_pd(mem_addr as *const f64, mask.0))
    }
    #[inline(always)]
    unsafe fn store_ps(mem_addr: &mut f32, a: Self::Vf32) {
        _mm512_store_ps(mem_addr as *mut f32, a.0);
    }
    #[inline(always)]
    unsafe fn store_pd(mem_addr: &mut f64, a: Self::Vf64) {
        _mm512_store_pd(mem_addr as *mut f64, a.0);
    }
    #[inline(always)]
    unsafe fn store_epi32(mem_addr: &mut i32, a: Self::Vi32) {
        _mm512_store_si512(mem_addr, a.0);
        // let mem_addr_512 = mem::transmute::<&mut i32, &mut __m512i>(mem_addr);
        // _mm512_store_si512(mem_addr_512, a.0);
    }
    #[inline(always)]
    unsafe fn store_epi64(mem_addr: &mut i64, a: Self::Vi64) {
        _mm512_store_si512(mem_addr as *mut _ as *mut i32, a.0);
        // let mem_addr_512 = mem::transmute::<&mut i64, &mut __m512i>(mem_addr);
        // _mm512_store_si512(mem_addr_512, a.0);
    }
    #[inline(always)]
    unsafe fn storeu_ps(mem_addr: &mut f32, a: Self::Vf32) {
        _mm512_storeu_ps(mem_addr as *mut f32, a.0);
    }
    #[inline(always)]
    unsafe fn storeu_pd(mem_addr: &mut f64, a: Self::Vf64) {
        _mm512_storeu_pd(mem_addr as *mut f64, a.0);
    }
    #[inline(always)]
    unsafe fn storeu_epi32(mem_addr: &mut i32, a: Self::Vi32) {
        _mm512_storeu_si512(mem_addr, a.0);
        // let mem_addr_512 = mem::transmute::<&mut i32, &mut __m512i>(mem_addr);
        // _mm512_storeu_si512(mem_addr_512, a.0);
    }
    #[inline(always)]
    unsafe fn storeu_epi64(mem_addr: &mut i64, a: Self::Vi64) {
        _mm512_storeu_si512(mem_addr as *mut _ as *mut i32, a.0);
        // let mem_addr_512 = mem::transmute::<&mut i64, &mut __m512i>(mem_addr);
        // _mm512_storeu_si512(mem_addr_512, a.0);
    }
    #[inline(always)]
    unsafe fn maskstore_epi32(mem_addr: &mut i32, mask: Self::Vi32, a: Self::Vi32) {
        let kmask = _mm512_cmpneq_epi32_mask(mask.0, _mm512_setzero_epi32());
        _mm512_mask_store_epi32(mem_addr, kmask, a.0)
        // _mm512_maskstore_epi32(mem_addr as *mut i32, mask.0, a.0)
    }
    #[inline(always)]
    unsafe fn maskstore_epi64(mem_addr: &mut i64, mask: Self::Vi64, a: Self::Vi64) {
        let kmask = _mm512_cmpneq_epi64_mask(mask.0, _mm512_setzero_epi32());
        _mm512_mask_store_epi64(mem_addr, kmask, a.0)
        // _mm512_maskstore_epi64(mem_addr as *mut i64, mask.0, a.0)
    }
    #[inline(always)]
    unsafe fn maskstore_ps(mem_addr: &mut f32, mask: Self::Vi32, a: Self::Vf32) {
        let kmask = _mm512_cmpneq_epi32_mask(mask.0, _mm512_setzero_epi32());
        _mm512_mask_store_ps(mem_addr, kmask, a.0)
        // _mm512_maskstore_ps(mem_addr as *mut f32, mask.0, a.0)
    }
    #[inline(always)]
    unsafe fn maskstore_pd(mem_addr: &mut f64, mask: Self::Vi64, a: Self::Vf64) {
        let kmask = _mm512_cmpneq_epi64_mask(mask.0, _mm512_setzero_epi32());
        _mm512_mask_store_pd(mem_addr, kmask, a.0)
        // _mm512_maskstore_pd(mem_addr as *mut f64, mask.0, a.0)
    }
    #[inline(always)]
    unsafe fn max_epi32(a: Self::Vi32, b: Self::Vi32) -> Self::Vi32 {
        I32x16(_mm512_max_epi32(a.0, b.0))
    }
    #[inline(always)]
    unsafe fn min_epi32(a: Self::Vi32, b: Self::Vi32) -> Self::Vi32 {
        I32x16(_mm512_min_epi32(a.0, b.0))
    }
    #[inline(always)]
    unsafe fn max_ps(a: Self::Vf32, b: Self::Vf32) -> Self::Vf32 {
        F32x16(_mm512_max_ps(a.0, b.0))
    }
    #[inline(always)]
    unsafe fn min_ps(a: Self::Vf32, b: Self::Vf32) -> Self::Vf32 {
        F32x16(_mm512_min_ps(a.0, b.0))
    }
    #[inline(always)]
    unsafe fn max_pd(a: Self::Vf64, b: Self::Vf64) -> Self::Vf64 {
        F64x8(_mm512_max_pd(a.0, b.0))
    }
    #[inline(always)]
    unsafe fn min_pd(a: Self::Vf64, b: Self::Vf64) -> Self::Vf64 {
        F64x8(_mm512_min_pd(a.0, b.0))
    }
    #[inline(always)]
    unsafe fn mullo_epi32(a: Self::Vi32, b: Self::Vi32) -> Self::Vi32 {
        I32x16(_mm512_mullo_epi32(a.0, b.0))
    }
    #[inline(always)]
    unsafe fn mullo_epi64(a: Self::Vi64, b: Self::Vi64) -> Self::Vi64 {
        let mut result = Self::setzero_epi64();
        result[0] = a[0] * b[0];
        result[1] = a[1] * b[1];
        result[2] = a[2] * b[2];
        result[3] = a[3] * b[3];
        result[4] = a[4] * b[4];
        result[5] = a[5] * b[5];
        result[6] = a[6] * b[6];
        result[7] = a[7] * b[7];
        result
    }
    #[inline(always)]
    unsafe fn rcp_ps(a: Self::Vf32) -> Self::Vf32 {
        F32x16(_mm512_rcp14_ps(a.0))
        // F32x16(_mm512_rcp_ps(a.0))
    }
    #[inline(always)]
    unsafe fn round_ps(a: Self::Vf32) -> Self::Vf32 {
        F32x16(_mm512_roundscale_ps(
            a.0,
            _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC,
        ))
        // F32x16(_mm512_round_ps(
        //     a.0,
        //     _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC,
        // ))
    }
    #[inline(always)]
    unsafe fn fast_round_ps(a: Self::Vf32) -> Self::Vf32 {
        F32x16(_mm512_roundscale_ps(
            a.0,
            _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC,
        ))
        // F32x16(_mm512_round_ps(
        //     a.0,
        //     _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC,
        // ))
    }
    #[inline(always)]
    unsafe fn round_pd(a: Self::Vf64) -> Self::Vf64 {
        F64x8(_mm512_roundscale_pd(
            a.0,
            _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC,
        ))
        // F64x8(_mm512_round_pd(
        //     a.0,
        //     _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC,
        // ))
    }
    #[inline(always)]
    unsafe fn set1_epi32(a: i32) -> Self::Vi32 {
        I32x16(_mm512_set1_epi32(a))
    }
    #[inline(always)]
    unsafe fn set1_epi64(a: i64) -> Self::Vi64 {
        I64x8(_mm512_set1_epi64(a))
        // I64x8(_mm512_set1_epi64x(a))
    }
    #[inline(always)]
    unsafe fn set1_ps(a: f32) -> Self::Vf32 {
        F32x16(_mm512_set1_ps(a))
    }
    #[inline(always)]
    unsafe fn set1_pd(a: f64) -> Self::Vf64 {
        F64x8(_mm512_set1_pd(a))
    }
    #[inline(always)]
    unsafe fn setzero_pd() -> Self::Vf64 {
        F64x8(_mm512_setzero_pd())
    }
    #[inline(always)]
    unsafe fn setzero_ps() -> Self::Vf32 {
        F32x16(_mm512_setzero_ps())
    }
    #[inline(always)]
    unsafe fn setzero_epi32() -> Self::Vi32 {
        I32x16(_mm512_setzero_si512())
    }
    #[inline(always)]
    unsafe fn setzero_epi64() -> Self::Vi64 {
        I64x8(_mm512_setzero_si512())
    }
    #[inline(always)]
    unsafe fn srai_epi64(a: Self::Vi64, amt_const: i32) -> Self::Vi64 {
        macro_rules! call {
            ($amt_const:expr) => {
                I64x8(_mm512_srai_epi64(a.0, $amt_const))
            };
        }
        constify_imm8!(amt_const, call)
        // instruction does not exist. Split into 32-bit shifts
        // if amt_const <= 32 {
        //     let bb = _mm_set_epi32(0, 0, 0, amt_const);
        //     let sra = _mm512_sra_epi32(a.0, bb); // a >> b signed dwords
        //     let srl = _mm512_srl_epi64(a.0, bb); // a >> b unsigned qwords
        //     let mask = _mm512_setr_epi32(0, -1, 0, -1, 0, -1, 0, -1); // mask for signed high part
        //     Self::blendv_epi64(I64x8(srl), I64x8(sra), I64x8(mask))
        // } else {
        //     // b > 32
        //     let bm32 = _mm_set_epi32(0, 0, 0, amt_const - 32);
        //     let sign = _mm512_srai_epi32(a.0, 31); // sign of a
        //     let sra2 = _mm512_sra_epi32(a.0, bm32); // a >> (b-32) signed dwords
        //     let sra3 = _mm512_srli_epi64(sra2, 32); // a >> (b-32) >> 32 (second shift unsigned qword)
        //     let mask = _mm512_setr_epi32(0, -1, 0, -1, 0, -1, 0, -1); // mask for high part containing only sign
        //     Self::blendv_epi64(I64x8(sra3), I64x8(sign), I64x8(mask))
        // }
    }
    #[inline(always)]
    unsafe fn srli_epi32(a: Self::Vi32, amt_const: i32) -> Self::Vi32 {
        macro_rules! call {
            ($amt_const:expr) => {
                I32x16(_mm512_srli_epi32(a.0, $amt_const))
            };
        }
        constify_imm8!(amt_const, call)
    }
    #[inline(always)]
    unsafe fn sra_epi32(a: Self::Vi32, amt: i32) -> Self::Vi32 {
        I32x16(_mm512_sra_epi32(a.0, _mm_set1_epi32(amt)))
    }
    #[inline(always)]
    unsafe fn srl_epi32(a: Self::Vi32, amt: i32) -> Self::Vi32 {
        I32x16(_mm512_srl_epi32(a.0, _mm_set1_epi32(amt)))
    }
    #[inline(always)]
    unsafe fn sll_epi32(a: Self::Vi32, amt: i32) -> Self::Vi32 {
        I32x16(_mm512_sll_epi32(a.0, _mm_set1_epi32(amt)))
    }

    #[inline(always)]
    unsafe fn sqrt_ps(a: Self::Vf32) -> Self::Vf32 {
        F32x16(_mm512_sqrt_ps(a.0))
    }
    #[inline(always)]
    unsafe fn rsqrt_ps(a: Self::Vf32) -> Self::Vf32 {
        F32x16(_mm512_rsqrt14_ps(a.0))
        // F32x16(_mm512_rsqrt_ps(a.0))
    }
    #[inline(always)]
    unsafe fn sqrt_pd(a: Self::Vf64) -> Self::Vf64 {
        F64x8(_mm512_sqrt_pd(a.0))
    }
    #[inline(always)]
    unsafe fn rsqrt_pd(a: Self::Vf64) -> Self::Vf64 {
        F64x8(_mm512_div_pd(_mm512_set1_pd(1.0), _mm512_sqrt_pd(a.0)))
    }
    #[inline(always)]
    unsafe fn shuffle_epi32(a: Self::Vi32, imm8: i32) -> I32x16 {
        macro_rules! call {
            ($imm8:expr) => {
                I32x16(_mm512_shuffle_epi32(a.0, $imm8))
            };
        }
        constify_imm8!(imm8, call)
    }
}
