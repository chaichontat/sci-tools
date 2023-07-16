use std::simd::{u8x32, SimdPartialEq};

#[inline]
fn core_comp(sub_seq: &[u8; 32], sequence: &[u8; 32]) -> [i8; 32] {
    let simd_sub_seq = u8x32::from_slice(sub_seq);
    let simd_sequence = u8x32::from_slice(sequence);
    simd_sub_seq.simd_eq(simd_sequence).to_int().to_array()
}

// #[cfg(target_arch = "x86_64")]
// pub unsafe fn test(sub_seq: &[u8; 32], sequence: &[u8; 32]) -> [i8; 32] {
//     use std::arch::x86_64::*;
//     let mut mismatch_count = [0; 32];
//     _mm256_storeu_si256(
//         mismatch_count.as_mut_ptr() as *mut __m256i,
//         _mm256_cmpeq_epi8(
//             _mm256_loadu_si256(sub_seq.as_ptr() as *const __m256i),
//             _mm256_loadu_si256(sequence.as_ptr() as *const __m256i),
//         ),
//     );
//     mismatch_count
// }

// #[cfg(target_arch = "aarch64")]
// /// # Safety
// /// Loads data from raw pointers and assumes that they are valid and properly aligned.
// pub unsafe fn test(sub_seq: &[u8; 32], sequence: &[u8; 32]) -> [i8; 32] {
//     use std::arch::aarch64::*;
//     let mut mismatch_count = [0; 32];
//     vst1q_s8(
//         mismatch_count.as_mut_ptr() as *mut i8,
//         vreinterpretq_s8_u8(vceqq_u8(
//             vld1q_u8(sub_seq.as_ptr() as *const u8),
//             vld1q_u8(sequence.as_ptr() as *const u8),
//         )),
//     );
//     mismatch_count
// }

#[inline]
pub fn mismatch_count_dual(sub_seq: &[u8; 32], sequence: &[u8; 32], tol: i8) -> (bool, bool) {
    let mismatch_count = core_comp(sub_seq, sequence);
    (
        mismatch_count[0..16].iter().sum::<i8>() + 16 <= tol,
        mismatch_count[16..32].iter().sum::<i8>() + 16 <= tol,
    )
}

#[inline]
pub fn mismatch_count(sub_seq: &[u8; 32], sequence: &[u8; 32], tol: i8) -> bool {
    let mismatch_count = core_comp(sub_seq, sequence);
    mismatch_count.iter().sum::<i8>() + 32 <= tol
}
