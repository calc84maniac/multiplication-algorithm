#include "bitwise.h"
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#ifndef CARRY_ONLY
#define CARRY_ONLY 1
#endif

#ifndef CARRY_DIRECT
#define CARRY_DIRECT 1
#endif

#ifndef USE_SSE2
#define USE_SSE2 __SSE2__
#endif

#ifndef USE_AVX2
#define USE_AVX2 0
#endif

#ifndef USE_SWAR64
#define USE_SWAR64 0
#endif

#ifndef USE_SWAR32
#define USE_SWAR32 0
#endif

#if defined(__ARM_ARCH_4T__)
// Takes a multiplier between -0x01000000 and 0x00FFFFFF, cycles between 0 and 2
bool booths_multiplication32_opt(u32 multiplicand, u32 multiplier, u32 accumulator, u32 cycles) {
    (void)cycles;
    // Set the low bit of the multiplicand to cause negation to invert the upper bits, this bit can't propagate to bit 31
    multiplicand |= 1;

    // Optimized first iteration
    u32 booth = (s32)(multiplier << 31) >> 31;
    u32 carry = booth & ~multiplicand;
    // Pre-populate accumulator for output
    u32 output = accumulator;

    u32 sum = output + carry;
    do {
        #pragma GCC unroll 4
        for (int shift = 0; shift < 8; shift += 2) {
            // Get next booth factor (-2 to 2, shifted left by shift)
            u32 next_booth = (s32)(multiplier << (29 - shift)) >> (30 - shift);
            u32 factor = next_booth - booth;
            booth = next_booth;
            // Get scaled value of booth addend
            u32 addend = multiplicand * factor;
            // Combine the addend with the CSA
            // Not performing any masking seems to work because the lower carries can't propagate to bit 31
            output ^= carry ^ (addend << 1);
            sum += addend << 1;
            carry = sum - output;
        }
        // Adjust next iteration to use booth factors between -128 and 128
        booth = (s32)booth >> 8;
        multiplier = (s32)multiplier >> 8;
        multiplicand <<= 8;
    } while (multiplier != booth);

    return carry >> 31;
}
#else
#if CARRY_DIRECT
#if USE_SWAR64
// Takes a multiplier between -0x01000000 and 0x00FFFFFF, cycles between 0 and 2
static inline bool booths_multiplication32_opt(u32 multiplicand, u32 multiplier, u32 accumulator, u32 cycles) {
    // Determine masks to apply to each iteration
    // Each mask is shifted left by 1 to account for the extra factor of 2
    u32 csa_shift = cycles * 4;
    u64 csa_masks = UINT64_C(0x0800040002000100) >> csa_shift;
    u32 csa_mask = (UINT32_C(1) << 27) >> csa_shift;

    // Optimized first iteration
    u32 carry = -(multiplier & 1) & ~multiplicand;
    // Pre-populate accumulator for output
    u32 output = accumulator;

    // Add the bits relevant to carry bit 31
    output = (output & -csa_mask) + (carry & -csa_mask);
    carry &= csa_mask;

    const u64 field_magic  = UINT64_C(0x0001000100010001);
    const u64 booth_magic  = UINT64_C(0x00C00030000C0003);
    const u64 buffer_magic = UINT64_C(0x0010000400010000);
    // Rotate the multiplicand left by 16+2 bits
    multiplicand = (multiplicand << 18) | (multiplicand >> 14);
    u64 outputs = 0;
    u64 carries = 0;
    do {
        // Repeat multiplier bits across a 64-bit value
        u64 booths = (u16)multiplier * field_magic;
        // Calculate 4 scaled booth factors (divided by 2), normalized within each 16 bits
        u64 factors = (booths & booth_magic) - ((booths >> 1) & booth_magic);
        // Multiply 16 bits of the multiplicand by the factors,
        // and simulate a 1 bit below the multiplicand (to invert for negative factors) by adding half the factors
        u64 addends = (s16)multiplicand * factors + ((s64)factors >> 1);
        // Add buffer bits to block positive or negative carries between fields
        // This works because the multiplicand bits are treated as signed,
        // so each multiplicand can at most contribute (-0x8000*2)+1 == (0x7FFF*-2)-1 or (0x7FFF*2)+1 == (-0x8000*-2)-1
        addends += buffer_magic;
        // Mask and accumulate the outputs and carries
        csa_masks <<= 4;
        carries += addends & csa_masks;
        outputs += addends & (field_magic + ~csa_masks);
        // Handle the next bytes of the multiplier and the multiplicand
        multiplier >>= 8;
        multiplicand = (multiplicand << 8) | (multiplicand >> 24);
    } while ((s64)csa_masks >= 0);

    // Reduce the outputs and carries and combine with the initial values, removing the factor of 2
    output += (outputs * field_magic) >> 33;
    carry += (carries * field_magic) >> 33 & UINT32_C(0xFFFF0000);
    // Detect carry bit 31 of the final iteration
    carry = output ^ (output - carry);
    return carry >> 31;
}
#elif USE_SWAR32
// For now, this is a simple translation of 64-bit SWAR with carries between 32-bit halves eliminated
// Takes a multiplier between -0x01000000 and 0x00FFFFFF, cycles between 0 and 2
static inline bool booths_multiplication32_opt(u32 multiplicand, u32 multiplier, u32 accumulator, u32 cycles) {
    // Determine masks to apply to each iteration
    // Each mask is shifted left by 1 to account for the extra factor of 2
    u32 csa_shift = cycles * 4;
    u32 csa_masks = UINT32_C(0x08000400) >> csa_shift;
    u32 csa_mask = (UINT32_C(1) << 27) >> csa_shift;

    // Optimized first iteration
    u32 carry = -(multiplier & 1) & ~multiplicand;
    // Pre-populate accumulator for output
    u32 output = accumulator;

    // Add the bits relevant to carry bit 31
    output = (output & -csa_mask) + (carry & -csa_mask);
    carry &= csa_mask;

    const u32 field_magic   = UINT32_C(0x00010001);
    const u32 booth_magic1  = UINT32_C(0x000C0003);
    const u32 booth_magic2  = UINT32_C(0x00C00030);
    const u32 buffer_magic1 = UINT32_C(0x00010000);
    const u32 buffer_magic2 = UINT32_C(0x00100000);
    // Rotate the multiplicand left by 16+2 bits
    multiplicand = (multiplicand << 18) | (multiplicand >> 14);
    u32 outputs = 0;
    u32 carries = 0;
    do {
        // Repeat multiplier bits across a 64-bit value
        u32 booths = (u16)multiplier * field_magic;
        // Calculate 4 scaled booth factors (divided by 2), normalized within each 16 bits
        u32 factors1 = (booths & booth_magic1) - ((booths >> 1) & booth_magic1);
        u32 factors2 = (booths & booth_magic2) - ((booths >> 1) & booth_magic2);
        // Multiply 16 bits of the multiplicand by the factors,
        // and simulate a 1 bit below the multiplicand (to invert for negative factors) by adding half the factors
        u32 addends1 = (s16)multiplicand * factors1 + ((s32)factors1 >> 1);
        u32 addends2 = (s16)multiplicand * factors2 + ((s32)factors2 >> 1);
        // Add buffer bits to block positive or negative carries between fields
        // This works because the multiplicand bits are treated as signed,
        // so each multiplicand can at most contribute (-0x8000*2)+1 == (0x7FFF*-2)-1 or (0x7FFF*2)+1 == (-0x8000*-2)-1
        addends1 += buffer_magic1;
        addends2 += buffer_magic2;
        // Mask and accumulate the outputs and carries
        csa_masks <<= 2;
        carries += addends1 & csa_masks;
        outputs += addends1 & (field_magic + ~csa_masks);
        csa_masks <<= 2;
        carries += addends2 & csa_masks;
        outputs += addends2 & (field_magic + ~csa_masks);
        // Handle the next bytes of the multiplier and the multiplicand
        multiplier >>= 8;
        multiplicand = (multiplicand << 8) | (multiplicand >> 24);
    } while ((s32)csa_masks >= 0);

    // Reduce the outputs and carries and combine with the initial values, removing the factor of 2
    output += (outputs * field_magic) >> 1;
    carry += (carries * field_magic) >> 1 & UINT32_C(0xFFFF0000);
    // Detect carry bit 31 of the final iteration
    carry = output ^ (output - carry);
    return carry >> 31;
}
#elif USE_AVX2
#include <immintrin.h>
// Takes a multiplier between -0x01000000 and 0x00FFFFFF, cycles between 0 and 2
static inline bool booths_multiplication32_opt(u32 multiplicand, u32 multiplier, u32 accumulator, u32 cycles) {
    // Calculate all scaled booth factors (divided by 2), normalized within 16 bits
    const __m256i factor_shuffle = _mm256_set_epi8(
        -1,-1,-1,-1,-1,-1,-1,-1,
         3, 2, 3, 2, 3, 2, 3, 2,
         2, 1, 2, 1, 2, 1, 2, 1,
         1, 0, 1, 0, 1, 0, 1, 0);
    const __m256i factor_mask = _mm256_set1_epi64x(UINT64_C(0x00C00030000C0003));

    __m256i factors = _mm256_set1_epi32(multiplier);
    factors = _mm256_shuffle_epi8(factors, factor_shuffle);
    __m256i factors_low = _mm256_and_si256(factors, factor_mask);
    __m256i factors_high = _mm256_and_si256(_mm256_srli_epi16(factors, 1), factor_mask);
    factors = _mm256_sub_epi16(factors_low, factors_high);
    __m256i factor_signs = _mm256_cmpgt_epi16(factors_high, factors_low);

    // Get the multiplicand (times 4) shifted relative to each booth factor
    const __m256i multiplicand_shuffle = _mm256_set_epi8(
        -1,-1,-1,-1,-1,-1,-1,-1,
         1, 0, 1, 0, 1, 0, 1, 0,
         2, 1, 2, 1, 2, 1, 2, 1,
         3, 2, 3, 2, 3, 2, 3, 2);
 
    __m256i multiplicands = _mm256_set1_epi32(multiplicand << 2);
    multiplicands = _mm256_shuffle_epi8(multiplicands, multiplicand_shuffle);

    // Calculate the 16 most significant bits of every booth addend (times 2)
    __m256i addends = _mm256_mullo_epi16(factors, multiplicands);
    // Subtract 1 if the factor was negative to remove the injected carries
    addends = _mm256_add_epi16(addends, factor_signs);

    // Determine masks to apply to each iteration
    // Each mask is shifted left by 1 to account for the extra factor of 2
    __m256i csa_masks = _mm256_set_epi16(
        (1 << 15), (1 << 14), (1 << 13), (1 << 12),
        (1 << 11), (1 << 10), (1 << 9), (1 << 8),
        (1 << 7), (1 << 6), (1 << 5), (1 << 4),
        (1 << 3), (1 << 2), (1 << 1), (1 << 0));
    u32 csa_shift = (3 - cycles) * 4;
    csa_masks = _mm256_sll_epi16(csa_masks, _mm_cvtsi32_si128(csa_shift));
    u32 csa_mask = (UINT32_C(1) << 15) << csa_shift;

    // Optimized first iteration
    u32 carry = -(multiplier & 1) & ~multiplicand;
    // Pre-populate accumulator for output
    u32 output = accumulator;

    // Add the bits relevant to carry bit 31
    output = (output & -csa_mask) + (carry & -csa_mask);
    carry &= csa_mask;

    // Mask and reduce the outputs and carries
    __m256i outputs = _mm256_and_si256(addends, _mm256_sub_epi16(_mm256_setzero_si256(), csa_masks));
    __m256i carries = _mm256_and_si256(addends, csa_masks);
    // Pack output into the high 64 bits and carry into the low 64 bits of each 128-bit lane
    __m256i reduce = _mm256_hadd_epi16(carries, outputs);
    // Combine lanes
    __m128i reduce2 = _mm_add_epi16(_mm256_castsi256_si128(reduce), _mm256_extracti128_si256(reduce, 1));
    // Pack output into the high 32 bits and carry into the low 32 bits
    reduce2 = _mm_add_epi16(_mm_shuffle_epi32(reduce2, _MM_SHUFFLE(3, 1, 2, 0)),
                            _mm_shuffle_epi32(reduce2, _MM_SHUFFLE(2, 0, 3, 1)));
    // Pack output into the high 16 bits and carry into the low 16 bits
    reduce2 = _mm_add_epi16(_mm_shufflelo_epi16(reduce2, _MM_SHUFFLE(3, 1, 2, 0)),
                            _mm_shufflelo_epi16(reduce2, _MM_SHUFFLE(2, 0, 3, 1)));
    // Get packed outputs and carries
    u32 reduce3 = _mm_cvtsi128_si32(reduce2);
    // Combine with the initial scalar outputs and carries, removing the factor of 2
    // No need to mask the low bits of the output, which don't affect carry detection
    output += reduce3 >> 1;
    carry += reduce3 << 15;
    // Detect carry bit 31 of the final iteration
    carry = output ^ (output - carry);
    return carry >> 31;
}
#elif USE_SSE2
#include <immintrin.h>
// Help the compiler interleave scalar and vector operations
static inline void booth_iter(u32 multiplicand, u32 multiplier, u32* sum, u32* output, u32* carry, u32* booth, int shift)
{
    // Get next booth factor (-2 to 2, shifted left by 30-shift)
    u32 next_booth = (s32)(multiplier << shift) >> shift;
    u32 factor = next_booth - *booth;
    *booth = next_booth;
    // Get scaled value of booth addend
    u32 addend = multiplicand * factor;
    // Combine the addend with the CSA
    // Not performing any masking seems to work because the lower carries can't propagate to bit 31
    *output ^= *carry ^ addend;
    *sum += addend;
    *carry = *sum - *output;
}

// Takes a multiplier between -0x01000000 and 0x00FFFFFF, cycles between 0 and 2
static inline bool booths_multiplication32_opt(u32 multiplicand, u32 multiplier, u32 accumulator, u32 cycles) {
    // Set the low bit of the multiplicand to cause negation to invert the upper bits, this bit can't propagate to bit 31
    multiplicand |= 1;

    // Optimized initial iteration
    u32 booth = (s32)(multiplier << 31) >> 31;
    u32 carry = booth & ~multiplicand;
    // Pre-populate accumulator for output
    u32 output = accumulator;
    u32 sum = output + carry;

    // Calculate the last 8 scaled booth factors, normalized within 16 bits
    const __m128i factor_mask = _mm_set_epi16(
        (3 << 7), (3 << 7), (3 << 5), (3 << 5),
        (3 << 3), (3 << 3), (3 << 1), (3 << 1));

    u32 packed_multiplier = (multiplier & UINT32_C(0xFFFF0000)) | ((multiplier >> 8) & UINT32_C(0x0000FFFF));
    __m128i factors = _mm_set1_epi32(packed_multiplier);
    __m128i factors_low = _mm_and_si128(_mm_slli_epi16(factors, 1), factor_mask);
    __m128i factors_high = _mm_and_si128(factors, factor_mask);
    factors = _mm_sub_epi16(factors_low, factors_high);

    // First iteration
    booth_iter(multiplicand, multiplier, &sum, &output, &carry, &booth, 29);

    // Get the multiplicand shifted relative to each booth factor, times 2
    u32 packed_multiplicand = (multiplicand << 17) | ((multiplicand >> 7 | 1) & UINT32_C(0x0000FFFF));
    __m128i multiplicands = _mm_set1_epi32(packed_multiplicand);

    // Second iteration
    booth_iter(multiplicand, multiplier, &sum, &output, &carry, &booth, 27);

    // Calculate bits 23-30 of the last 8 booth addends
    __m128i addends = _mm_mullo_epi16(factors, multiplicands);
    addends = _mm_srli_epi16(addends, 8);

    // Third iteration
    booth_iter(multiplicand, multiplier, &sum, &output, &carry, &booth, 25);

    // Determine masks to apply to each 8-bit addend slice
    __m128i csa_masks = _mm_set_epi16(
        -(1 << 7), -(1 << 3), -(1 << 6), -(1 << 2), -(1 << 5), -(1 << 1), -(1 << 4), -(1 << 0));
    u32 csa_shift = (2 - cycles) * 4;
    csa_masks = _mm_sll_epi16(csa_masks, _mm_cvtsi32_si128(csa_shift));
    u32 csa_mask = -(UINT32_C(1) << 23) << csa_shift;

    // Fourth iteration
    booth_iter(multiplicand, multiplier, &sum, &output, &carry, &booth, 23);

    // Add the bits relevant to carry bit 31
    sum = output + (carry & csa_mask);
    output += (carry & (csa_mask << 1));

    // Mask the sums and the outputs
    __m128i sums = _mm_and_si128(addends, csa_masks);
    __m128i outputs = _mm_and_si128(addends, _mm_slli_epi16(csa_masks, 1));

    // Pack the sums into the top 64 bits and the outputs into the bottom 64 bits
    __m128i reduce = _mm_packus_epi16(outputs, sums);
    // Reduce each group into a single sum
    reduce = _mm_sad_epu8(reduce, _mm_setzero_si128());
    // Combine with the initial scalar sum and output, shifting back up to bits 23-30
    sum += _mm_extract_epi16(reduce, 4) << 23;
    output += _mm_cvtsi128_si32(reduce) << 23;

    // Detect carry bit 31 of the final iteration
    carry = sum ^ output;
    return carry >> 31;
}
#else
// Takes a multiplier between -0x01000000 and 0x00FFFFFF, cycles between 0 and 2
static inline bool booths_multiplication32_opt(u32 multiplicand, u32 multiplier, u32 accumulator, u32 cycles) {
    // Set the low bit of the multiplicand to cause negation to invert the upper bits, this bit can't propagate to bit 31
    multiplicand |= 1;

    // Optimized first iteration
    u32 booth = (s32)(multiplier << 31) >> 31;
    u32 carry = booth * multiplicand;
    // Pre-populate accumulator for output
    u32 output = accumulator;

    // Determine bit to start at
    u32 csa_shift = cycles * 4;
    u32 csa_mask = UINT32_C(0x08000000) >> csa_shift;

    // Add the bits relevant to carry bit 31
    output = (output & -csa_mask) + (carry & -csa_mask);
    carry &= csa_mask;
    int shift = 29;
    do {
        // Let the compiler unroll
        for (int i = 0; i < 4; i++) {
            // Get next booth factor (-2 to 2, shifted left by 30-shift)
            u32 next_booth = (s32)(multiplier << shift) >> shift;
            u32 factor = next_booth - booth;
            booth = next_booth;
            // Get scaled value of booth addend
            u32 addend = multiplicand * factor;
            // Aggregate bits from each iteration that propagate to carry bit 31
            output += addend & -csa_mask;
            carry += addend & csa_mask;
            csa_mask <<= 1;
            shift -= 2;
        }
    } while ((s32)csa_mask >= 0);
    // Detect carry bit 31 of the final iteration
    carry = output ^ (output - carry);
    return carry >> 31;
}
#endif
#else
// Takes a multiplier between -0x01000000 and 0x00FFFFFF, cycles between 0 and 2
static inline bool booths_multiplication32_opt(u32 multiplicand, u32 multiplier, u32 accumulator, u32 cycles) {
    (void)cycles;
    // Set the low bit of the multiplicand to cause negation to invert the upper bits, this bit can't propagate to bit 31
    multiplicand |= 1;

    // Optimized first iteration
    u32 booth = (s32)(multiplier << 31) >> 31;
    u32 carry = booth * multiplicand;
    // Pre-populate accumulator for output
    u32 output = accumulator;

    u32 sum = output + carry;
    int shift = 29;
    do {
        for (int i = 0; i < 4; i++, shift -= 2) {
            // Get next booth factor (-2 to 2, shifted left by 30-shift)
            u32 next_booth = (s32)(multiplier << shift) >> shift;
            u32 factor = next_booth - booth;
            booth = next_booth;
            // Get scaled value of booth addend
            u32 addend = multiplicand * factor;
            // Combine the addend with the CSA
            // Not performing any masking seems to work because the lower carries can't propagate to bit 31
            output ^= carry ^ addend;
            sum += addend;
            carry = sum - output;
        }
    } while (booth != multiplier);

    return carry >> 31;
}
#endif
#endif

// Takes a multiplicand shifted right by 6 and a multiplier shifted right by 26 (zero or sign extended)
static inline bool booths_multiplication64_opt(u32 multiplicand, u32 multiplier, u32 accum_hi) {
    // Skipping the first 14 iterations seems to work because the lower carries can't propagate to bit 63
    // This means only magic bits 62-61 are needed (which requires decoding 3 booth chunks),
    // and only the last two booth iterations are needed

    // Set the low bit of the multiplicand to cause negation to invert the upper bits
    multiplicand |= 1;

    // Pre-populate magic bit 61 for carry
    u32 carry = ~accum_hi & UINT32_C(0x20000000);
    // Pre-populate magic bits 63-60 for output (with carry magic pre-added in)
    u32 output = accum_hi - UINT32_C(0x08000000);

    // Get factors from the top 3 booth chunks
    u32 booth0 = (s32)(multiplier << 27) >> 27;
    u32 booth1 = (s32)(multiplier << 29) >> 29;
    u32 booth2 = (s32)(multiplier << 31) >> 31;
    u32 factor0 = multiplier - booth0;
    u32 factor1 = booth0 - booth1;
    u32 factor2 = booth1 - booth2;

    // Get scaled value of the 3rd top booth addend
    u32 addend = multiplicand * factor2;
    // Finalize bits 61-60 of output magic using its sign
    output -= addend & UINT32_C(0x10000000);
    // Get scaled value of the 2nd top booth addend
    addend = multiplicand * factor1;
    // Finalize bits 63-62 of output magic using its sign
    output -= addend & UINT32_C(0x40000000);

    // Get the carry from the CSA in bit 61 and propagate it to bit 62, which is not processed in this iteration
    u32 sum = output + (addend & UINT32_C(0x20000000));
    // Subtract out the carry magic to get the actual output magic
    output -= carry;

    // Get scaled value of the 1st top booth addend
    addend = multiplicand * factor0;
    // Add to bit 62 and propagate the carry
    sum += addend & UINT32_C(0x40000000);

    // Cancel out the output magic bit 63 to get the carry bit 63
    return (sum ^ output) >> 31;
}

static inline struct MultiplicationOutput booths_multiplication_opt(enum MultiplicationFlavor flavor, s64 multiplicand, s64 multiplier, u64 accumulator) {
    bool alu_carry_in = multiplier & 1;

    // Optimized first iteration
    u64 carry = (multiplier & 1) ? (~multiplicand & UINT64_C(0x3FFFFFFFF)) : 0;
    // Pre-populate magic bits for carry
    carry |= ~accumulator & UINT64_C(0xAAAAAAA800000000);
    // Pre-populate magic bits and lower accumulator for output
    u64 output = accumulator & UINT64_C(0x55555557FFFFFFFF);
    // Add inverted upper bits of carry into output magic bits
    output += (~carry & UINT64_C(0xAAAAAAAA00000000)) << 1;

    // Extract signs from each booth chunk
    u32 booth_signs = (multiplier >> 2) & UINT32_C(0x55555555);
    // Take the absolute value of each chunk
    u32 booth_chunks = multiplier ^ (booth_signs * 3);
    // Determine which chunks are non-zero
    u32 booth_nonzero = booth_chunks | booth_chunks >> 1;
    // Finalize output magic using the sign of each booth addend
    output += (~(booth_nonzero & (booth_signs ^ (multiplicand >> 63))) & UINT32_C(0x55555555)) << 34;

    // Determine which booth chunks inject a carry
    booth_signs &= booth_nonzero;
    // Translate chunk absolute values into direct factors of 0, 1, 2
    booth_chunks -= (booth_chunks >> 1) & UINT32_C(0x55555555);

    // Start at bit 1 to account for the optimized first iteration
    u64 csa_mask = UINT64_C(0x1FFFFFFFF) << 1;
    multiplicand <<= 1;
    u32 booth_mask = 3;
    multiplier ^= multiplier >> 63;
    u64 sum = output + carry;
    do {
        for (int i = 0; i < 4; i++, csa_mask <<= 2, booth_mask <<= 2) {
            // Get absolute value of booth addend, pre-scaled
            u64 addend = (u64)multiplicand * (booth_chunks & booth_mask);
            // Invert the addend if there's an injected carry
            addend ^= -(u64)(booth_signs & booth_mask);
            // Combine the addend with the CSA, within the current 33-bit mask
            addend &= csa_mask;
            output ^= (carry & csa_mask) ^ addend;
            sum += addend;
            carry = sum - output;
        }
    } while (multiplier >>= 8);

    // Inject booth carries
    carry += booth_signs << 1;

    u64 product = output + carry + alu_carry_in;
    if (is_long(flavor)) {
        return (struct MultiplicationOutput) {
            product,
            bit(carry, booth_mask ? 31 : 63)
        };
    } else {
        return (struct MultiplicationOutput) { 
            (u32) product,
            bit(carry, 31)
        };
    }
}

struct MultiplicationOutput mla_opt(s32 rm, s32 rs, u32 rn) {
#if CARRY_ONLY
    u32 product = (u32) rm * (u32) rs + rn;
    u32 cycles = (u32) (31 - __builtin_clz((rs ^ (rs >> 31)) | 1)) / 8;
    bool carry;
    if (cycles >= 3) {
        carry = (rs >> 30) == -2;
    } else {
        carry = booths_multiplication32_opt(rm, rs, rn, cycles);
    }
    return (struct MultiplicationOutput) { product, carry };
#else
    return booths_multiplication_opt(SHORT, rm, rs, rn);
#endif
}

struct MultiplicationOutput mul_opt(s32 rm, s32 rs) {
    return mla_opt(rm, rs, 0);
}

struct MultiplicationOutput umlal_opt(u32 rdlo, u32 rdhi, u32 rm, u32 rs) {
#if CARRY_ONLY
    u64 product = (u64) rm * (u64) rs + ((u64) rdhi << 32 | rdlo);
    u32 cycles = (u32) (31 - __builtin_clz(rs | 1)) / 8;
    bool carry;
    if (cycles >= 3) {
        carry = booths_multiplication64_opt(rm >> 6, rs >> 26, rdhi);
    } else {
        carry = booths_multiplication32_opt(rm, rs, rdlo, cycles);
    }
    return (struct MultiplicationOutput) { product, carry };
#else
    return booths_multiplication_opt(LONG_UNSIGNED, rm, rs, ((u64) rdhi << 32) | rdlo);
#endif
}

struct MultiplicationOutput umull_opt(u32 rm, u32 rs) {
    return umlal_opt(0, 0, rm, rs);
}

struct MultiplicationOutput smlal_opt(u32 rdlo, u32 rdhi, s32 rm, s32 rs) {
#if CARRY_ONLY
    u64 product = (u64) rm * (u64) rs + ((u64) rdhi << 32 | rdlo);
    u32 cycles = (u32) (31 - __builtin_clz((rs ^ (rs >> 31)) | 1)) / 8;
    bool carry;
    if (cycles >= 3) {
        carry = booths_multiplication64_opt(rm >> 6, rs >> 26, rdhi);
    } else {
        carry = booths_multiplication32_opt(rm, rs, rdlo, cycles);
    }
    return (struct MultiplicationOutput) { product, carry };
#else
    return booths_multiplication_opt(LONG_SIGNED, rm, rs, (u64) rdhi << 32 | rdlo);
#endif
}

struct MultiplicationOutput smull_opt(s32 rm, s32 rs) {
    return smlal_opt(0, 0, rm, rs);
}

