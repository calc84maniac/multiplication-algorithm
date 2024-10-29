#include "bitwise.h"
#include <stdbool.h>
#include <stdint.h>

#ifndef CARRY_ONLY
#define CARRY_ONLY 1
#endif

#ifndef CARRY_DIRECT
#define CARRY_DIRECT 1
#endif

#if CARRY_DIRECT
// Takes a multiplier between -0x01000000 and 0x00FFFFFF
static inline bool booths_multiplication32_opt(u32 multiplicand, u32 multiplier, u32 accumulator) {
    // Set the low bit of the multiplicand to cause negation to invert the upper bits, this bit can't propagate to bit 31
    multiplicand |= 1;

    // Optimized first iteration
    u32 booth = (s32)(multiplier << 31) >> 31;
    u32 carry = booth * multiplicand;
    // Pre-populate accumulator for output
    u32 output = accumulator;

    // Determine number of iterations and bit to start at
    u32 csa_mask = UINT32_C(0x80000000);
    u32 temp = multiplier ^ ((s32)multiplier >> 31);
    do {
        csa_mask >>= 4;
        temp >>= 8;
    } while (temp != 0);

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
#else
// Takes a multiplier between -0x01000000 and 0x00FFFFFF
static inline bool booths_multiplication32_opt(u32 multiplicand, u32 multiplier, u32 accumulator) {
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
    bool carry;
    if (rs < -0x01000000 || rs >= 0x01000000) {
        carry = (rs >> 30) == -2;
    } else {
        carry = booths_multiplication32_opt(rm, rs, rn);
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
    bool carry;
    if (rs >= 0x01000000) {
        carry = booths_multiplication64_opt(rm >> 6, rs >> 26, rdhi);
    } else {
        carry = booths_multiplication32_opt(rm, rs, rdlo);
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
    bool carry;
    if (rs < -0x01000000 || rs >= 0x01000000) {
        carry = booths_multiplication64_opt(rm >> 6, rs >> 26, rdhi);
    } else {
        carry = booths_multiplication32_opt(rm, rs, rdlo);
    }
    return (struct MultiplicationOutput) { product, carry };
#else
    return booths_multiplication_opt(LONG_SIGNED, rm, rs, (u64) rdhi << 32 | rdlo);
#endif
}

struct MultiplicationOutput smull_opt(s32 rm, s32 rs) {
    return smlal_opt(0, 0, rm, rs);
}

