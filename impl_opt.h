#include "bitwise.h"
#include <stdbool.h>
#include <stdint.h>

static inline struct CSAOutput perform_csa_opt(struct CSAOutput prev_output, u64 addend, u64 mask) {
    addend &= mask;
    u64 output = prev_output.output ^ (prev_output.carry & mask) ^ addend;
    u64 carry  = (prev_output.output + prev_output.carry + addend) - output;
    return (struct CSAOutput) { output, carry };
}

static inline struct MultiplicationOutput booths_multiplication_opt(enum MultiplicationFlavor flavor, s64 multiplicand, s64 multiplier, u64 accumulator) {
    struct CSAOutput csa_output;
    bool alu_carry_in = multiplier & 1;

    // Optimized first iteration
    csa_output.carry = (multiplier & 1) ? (~multiplicand & UINT64_C(0x3FFFFFFFF)) : 0;
    // Pre-populate magic bits for carry
    csa_output.carry |= ~accumulator & UINT64_C(0xAAAAAAA800000000);
    // Pre-populate magic bits and lower accumulator for output
    csa_output.output = accumulator & UINT64_C(0x55555557FFFFFFFF);
    // Add inverted upper bits of carry into output magic bits
    csa_output.output += (~csa_output.carry & UINT64_C(0xAAAAAAAA00000000)) << 1;

    // Extract signs from each booth chunk
    u32 booth_signs = (multiplier >> 2) & UINT32_C(0x55555555);
    // Take the absolute value of each chunk
    u32 booth_chunks = multiplier ^ (booth_signs * 3);
    // Determine which chunks are non-zero
    u32 booth_nonzero = booth_chunks | booth_chunks >> 1;
    // Finalize output magic using the sign of each booth addend
    csa_output.output += (~(booth_nonzero & (booth_signs ^ (multiplicand >> 63))) & UINT32_C(0x55555555)) << 34;

    // Determine which booth chunks inject a carry
    booth_signs &= booth_nonzero;
    // Translate chunk absolute values into direct factors of 0, 1, 2
    booth_chunks -= (booth_chunks >> 1) & UINT32_C(0x55555555);

    // Start at bit 1 to account for the optimized first iteration
    u64 csa_mask = UINT64_C(0x1FFFFFFFF) << 1;
    multiplicand <<= 1;
    u32 booth_mask = 3;
    multiplier ^= multiplier >> 63;
    do {
        for (int i = 0; i < 4; i++, csa_mask <<= 2, booth_mask <<= 2) {
            // Get absolute value of booth addend, pre-scaled
            u64 addend = multiplicand * (booth_chunks & booth_mask);
            // Invert the addend if there's an injected carry
            addend ^= -(u64)(booth_signs & booth_mask);
            // Combine the addend with the CSA, within the current 33-bit mask
            csa_output = perform_csa_opt(csa_output, addend, csa_mask);
        }
    } while (multiplier >>= 8);

    // Inject booth carries
    csa_output.carry += booth_signs << 1;

    u64 product = csa_output.output + csa_output.carry + alu_carry_in;
    if (is_long(flavor)) {
        return (struct MultiplicationOutput) {
            product,
            bit(csa_output.carry, booth_mask ? 31 : 63)
        };
    } else {
        return (struct MultiplicationOutput) { 
            (u32) product,
            bit(csa_output.carry, 31)
        };
    }
}

struct MultiplicationOutput mul_opt(s32 rm, s32 rs) {
    return booths_multiplication_opt(SHORT, rm, rs, 0);
}

struct MultiplicationOutput mla_opt(s32 rm, s32 rs, u32 rn) {
    return booths_multiplication_opt(SHORT, rm, rs, rn);
}

struct MultiplicationOutput umull_opt(u32 rm, u32 rs) {
    return booths_multiplication_opt(LONG_UNSIGNED, rm, rs, 0);
}

struct MultiplicationOutput umlal_opt(u32 rdlo, u32 rdhi, u32 rm, u32 rs) {
    return booths_multiplication_opt(LONG_UNSIGNED, rm, rs, ((u64) rdhi << 32) | (u64) rdlo);
}

struct MultiplicationOutput smull_opt(s32 rm, s32 rs) {
    return booths_multiplication_opt(LONG_SIGNED, rm, rs, 0);
}

struct MultiplicationOutput smlal_opt(u32 rdlo, u32 rdhi, s32 rm, s32 rs) {
    return booths_multiplication_opt(LONG_SIGNED, rm, rs, (u64) rdhi << 32 | (u64) rdlo);
}