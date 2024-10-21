#include "bitwise.h"
#include <stdbool.h>
#include <stdint.h>

// dont define PC_BUILD if you want to use this file in a GBA rom.
#ifdef PC_BUILD
#include <stdio.h>
#define printf(...) printf(__VA_ARGS__)
#else
#define printf(...)
#endif

// realistically this can only be a 3-bit value.
typedef u8 BoothChunk;

struct BoothRecodingOutput {
    u64  recoded_output;
    bool carry;
};

struct RecodedMultiplicands {
    struct BoothRecodingOutput m[4];
};

struct BoothRecodingOutput booth_recode(u64 input, BoothChunk booth_chunk) {
    struct BoothRecodingOutput output;
    switch (booth_chunk) {
        case 0: output = (struct BoothRecodingOutput) {            0, 0 }; break;
        case 1: output = (struct BoothRecodingOutput) {        input, 0 }; break;
        case 2: output = (struct BoothRecodingOutput) {        input, 0 }; break;
        case 3: output = (struct BoothRecodingOutput) {    2 * input, 0 }; break;
        case 4: output = (struct BoothRecodingOutput) { ~(2 * input), 1 }; break;
        case 5: output = (struct BoothRecodingOutput) {       ~input, 1 }; break;
        case 6: output = (struct BoothRecodingOutput) {       ~input, 1 }; break;
        case 7: output = (struct BoothRecodingOutput) {            0, 0 }; break;
    }

    output.recoded_output &= 0x3FFFFFFFFULL;
    return output;
}

struct CSAOutput {
    u64 output;
    u64 carry;
};

struct CSAOutput perform_csa(u64 a, u64 b, u64 c) {
    u64 output = a ^ b ^ c;
    u64 carry  = (a & b) | (b & c) | (c & a);
    return (struct CSAOutput) { output, carry };
}

// contains the current high 31 bits of the acc. this is shifted by 2 after each CSA.
u64 acc_shift_register = 0;

struct CSAOutput perform_csa_array(u64 partial_sum, u64 partial_carry, struct RecodedMultiplicands addends) {
    struct CSAOutput csa_output = { partial_sum, partial_carry };
    struct CSAOutput final_csa_output = { 0, 0 };

    for (int i = 0; i < 4; i++) {
        csa_output.output &= 0x1FFFFFFFFULL;
        csa_output.carry  &= 0x1FFFFFFFFULL;

        struct CSAOutput result = perform_csa(csa_output.output, addends.m[i].recoded_output & 0x1FFFFFFFFULL, csa_output.carry);

        // Inject the carry caused by booth recoding
        result.carry <<= 1;
        result.carry |= addends.m[i].carry;

        // Take the bottom two bits and inject them into the final output.
        // The value of the bottom two bits will not be changed by future
        // addends, because those addends must be at least 4 times as big
        // as the current addend. By directly injecting these two bits, the
        // hardware saves some space on the chip.
        final_csa_output.output |= (result.output & 3) << (2 * i);
        final_csa_output.carry  |= (result.carry  & 3) << (2 * i);
        
        // The next CSA will only operate on the upper bits - as explained
        // in the previous comment.
        result.output >>= 2;
        result.carry  >>= 2;

        // Perform the magic described in the tables for the handling of TransH
        // and High. acc_shift_register contains the upper 31 bits of the acc
        // in its lower bits.
        u64 magic = bit(acc_shift_register, 0) + !bit(csa_output.carry, 32) + !bit(addends.m[i].recoded_output, 33);
        result.output |= magic << 31;
        result.carry |= (u64) !bit(acc_shift_register, 1) << 32;        
        acc_shift_register >>= 2;

        csa_output = result;
    }

    final_csa_output.output |= csa_output.output << 8;
    final_csa_output.carry  |= csa_output.carry  << 8;

    return final_csa_output;
}

struct RecodedMultiplicands get_recoded_multiplicands(u64 multiplicand, u64 multiplier) {
    struct RecodedMultiplicands recoded_multiplicands;

    for (int i = 0; i < 4; i++) {
        recoded_multiplicands.m[i] = booth_recode(multiplicand, (multiplier >> (2 * i)) & 0b111);
    }

    return recoded_multiplicands;
}

struct CSAOutput perform_one_cycle_of_booths_mutliplication(struct CSAOutput previous_output, u64 multiplicand, u64 multiplier) {
    struct RecodedMultiplicands recoded_multiplicands = get_recoded_multiplicands(multiplicand, multiplier);
    return perform_csa_array(previous_output.output, previous_output.carry, recoded_multiplicands);
}

enum MultiplicationFlavor {
    SHORT,
    LONG_SIGNED,
    LONG_UNSIGNED,
};

bool is_long(enum MultiplicationFlavor flavor) {
    return flavor == LONG_SIGNED || flavor == LONG_UNSIGNED;
}

bool is_signed(enum MultiplicationFlavor flavor) {
    return flavor == LONG_SIGNED || flavor == SHORT;
}

bool should_terminate(u64 multiplier, enum MultiplicationFlavor flavor) {
    if (is_signed(flavor)) {
        return multiplier == 0x1FFFFFFFF || multiplier == 0;
    } else {
        return multiplier == 0;
    }
}

struct AdderOutput {
    u32 output;
    bool carry;
};

struct AdderOutput adder(u32 a, u32 b, bool carry) {
    u32 output = a + b + carry;
    bool overflow = output < a || output < b;
    return (struct AdderOutput) { output, overflow };
}

struct MultiplicationOutput {
    u64 output;
    bool carry;
};

struct u128 {
    u64 lo;
    u64 hi;
};

struct u128 u128_ror(struct u128 input, int shift) {
    return (struct u128) {
        (input.lo >> shift) | (input.hi << (64 - shift)),
        (input.hi >> shift) | (input.lo << (64 - shift)),
    };
}

struct MultiplicationOutput booths_multiplication(enum MultiplicationFlavor flavor, u64 multiplicand, u64 multiplier, u64 accumulator) {
    struct CSAOutput csa_output = { 0, 0 };

    bool alu_carry_in = multiplier & 1;

    if (is_signed(flavor)) {
        multiplier = sign_extend(multiplier, 32, 34);
    } else {
        multiplier = multiplier & 0x1FFFFFFFFull;
    }

    if (is_signed(flavor)) {
        multiplicand = sign_extend(multiplicand, 32, 34);
    } else {
        multiplicand = multiplicand & 0x1FFFFFFFFull;
    }

    csa_output.carry = (multiplier & 1) ? ~(multiplicand) : 0;
    csa_output.output = accumulator;
    acc_shift_register = accumulator >> 34;

    struct u128 partial_sum   = { 0, 0 };
    struct u128 partial_carry = { 0, 0 };
    partial_sum.lo   = csa_output.output & 1;
    partial_carry.lo = csa_output.carry  & 1;

    csa_output.output >>= 1;
    csa_output.carry >>= 1;
    partial_sum   = u128_ror(partial_sum, 1);
    partial_carry = u128_ror(partial_carry, 1);

    int num_iterations = 0;
    do {
        csa_output = perform_one_cycle_of_booths_mutliplication(csa_output, multiplicand, multiplier);

        partial_sum.lo   |= csa_output.output & 0xFF;
        partial_carry.lo |= csa_output.carry  & 0xFF;

        csa_output.output >>= 8;
        csa_output.carry >>= 8;

        partial_sum = u128_ror(partial_sum, 8);
        partial_carry = u128_ror(partial_carry, 8);

        multiplier = asr(multiplier, 8, 33);
        num_iterations++;
    } while (!should_terminate(multiplier, flavor));
    partial_sum.lo |= csa_output.output;
    partial_carry.lo |= csa_output.carry;

    // we have ror'd partial_sum and partial_carry by 8 * num_iterations + 1
    // we now need to ror backwards, i tried my best to mimic the table, but
    // i'm off by one for whatever reason.
    int correction_ror;
    if (num_iterations == 1) correction_ror = 23;
    if (num_iterations == 2) correction_ror = 15;
    if (num_iterations == 3) correction_ror = 7;
    if (num_iterations == 4) correction_ror = 31;

    partial_sum   = u128_ror(partial_sum, correction_ror);
    partial_carry = u128_ror(partial_carry, correction_ror);

    if (is_long(flavor)) {
        if (num_iterations == 4) {
            struct AdderOutput adder_output_lo = 
                adder(partial_sum.hi, partial_carry.hi, alu_carry_in);
            struct AdderOutput adder_output_hi = 
                adder(partial_sum.hi >> 32, partial_carry.hi >> 32, 
                    adder_output_lo.carry);

            return (struct MultiplicationOutput) {
                ((u64) adder_output_hi.output << 32) | adder_output_lo.output,
                (partial_carry.hi >> 63) & 1
            };
        } else {
            struct AdderOutput adder_output_lo = 
                adder(partial_sum.hi >> 32, partial_carry.hi >> 32, alu_carry_in);

            int shift_amount = 1 + 8 * num_iterations;

            // why this is needed is unknown, but the multiplication doesn't work
            // without it
            shift_amount++;

            partial_carry.lo = sign_extend(partial_carry.lo, shift_amount, 64);
            partial_sum.lo |= acc_shift_register << (shift_amount);

            struct AdderOutput adder_output_hi = 
                adder(partial_sum.lo, partial_carry.lo, adder_output_lo.carry);
            return (struct MultiplicationOutput) { 
                ((u64) adder_output_hi.output << 32) | adder_output_lo.output,
                (partial_carry.hi >> 63) & 1
            };
        }
    } else {
        if (num_iterations == 4) {
            struct AdderOutput adder_output = 
                adder(partial_sum.hi, partial_carry.hi, alu_carry_in);
            return (struct MultiplicationOutput) { 
                adder_output.output,
                (partial_carry.hi >> 31) & 1
            };
        } else {
            struct AdderOutput adder_output = 
                adder(partial_sum.hi >> 32, partial_carry.hi >> 32, alu_carry_in);
            return (struct MultiplicationOutput) { 
                adder_output.output,
                (partial_carry.hi >> 63) & 1
            };
        }
    }
}

struct MultiplicationOutput mul(u32 rm, u32 rs) {
    return booths_multiplication(SHORT, rm, rs, 0);
}

struct MultiplicationOutput mla(u32 rm, u32 rs, u32 rn) {
    return booths_multiplication(SHORT, rm, rs, rn);
}

struct MultiplicationOutput umull(u32 rm, u32 rs) {
    return booths_multiplication(LONG_UNSIGNED, rm, rs, 0);
}

struct MultiplicationOutput umlal(u32 rdlo, u32 rdhi, u32 rm, u32 rs) {
    return booths_multiplication(LONG_UNSIGNED, rm, rs, ((u64) rdhi << 32) | (u64) rdlo);
}

struct MultiplicationOutput smull(u32 rm, u32 rs) {
    return booths_multiplication(LONG_SIGNED, rm, rs, 0);
}

struct MultiplicationOutput smlal(u32 rdlo, u32 rdhi, u32 rm, u32 rs) {
    return booths_multiplication(LONG_SIGNED, rm, rs, (u64) rdhi << 32 | (u64) rdlo);
}