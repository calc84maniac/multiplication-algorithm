
#define PC_BUILD 1
#include "impl.h"
#undef printf
#include "impl_opt.h"

#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <inttypes.h>

// generate interesting random numbers
unsigned int get_rand_32() {
    retry:
    switch (rand() & 0b1111) {
        case 0: return rand() & 0xFF;
        case 1: return rand() & 0xFFFF;
        case 2: return rand() & 0xFFFFFF;
        case 3: return 0xFFFFFFFF - (rand() & 0xFF);
        case 4: return 0xFFFFFFFF - (rand() & 0xFFFF);
        case 5: return 0xFFFFFFFF - (rand() & 0xFFFFFF);
        case 6: return (0xAAAAAAAA & (rand() & 0xFFFF)) | ((rand() & 0xFFFF) << 16);
        case 7: return (0x55555555 & (rand() & 0xFFFF)) | ((rand() & 0xFFFF) << 16);
        case 8: return 0;
        case 9: return ((rand() & 0xFFFF) | ((rand() & 0xFFFF) << 16));

        default:
            goto retry;
            // return (rand() & 0xFFFF) | ((rand() & 0xFFFF) << 16);   
    }
}

bool guess_mul_zero(uint32_t x) {
    // the carry flag's behavior for short mul when the multiplicand is 0:
    // used for some crude testing of the mul instruction, before doing full
    // out fuzzing on a GBA
    
    if (x >> 8 == 0xFFFFFF) {
        unsigned int masked_x = x & 0xFF;
        if (masked_x >= 0xC0) return false;
        return (masked_x & 0x55) != 0;
    }

    if (x >> 16 == 0xFFFF) {
        unsigned int masked_x = x & 0xFFFF;
        if (masked_x >= 0xC000) return false;
        return (masked_x & 0x5555) != 0;
    }

    if (x >> 24 == 0xFF) {
        unsigned int masked_x = x & 0xFFFFFF;
        if (masked_x >= 0xC00000) return false;
        return (masked_x & 0x555555) != 0;
    }

    else {
        return (x >> 30) == 2;
    }

}

int main() {
    srand(time(NULL));

    printf("Running mul carry regression tests...\n");
    for (int i = 0; i < 256; i++) {
        unsigned int multiplicand = 0;
        unsigned int multiplier = 0xFFFFFFFF - i;

        unsigned int guess = mul(multiplicand, multiplier).carry;
        unsigned int actual = guess_mul_zero(multiplier);

        if (guess != actual) {
            printf("[MUL CARRY REGRESSION] Failed: %x * %x = %x, got %x\n", multiplicand, multiplier, actual, guess);
            return 1;
        } else {
            printf("[MUL CARRY REGRESSION] Passed: %x * %x = %x, got %x\n", multiplicand, multiplier, actual, guess);
        }
    }

    printf("Running mul regression tests...\n");
    for (int i = 0; i < 10000; i++) {
        unsigned int multiplicand = get_rand_32();
        unsigned int multiplier = get_rand_32();

        struct MultiplicationOutput guess = mul_opt(multiplicand, multiplier);
        struct MultiplicationOutput actual = mul(multiplicand, multiplier);

        if (guess.output != actual.output || guess.carry != actual.carry) {
            printf("[MUL REGRESSION] Failed: %x * %x = %"PRIx64":%x, got %"PRIx64":%x\n", multiplicand, multiplier, actual.output, actual.carry, guess.output, guess.carry);
            return 1;
        } else {
            printf("[MUL REGRESSION] Passed: %x * %x = %"PRIx64":%x\n", multiplicand, multiplier, actual.output, actual.carry);
        }
    }

    printf("Running mla regression tests...\n");
    for (int i = 0; i < 10000; i++) {
        unsigned int multiplicand = get_rand_32();
        unsigned int multiplier = get_rand_32();
        unsigned int accumulate = get_rand_32();

        struct MultiplicationOutput guess = mla_opt(multiplicand, multiplier, accumulate);
        struct MultiplicationOutput actual = mla(multiplicand, multiplier, accumulate);

        if (guess.output != actual.output || guess.carry != actual.carry) {
            printf("[MLA REGRESSION] Failed: %x * %x + %x = %"PRIx64":%x, got %"PRIx64":%x\n", multiplicand, multiplier, accumulate, actual.output, actual.carry, guess.output, guess.carry);
            return 1;
        } else {
            printf("[MLA REGRESSION] Passed: %x * %x + %x = %"PRIx64":%x\n", multiplicand, multiplier, accumulate, actual.output, actual.carry);
        }
    }

    printf("Running umull regression tests...\n");
    for (int i = 0; i < 10000; i++) {
        unsigned int multiplicand = get_rand_32();
        unsigned int multiplier = get_rand_32();

        struct MultiplicationOutput guess = umull_opt(multiplicand, multiplier);
        struct MultiplicationOutput actual = umull(multiplicand, multiplier);

        if (guess.output != actual.output || guess.carry != actual.carry) {
            printf("[UMULL REGRESSION] Failed #%d: %x * %x = %"PRIx64":%x, got %"PRIx64":%x\n", i, multiplicand, multiplier, actual.output, actual.carry, guess.output, guess.carry);
            return 1;
        } else {
            printf("[UMULL REGRESSION] Passed #%d: %x * %x = %"PRIx64":%x", i, multiplicand, multiplier, actual.output, actual.carry);
        }
    }

    printf("Running umlal regression tests...\n");
    for (int i = 0; i < 10000; i++) {
        unsigned int multiplicand = get_rand_32();
        unsigned int multiplier = get_rand_32();
        unsigned int accumulate = get_rand_32();
        unsigned int accumulate2 = get_rand_32();

        struct MultiplicationOutput guess = umlal_opt(accumulate, accumulate2, multiplicand, multiplier);
        struct MultiplicationOutput actual = umlal(accumulate, accumulate2, multiplicand, multiplier);
        u64 actual_acc = (u64) accumulate + ((u64) accumulate2 << 32);

        if (guess.output != actual.output || guess.carry != actual.carry) {
            printf("[UMLAL REGRESSION] Failed: %x * %x + %"PRIx64" = %"PRIx64":%x, got %"PRIx64":%x\n", multiplicand, multiplier, actual_acc, actual.output, actual.carry, guess.output, guess.carry);
            return 1;
        } else {
            printf("[UMLAL REGRESSION] Passed: %x * %x + %"PRIx64" = %"PRIx64":%x\n", multiplicand, multiplier, actual_acc, actual.output, actual.carry);
        }
    }

    printf("Running smull regression tests...\n");
    for (int i = 0; i < 10000; i++) {
        unsigned int multiplicand = get_rand_32();
        unsigned int multiplier = get_rand_32();

        struct MultiplicationOutput guess = smull_opt(multiplicand, multiplier);
        struct MultiplicationOutput actual = smull(multiplicand, multiplier);

        if (guess.output != actual.output || guess.carry != actual.carry) {
            printf("[SMULL REGRESSION] Failed: %x * %x = %"PRIx64":%x, got %"PRIx64":%x\n", multiplicand, multiplier, actual.output, actual.carry, guess.output, guess.carry);
            return 1;
        } else {
            printf("[SMULL REGRESSION] Passed: %x * %x = %"PRIx64":%x\n", multiplicand, multiplier, actual.output, actual.carry);
        }
    }

    printf("Running smlal regression tests...\n");
    for (int i = 0; i < 10000; i++) {
        unsigned int multiplicand = get_rand_32();
        unsigned int multiplier = get_rand_32();
        unsigned int accumulate = get_rand_32();
        unsigned int accumulate2 = get_rand_32();

        struct MultiplicationOutput guess = smlal_opt(accumulate, accumulate2, multiplicand, multiplier);
        struct MultiplicationOutput actual = smlal(accumulate, accumulate2, multiplicand, multiplier);
        u64 actual_acc = (u64) accumulate + ((u64) accumulate2 << 32);

        if (guess.output != actual.output || guess.carry != actual.carry) {
            printf("[SMLAL REGRESSION] Failed: %x * %x + %"PRIx64" = %"PRIx64":%x, got %"PRIx64":%x\n", multiplicand, multiplier, actual_acc, actual.output, actual.carry, guess.output, guess.carry);
            return 1;
        } else {
            printf("[SMLAL REGRESSION] Passed: %x * %x + %"PRIx64" = %"PRIx64":%x\n", multiplicand, multiplier, actual_acc, actual.output, actual.carry);
        }
    }

    printf("All tests passed!\n");
    return 0;
}