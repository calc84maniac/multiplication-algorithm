
#define PC_BUILD 1
#include "impl.h"
#undef printf

#include <stdio.h>
#include <time.h>
#include <stdlib.h>

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
        case 6: return 0xAAAAAAAA & (rand() & 0xFFFF) | ((rand() & 0xFFFF) << 16);
        case 7: return 0x55555555 & (rand() & 0xFFFF) | ((rand() & 0xFFFF) << 16);
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
            printf("[MUL CARRY REGRESSION] Failed: %llx * %llx = %llx, got %llx\n", multiplicand, multiplier, actual, guess);
            return 1;
        } else {
            printf("[MUL CARRY REGRESSION] Passed: %llx * %llx = %llx, got %llx\n", multiplicand, multiplier, actual, guess);
        }
    }

    printf("Running mul regression tests...\n");
    for (int i = 0; i < 10000; i++) {
        unsigned int multiplicand = get_rand_32();
        unsigned int multiplier = get_rand_32();

        unsigned int guess = mul(multiplicand, multiplier).output;
        unsigned int actual = (unsigned int) multiplicand * (unsigned int) multiplier;

        if (guess != actual) {
            printf("[MUL REGRESSION] Failed: %llx * %llx = %llx, got %llx\n", multiplicand, multiplier, actual, guess);
            return 1;
        } else {
            printf("[MUL REGRESSION] Passed: %llx * %llx = %llx, got %llx\n", multiplicand, multiplier, actual, guess);
        }
    }

    printf("Running mla regression tests...\n");
    for (int i = 0; i < 10000; i++) {
        unsigned int multiplicand = get_rand_32();
        unsigned int multiplier = get_rand_32();
        unsigned int accumulate = get_rand_32();

        unsigned int guess = mla(multiplicand, multiplier, accumulate).output;
        unsigned int actual = (unsigned int) multiplicand * (unsigned int) multiplier + (unsigned int) accumulate;

        if (guess != actual) {
            printf("[MLA REGRESSION] Failed: %llx * %llx + %llx = %llx, got %llx\n", multiplicand, multiplier, accumulate, actual, guess);
            return 1;
        }
    }

    printf("Running umull regression tests...\n");
    for (int i = 0; i < 10000; i++) {
        unsigned int multiplicand = get_rand_32();
        unsigned int multiplier = get_rand_32();
        unsigned int accumulate = get_rand_32();

        unsigned long long guess = umull(multiplicand, multiplier).output;
        unsigned long long actual = (unsigned long long) multiplicand * (unsigned long long) multiplier;

        if (guess != actual) {
            printf("[UMULL REGRESSION] Failed #%d: %llx * %llx = %llx, got %llx\n", i, multiplicand, multiplier, actual, guess);
            return 1;
        } else {
            printf("[UMULL REGRESSION] Passed #%d: %llx * %llx = %llx, got %llx\n", i, multiplicand, multiplier, actual, guess);
        }
    }

    printf("Running umlal regression tests...\n");
    for (int i = 0; i < 10000; i++) {
        unsigned int multiplicand = get_rand_32();
        unsigned int multiplier = get_rand_32();
        unsigned int accumulate = get_rand_32();
        unsigned int accumulate2 = get_rand_32();

        unsigned long long guess = umlal(accumulate, accumulate2, multiplicand, multiplier).output;
        unsigned long long actual_acc = (unsigned long long) accumulate + ((unsigned long long) accumulate2 << 32);
        unsigned long long actual =
            (unsigned long long) multiplicand * (unsigned long long) multiplier + actual_acc;

        if (guess != actual) {
            printf("[UMLAL REGRESSION] Failed: %llx * %llx + %llx = %llx, got %llx\n", multiplicand, multiplier, actual_acc, actual, guess);
            return 1;
        } else {
            printf("[UMLAL REGRESSION] Passed: %llx * %llx + %llx = %llx, got %llx\n", multiplicand, multiplier, actual_acc, actual, guess);
        }
    }

    printf("Running smull regression tests...\n");
    for (int i = 0; i < 10000; i++) {
        unsigned int multiplicand = get_rand_32();
        unsigned int multiplier = get_rand_32();

        long long guess = smull(multiplicand, multiplier).output;
        long long actual = (long long) (int) multiplicand * (long long) (int) multiplier;

        if (guess != actual) {
            printf("[SMULL REGRESSION] Failed: %llx * %llx = %llx, got %llx\n", multiplicand, multiplier, actual, guess);
            return 1;
        } else {
            printf("[SMULL REGRESSION] Passed: %llx * %llx = %llx, got %llx\n", multiplicand, multiplier, actual, guess);
        }
    }

    printf("Running smlal regression tests...\n");
    for (int i = 0; i < 10000; i++) {
        unsigned int multiplicand = get_rand_32();
        unsigned int multiplier = get_rand_32();
        unsigned int accumulate = get_rand_32();
        unsigned int accumulate2 = get_rand_32();

        long long guess = smlal(accumulate, accumulate2, multiplicand, multiplier).output;
        long long actual_acc = (long long) accumulate + ((long long) accumulate2 << 32);
        long long actual = (long long) (int) multiplicand * (long long) (int) multiplier + actual_acc;

        if (guess != actual) {
            printf("[SMLAL REGRESSION] Failed: %llx * %llx + %llx = %llx, got %llx\n", multiplicand, multiplier, actual_acc, actual, guess);
            return 1;
        } else {
            printf("[SMLAL REGRESSION] Passed: %llx * %llx + %llx = %llx, got %llx\n", multiplicand, multiplier, actual_acc, actual, guess);
        }
    }

    printf("All tests passed!\n");
}