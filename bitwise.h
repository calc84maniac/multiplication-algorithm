#ifndef BITWISE_H
#define BITWISE_H

#include <stdbool.h>
#include <stdint.h>

typedef uint64_t    u64;
typedef uint32_t    u32;
typedef uint16_t    u16;
typedef uint8_t     u8;

typedef int64_t     s64;
typedef int32_t     s32;

static inline u64 mask(int lo, int hi) {
    u64 size = hi - lo;
    return ((1ull << size) - 1) << lo;
}

static inline bool bit(u64 x, int n) {
    return (x >> n) & 1;
}

static inline u64 sign_extend(u64 value, int from_size, int to_size) {
    bool negative = bit(value, from_size - 1);
    if (negative) value |= mask(from_size, to_size);
    return value;
};

static inline u64 asr(u64 value, int shift, int size) {
    s64 sign_extended = (s64) sign_extend(value, size, 64);
    sign_extended >>= shift;
    return sign_extended & mask(0, size);
}

#endif