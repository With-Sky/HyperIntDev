#ifndef HINT_HPP
#define HINT_HPP

#include <algorithm>
#include <iostream>
#include <random>
#include <string>
#include <vector>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include <cassert>
#include <climits>
#include <cstdint>

#include <chrono>

#if defined(_WIN64) // Windows MSVC X64
#include <intrin.h>
#define HINT_WIN64
#endif

#if defined(__x86_64__) || defined(__amd64__) // None MSVC X64
#include <x86gprintrin.h>
#define HINT_X86_64
#endif

#if defined(__SIZEOF_INT128__) // GCC int128
#define HINT_INT128
#endif

#if defined(__GNUC__) || defined(__clang__)
#define HINT_ASM
#endif

#if defined(HINT_WIN64) || defined(HINT_X86_64)
#define HINT_ADC64
#endif

#if defined(HINT_WIN64) || defined(HINT_INT128)
#define HINT_MUL64X64
#endif

namespace hint
{
    using HintULL = unsigned long long;
    static_assert(sizeof(HintULL) == sizeof(uint64_t), "HintULL is not 64bits");

    template <typename IntTy>
    constexpr bool is_2pow(IntTy n)
    {
        return (n > 0) && (n & (n - 1)) == 0;
    }

    /// @brief Floor round to the nearest power of 2
    /// @tparam T
    /// @param n Integer that not negative
    /// @return Power of 2 that is not larger than n
    template <typename T>
    constexpr T int_floor2(T n)
    {
        constexpr int bits = sizeof(n) * 8;
        for (int i = 1; i < bits; i *= 2)
        {
            n |= (n >> i);
        }
        return (n >> 1) + 1;
    }

    /// @brief Ceiling round to the nearest power of 2
    /// @tparam T
    /// @param n Integer that not negative
    /// @return Power of 2 that is not smaller than n
    template <typename T>
    constexpr T int_ceil2(T n)
    {
        constexpr int bits = sizeof(n) * 8;
        if (n > 0)
        {
            n--;
        }
        for (int i = 1; i < bits; i *= 2)
        {
            n |= (n >> i);
        }
        return n + 1;
    }

    template <typename T>
    constexpr T all_one(int bits)
    {
        T temp = T(1) << (bits - 1);
        return temp - 1 + temp;
    }
    constexpr int hint_ctz(uint32_t x)
    {
        int r0 = 16, r1 = 8, r2 = 4, r3 = 2, r4 = 1;
        x &= (0 - x);
        if (x & 0x0000FFFF)
        {
            r0 = 0;
        }
        if (x & 0x00FF00FF)
        {
            r1 = 0;
        }
        if (x & 0x0F0F0F0F)
        {
            r2 = 0;
        }
        if (x & 0x33333333)
        {
            r3 = 0;
        }
        if (x & 0x55555555)
        {
            r4 = 0;
        }
        return r0 + r1 + r2 + r3 + r4 + (x == 0);
    }

    constexpr int hint_ctz(uint64_t x)
    {
        if (uint32_t(x))
        {
            return hint_ctz(uint32_t(x));
        }
        return hint_ctz(uint32_t(x >> 32)) + 32;
    }

    // Leading zeros
    constexpr int hint_clz(uint32_t x)
    {
        constexpr uint32_t MASK32 = uint32_t(0xFFFF) << 16;
        int res = sizeof(uint32_t) * CHAR_BIT;
        if (x & MASK32)
        {
            res -= 16;
            x >>= 16;
        }
        if (x & (MASK32 >> 8))
        {
            res -= 8;
            x >>= 8;
        }
        if (x & (MASK32 >> 12))
        {
            res -= 4;
            x >>= 4;
        }
        if (x & (MASK32 >> 14))
        {
            res -= 2;
            x >>= 2;
        }
        if (x & (MASK32 >> 15))
        {
            res -= 1;
            x >>= 1;
        }
        return res - x;
    }
    // Leading zeros
    constexpr int hint_clz(uint64_t x)
    {
        if (x & (uint64_t(0xFFFFFFFF) << 32))
        {
            return hint_clz(uint32_t(x >> 32));
        }
        return hint_clz(uint32_t(x)) + 32;
    }
    template <typename T>
    constexpr T hint_gcd(T a, T b)
    {
        if (a < b)
        {
            std::swap(a, b);
        }
        if (0 == b)
        {
            return a;
        }
        const int i = hint::hint_ctz(a);
        a >>= i;
        const int j = hint::hint_ctz(b);
        b >>= j;
        const int k = std::min(i, j);
        while (true)
        {
            if (a < b)
            {
                std::swap(a, b);
            }
            if (0 == b)
            {
                break;
            }
            a -= b;
            a >>= hint::hint_ctz(a);
        }
        return a << k;
    }
    template <typename T>
    constexpr T hint_lcm(T a, T b)
    {
        return a / hint::hint_gcd(a, b) * b;
    }
    template <typename IntTy>
    constexpr IntTy exgcd(IntTy a, IntTy b, IntTy &x, IntTy &y)
    {
        if (0 == b)
        {
            x = 1;
            y = 0;
            return a;
        }
        IntTy k = a / b;
        IntTy g = exgcd(b, a - k * b, y, x);
        y -= k * x;
        return g;
    }

    template <typename IntTy>
    inline IntTy exgcd_iter(IntTy a, IntTy b, IntTy &x, IntTy &y)
    {
        auto exec = [](IntTy &m, IntTy &n, IntTy q)
        {
            IntTy temp = m - n * q;
            m = n, n = temp;
        };
        x = 1, y = 0;
        IntTy x1 = 0, y1 = 1;
        while (b > 0)
        {
            IntTy q = a / b;
            exec(x, x1, q);
            exec(y, y1, q);
            exec(a, b, q);
        }
        return a;
    }
    // return n^-1 mod mod
    template <typename IntTy>
    constexpr IntTy mod_inv(IntTy n, IntTy mod)
    {
        n %= mod;
        IntTy x = 0, y = 0;
        exgcd(n, mod, x, y);
        if (x < 0)
        {
            x += mod;
        }
        else if (x >= mod)
        {
            x -= mod;
        }
        return x;
    }
    // Integer bit length
    template <typename IntTy>
    constexpr int hint_bit_length(IntTy x)
    {
        if (0 == x)
        {
            return 0;
        }
        return sizeof(IntTy) * CHAR_BIT - hint_clz(x);
    }

    // Integer log2
    template <typename IntTy>
    constexpr int hint_log2(IntTy x)
    {
        return (sizeof(IntTy) * CHAR_BIT - 1) - hint_clz(x);
    }
    // return n^-1 mod 2^pow, Newton iteration
    constexpr uint64_t inv_mod2pow(uint64_t n, int pow)
    {
        const uint64_t mask = all_one<uint64_t>(pow);
        uint64_t xn = 1, t = n & mask;
        while (t != 1)
        {
            xn = (xn * (2 - t));
            t = (xn * n) & mask;
        }
        return xn & mask;
    }

    template <typename T1, typename T2>
    bool mem_overlap(T1 begin1, T1 end1, T2 begin2, T2 end2)
    {
        return begin1 <= end2 && begin2 <= end1;
    }
    // Fast power
    template <typename T, typename T1>
    constexpr T qpow(T m, T1 n)
    {
        T result = 1;
        while (true)
        {
            if (n & 1)
            {
                result *= m;
            }
            if (0 == n)
            {
                break;
            }
            m *= m;
            n >>= 1;
        }
        return result;
    }
    // Fast power with mod
    template <typename T, typename T1>
    constexpr T qpow(T m, T1 n, T mod)
    {
        T result = 1;
        while (true)
        {
            if (n & 1)
            {
                result *= m;
                result %= mod;
            }
            if (0 == n)
            {
                break;
            }
            m *= m;
            m %= mod;
            n >>= 1;
        }
        return result;
    }
    template <int BITS>
    struct Uint
    {
        static constexpr int MAX_BITS = 64;
        static constexpr int MIN_BITS = 8;
        static constexpr int NORM_BITS = hint::int_ceil2(std::min(MAX_BITS, std::max(MIN_BITS, BITS)));
        static constexpr int NEXT_BITS = std::min(MAX_BITS, NORM_BITS * 2);
        static constexpr int LAST_BITS = std::max(MIN_BITS, NORM_BITS / 2);
        using SignType = typename Uint<NORM_BITS>::SignType;
        using Type = typename Uint<NORM_BITS>::Type;
        using NextType = Uint<NEXT_BITS>;
        using LastType = Uint<LAST_BITS>;
    };

    template <int BITS>
    constexpr int Uint<BITS>::MAX_BITS;
    template <int BITS>
    constexpr int Uint<BITS>::MIN_BITS;
    template <int BITS>
    constexpr int Uint<BITS>::NORM_BITS;
    template <int BITS>
    constexpr int Uint<BITS>::NEXT_BITS;
    template <int BITS>
    constexpr int Uint<BITS>::LAST_BITS;

    template <>
    struct Uint<8>
    {
        using SignType = int8_t;
        using Type = uint8_t;
    };

    template <>
    struct Uint<16>
    {
        using SignType = int16_t;
        using Type = uint16_t;
    };

    template <>
    struct Uint<32>
    {
        using SignType = int32_t;
        using Type = uint32_t;
    };

    template <>
    struct Uint<64>
    {
        using SignType = int64_t;
        using Type = uint64_t;
    };

    // template <int BITS>
    // struct UintType
    // {
    //     static constexpr int MAX_BITS = 128;
    //     static constexpr int MIN_BITS = 8;
    //     static constexpr int NORM_BITS = hint::int_ceil2(std::min(MAX_BITS, std::max(MIN_BITS, BITS)));
    //     static constexpr int NEXT_BITS = std::min(MAX_BITS, NORM_BITS * 2);
    //     static constexpr int LAST_BITS = std::max(MIN_BITS, NORM_BITS / 2);

    //     using SignType = typename Uint<NORM_BITS>::SignType;
    //     using Type = typename Uint<NORM_BITS>::Type;

    //     using NextType = UintType<NEXT_BITS>;
    //     using LastType = UintType<LAST_BITS>;
    // };
    // template <int BITS>
    // constexpr int UintType<BITS>::MAX_BITS;
    // template <int BITS>
    // constexpr int UintType<BITS>::MIN_BITS;
    // template <int BITS>
    // constexpr int UintType<BITS>::NORM_BITS;
    // template <int BITS>
    // constexpr int UintType<BITS>::NEXT_BITS;
    // template <int BITS>
    // constexpr int UintType<BITS>::LAST_BITS;

    namespace utility
    {
        template <typename T>
        void ary_print(const T a[], size_t len, bool is_rev = false)
        {
            if (0 == len)
            {
                std::cout << "0:[]\n";
                return;
            }
            std::cout << (is_rev ? "h->l " : "l->h ");
            std::cout << len << ":[";
            if (is_rev)
            {
                std::cout << *(a + len - 1);
                for (auto it = a + len - 1; it > a;)
                {
                    std::cout << ", " << *--it;
                }
            }
            else
            {
                std::cout << *a;
                for (auto it = a + 1; it < a + len; ++it)
                {
                    std::cout << ", " << *it;
                }
            }
            std::cout << "]" << std::endl;
        }
        size_t get_add_len(size_t l_len, size_t r_len)
        {
            return std::max(l_len, r_len) + 1;
        }

        size_t get_sub_len(size_t l_len, size_t r_len)
        {
            return std::max(l_len, r_len);
        }

        size_t get_mul_len(size_t l_len, size_t r_len)
        {
            if (0 == l_len || 0 == r_len)
            {
                return 0;
            }
            return l_len + r_len;
        }

        size_t get_div_len(size_t l_len, size_t r_len)
        {
            assert(r_len > 0);
            if (l_len < r_len)
            {
                return 0;
            }
            return l_len - r_len + 1;
        }
        // x + y = sum with carry
        template <typename UintTy>
        constexpr UintTy add_half(UintTy x, UintTy y, bool &cf)
        {
            x = x + y;
            cf = (x < y);
            return x;
        }

        // x - y = diff with borrow
        template <typename UintTy>
        constexpr UintTy sub_half(UintTy x, UintTy y, bool &bf)
        {
            y = x - y;
            bf = (y > x);
            return y;
        }

        // x + y + cf = sum with carry
        template <typename UintTy>
        constexpr UintTy add_carry(UintTy x, UintTy y, bool &cf)
        {
            UintTy sum = x + cf;
            cf = (sum < x);
            sum += y;             // carry
            cf = cf || (sum < y); // carry
            return sum;
        }

        // x - y - bf = diff with borrow
        template <typename UintTy>
        constexpr UintTy sub_borrow(UintTy x, UintTy y, bool &bf)
        {
            UintTy diff = x - bf;
            bf = (diff > x);
            y = diff - y;          // borrow
            bf = bf || (y > diff); // borrow
            return y;
        }

        template <typename Ui64>
        constexpr void mul64x64to128_base(uint64_t a, uint64_t b, Ui64 &low, Ui64 &high)
        {
            static_assert(sizeof(Ui64) == sizeof(uint64_t), "mul64x64to128_base: low and high must be 64bit");
            uint64_t ah = a >> 32, bh = b >> 32;
            a = uint32_t(a), b = uint32_t(b);
            uint64_t r0 = a * b, r1 = a * bh, r2 = ah * b, r3 = ah * bh;
            r3 += (r1 >> 32) + (r2 >> 32);
            r1 = uint32_t(r1), r2 = uint32_t(r2);
            r1 += r2;
            r1 += (r0 >> 32);
            high = r3 + (r1 >> 32);
            low = (r1 << 32) | uint32_t(r0);
        }

        template <typename Ui64>
        inline void mul64x64to128(uint64_t a, uint64_t b, Ui64 &low, Ui64 &high)
        {
            static_assert(sizeof(Ui64) == sizeof(uint64_t), "mul64x64to128: low and high must be 64bit");
#if defined(HINT_INT128) // Has __uint128_t
#pragma message("Using __uint128_t to compute 64bit x 64bit to 128bit")
            __uint128_t x = __uint128_t(a) * b;
            low = uint64_t(x), high = uint64_t(x >> 64);
#elif defined(HINT_WIN64) // Has _umul128
#pragma message("Using _umul128 to compute 64bit x 64bit to 128bit")
            HintULL lo, hi;
            lo = _umul128(a, b, &hi);
            low = lo, high = hi;
#else // No _umul128 or __uint128_t
#pragma message("Using basic function to compute 64bit x 64bit to 128bit")
            mul64x64to128_base(a, b, low, high);
#endif
        }

        template <typename UintTy>
        inline void mul_binary(UintTy a, UintTy b, UintTy &low, UintTy &high)
        {
            constexpr int BITS = sizeof(UintTy) * CHAR_BIT;
            uint64_t prod = uint64_t(a) * uint64_t(b);
            low = UintTy(prod), high = UintTy(prod >> BITS);
        }
        inline void mul_binary(uint64_t a, uint64_t b, uint64_t &low, uint64_t &high)
        {
            mul64x64to128(a, b, low, high);
        }

        constexpr uint32_t div128by32(uint64_t &dividend_hi64, uint64_t &dividend_lo64, uint32_t divisor)
        {
            uint32_t quot_hi32 = 0, quot_lo32 = 0;
            uint64_t dividend = dividend_hi64 >> 32;
            quot_hi32 = dividend / divisor;
            dividend %= divisor;

            dividend = (dividend << 32) | uint32_t(dividend_hi64);
            quot_lo32 = dividend / divisor;
            dividend %= divisor;
            dividend_hi64 = (uint64_t(quot_hi32) << 32) | quot_lo32;

            dividend = (dividend << 32) | uint32_t(dividend_lo64 >> 32);
            quot_hi32 = dividend / divisor;
            dividend %= divisor;

            dividend = (dividend << 32) | uint32_t(dividend_lo64);
            quot_lo32 = dividend / divisor;
            dividend %= divisor;
            dividend_lo64 = (uint64_t(quot_hi32) << 32) | quot_lo32;
            return dividend;
        }

        // 96bit integer divided by 64bit integer, input make sure the quotient smaller than 2^32.
        constexpr uint32_t div96by64to32(uint32_t dividend_hi32, uint64_t &dividend_lo64, uint64_t divisor)
        {
            if (0 == dividend_hi32)
            {
                uint32_t quotient = dividend_lo64 / divisor;
                dividend_lo64 %= divisor;
                return quotient;
            }
            uint64_t divid2 = (uint64_t(dividend_hi32) << 32) | (dividend_lo64 >> 32);
            uint64_t divis1 = divisor >> 32;
            divisor = uint32_t(divisor);
            uint64_t qhat = divid2 / divis1;
            divid2 %= divis1;
            divid2 = (divid2 << 32) | uint32_t(dividend_lo64);
            uint64_t prod = qhat * divisor;
            divis1 <<= 32;
            if (prod > divid2)
            {
                qhat--;
                prod -= divisor;
                divid2 += divis1;
                // if divid2 <= divis1, the addtion of divid2 is overflow, so prod must not be larger than divid2.
                if ((divid2 > divis1) && (prod > divid2)) [[unlikely]]
                {
                    qhat--;
                    prod -= divisor;
                    divid2 += divis1;
                }
            }
            dividend_lo64 = divid2 - prod;
            return uint32_t(qhat);
        }
        // 128bit integer divided by 64bit integer, input make sure the quotient smaller than 2^64.
        constexpr uint64_t div128by64to64(uint64_t dividend_hi64, uint64_t &dividend_lo64, uint64_t divisor)
        {
            int k = 0;
            if (divisor < (uint64_t(1) << 63))
            {
                k = hint::hint_clz(divisor);
                divisor <<= k; // Normalization.
                dividend_hi64 = (dividend_hi64 << k) | (dividend_lo64 >> (64 - k));
                dividend_lo64 <<= k;
            }
            uint32_t divid_hi32 = dividend_hi64 >> 32;
            uint64_t divid_lo64 = (dividend_hi64 << 32) | (dividend_lo64 >> 32);
            uint64_t quotient = div96by64to32(divid_hi32, divid_lo64, divisor);

            divid_hi32 = divid_lo64 >> 32;
            dividend_lo64 = uint32_t(dividend_lo64) | (divid_lo64 << 32);
            quotient = (quotient << 32) | div96by64to32(divid_hi32, dividend_lo64, divisor);
            dividend_lo64 >>= k;
            return quotient;
        }

        constexpr uint64_t div_dword_word(uint64_t dividend_hi, uint64_t dividend_lo, uint64_t divisor, uint64_t &remainder)
        {
            uint64_t quot = div128by64to64(dividend_hi, dividend_lo, divisor);
            remainder = dividend_lo;
            return quot;
        }

        template <typename UintTy>
        constexpr UintTy div_dword_word(UintTy dividend_hi, UintTy dividend_lo, UintTy divisor, UintTy &remainder)
        {
            constexpr int BITS = sizeof(UintTy) * CHAR_BIT;
            uint64_t dividend = (uint64_t(dividend_hi) << BITS) | dividend_lo;
            remainder = dividend % divisor;
            return dividend / divisor;
        }

        constexpr uint64_t mulMod(uint64_t a, uint64_t b, uint64_t mod)
        {
            uint64_t prod_lo = 0, prod_hi = 0;
            mul64x64to128_base(a, b, prod_lo, prod_hi);
            div128by64to64(prod_hi, prod_lo, mod);
            return prod_lo;
        }
        template <typename UintTy>
        constexpr UintTy mulMod(UintTy a, UintTy b, UintTy mod)
        {
            return uint64_t(a) * b % mod;
        }
        template <typename UintTy>
        constexpr void add_carry_x4_basic(UintTy in1[4], UintTy in2[4], UintTy sum[4], bool &carry)
        {
            sum[0] = add_carry(in1[0], in2[0], carry);
            sum[1] = add_carry(in1[1], in2[1], carry);
            sum[2] = add_carry(in1[2], in2[2], carry);
            sum[3] = add_carry(in1[3], in2[3], carry);
        }

        template <typename UintTy>
        inline void add_carry_x4(const UintTy in1[4], const UintTy in2[4], UintTy sum[4], bool &carry)
        {
            static_assert(std::is_integral<UintTy>::value, "Input type must be integral type");
            add_carry_x4_basic(in1, in2, sum, carry);
        }

        inline void add_carry_x4(const uint64_t in1[4], const uint64_t in2[4], uint64_t sum[4], bool &carry)
        {
#if defined(HINT_ADC64)
#pragma message("Using _addcarry_u64 to compute 256bit + 256bit")
            unsigned char cf = carry;
            auto p = (HintULL *)sum;
            cf = _addcarry_u64(cf, in1[0], in2[0], p);
            cf = _addcarry_u64(cf, in1[1], in2[1], p + 1);
            cf = _addcarry_u64(cf, in1[2], in2[2], p + 2);
            carry = _addcarry_u64(cf, in1[3], in2[3], p + 3);
#else
#pragma message("Using basic function to compute 256bit + 256bit")
            add_carry_x4_basic(in1, in2, sum, carry);
#endif
        }

        template <typename UintTy>
        constexpr void sub_borrow_x4_basic(const UintTy in1[4], const UintTy in2[4], UintTy diff[4], bool &borrow)
        {
            diff[0] = sub_borrow(in1[0], in2[0], borrow);
            diff[1] = sub_borrow(in1[1], in2[1], borrow);
            diff[2] = sub_borrow(in1[2], in2[2], borrow);
            diff[3] = sub_borrow(in1[3], in2[3], borrow);
        }
        template <typename UintTy>
        inline void sub_borrow_x4(const UintTy in1[4], const UintTy in2[4], UintTy diff[4], bool &borrow)
        {
            static_assert(std::is_integral<UintTy>::value, "Input type must be integral type");
            sub_borrow_x4_basic(in1, in2, diff, borrow);
        }

        inline void sub_borrow_x4(const uint64_t in1[4], const uint64_t in2[4], uint64_t diff[4], bool &borrow)
        {
#if defined(HINT_ADC64)
#pragma message("Using _subborrow_u64 to compute 256bit - 256bit")
            unsigned char bf = borrow;
            auto p = (HintULL *)diff;
            bf = _subborrow_u64(bf, in1[0], in2[0], p);
            bf = _subborrow_u64(bf, in1[1], in2[1], p + 1);
            bf = _subborrow_u64(bf, in1[2], in2[2], p + 2);
            borrow = _subborrow_u64(bf, in1[3], in2[3], p + 3);
#else
#pragma message("Using basic function to compute 256bit - 256bit")
            sub_borrow_x4_basic(in1, in2, diff, borrow);
#endif
        }
        template <typename UintTy>
        constexpr void add_cf_x4_basic(const UintTy in[4], UintTy sum[4], bool &carry)
        {
            sum[0] = add_half(in[0], UintTy(carry), carry);
            sum[1] = add_half(in[1], UintTy(carry), carry);
            sum[2] = add_half(in[2], UintTy(carry), carry);
            sum[3] = add_half(in[3], UintTy(carry), carry);
        }
        template <typename UintTy>
        inline void add_cf_x4(const UintTy in[4], UintTy sum[4], bool &carry)
        {
            static_assert(std::is_integral<UintTy>::value, "Input type must be integral type");
            add_cf_x4_basic(in, sum, carry);
        }
        inline void add_cf_x4(const uint64_t in[4], uint64_t sum[4], bool &carry)
        {
#if defined(HINT_ADC64)
#pragma message("Using _addcarry_u64 to compute 256bit + carry")
            unsigned char cf = carry;
            auto p = (HintULL *)sum;
            cf = _addcarry_u64(cf, in[0], 0, p);
            cf = _addcarry_u64(cf, in[1], 0, p + 1);
            cf = _addcarry_u64(cf, in[2], 0, p + 2);
            carry = _addcarry_u64(cf, in[3], 0, p + 3);
#else
#pragma message("Using basic function to compute 256bit + carry")
            add_cf_x4_basic(in, sum, carry);
#endif
        }
        template <typename UintTy>
        constexpr void sub_bf_x4_basic(const UintTy in[4], UintTy diff[4], bool &borrow)
        {
            diff[0] = sub_half(in[0], UintTy(borrow), borrow);
            diff[1] = sub_half(in[1], UintTy(borrow), borrow);
            diff[2] = sub_half(in[2], UintTy(borrow), borrow);
            diff[3] = sub_half(in[3], UintTy(borrow), borrow);
        }
        template <typename UintTy>
        inline void sub_bf_x4(const UintTy in[4], UintTy diff[4], bool &borrow)
        {
            static_assert(std::is_integral<UintTy>::value, "Input type must be integral type");
            sub_bf_x4_basic(in, diff, borrow);
        }
        inline void sub_bf_x4(const uint64_t in[4], uint64_t diff[4], bool &borrow)
        {
#if defined(HINT_ADC64)
#pragma message("Using _subborrow_u64 to compute 256bit - borrow")
            unsigned char bf = borrow;
            auto p = (HintULL *)diff;
            bf = _subborrow_u64(bf, in[0], 0, p);
            bf = _subborrow_u64(bf, in[1], 0, p + 1);
            bf = _subborrow_u64(bf, in[2], 0, p + 2);
            borrow = _subborrow_u64(bf, in[3], 0, p + 3);
#else
#pragma message("Using basic function to compute 256bit - borrow")
            sub_bf_x4_basic(in, diff, borrow);
#endif
        }

        template <typename UintTy>
        constexpr void mul_add(UintTy a, UintTy b, UintTy c, UintTy &lo, UintTy &hi)
        {
            constexpr int NUM_BITS = sizeof(UintTy) * CHAR_BIT;
            static_assert(NUM_BITS <= 32, "Only for 32bit or lower");
            using ProdTy = typename Uint<NUM_BITS>::NextType::Type;
            ProdTy prod = ProdTy(a) * b;
            prod += c;
            hi = UintTy(prod >> NUM_BITS);
            lo = UintTy(prod);
        }

        // [hi,lo] = a * b + c
        inline void mul_add(uint64_t a, uint64_t b, uint64_t c, uint64_t &lo, uint64_t &hi)
        {
            mul64x64to128(a, b, lo, hi);
            lo += c;
            hi += (lo < c);
        }

        template <typename UintTy>
        inline UintTy mul_add_x4_basic(const UintTy in[4], UintTy num_mul, UintTy num_add, UintTy out[4])
        {
            mul_add(in[0], num_mul, num_add, out[0], num_add);
            mul_add(in[1], num_mul, num_add, out[1], num_add);
            mul_add(in[2], num_mul, num_add, out[2], num_add);
            mul_add(in[3], num_mul, num_add, out[3], num_add);
            return num_add;
        }
        template <typename UintTy>
        inline UintTy mul_add_x4(const UintTy in[4], UintTy num_mul, UintTy num_add, UintTy out[4])
        {
            static_assert(std::is_integral<UintTy>::value, "Input type must be integral type");
            return mul_add_x4_basic(in, num_mul, num_add, out);
        }

        inline uint64_t mul_add_x4_basic(const uint64_t in[4], uint64_t num_mul, uint64_t num_add, uint64_t out[4])
        {
#if defined(HINT_ADC64) && defined(HINT_MUL64X64)
#pragma message("Using _subborrow_u64 to compute 256bit * num_mul + num_add")
            unsigned char cf = 0;
            auto p = (HintULL *)out;
            HintULL lo0, hi0, lo1, hi1, lo2, hi2, lo3, hi3;
            mul64x64to128(in[0], num_mul, lo0, hi0);
            mul64x64to128(in[1], num_mul, lo1, hi1);
            mul64x64to128(in[2], num_mul, lo2, hi2);
            mul64x64to128(in[3], num_mul, lo3, hi3);
            cf = _addcarry_u64(0, lo0, num_add, p);
            cf = _addcarry_u64(cf, lo1, hi0, p + 1);
            cf = _addcarry_u64(cf, lo2, hi1, p + 2);
            cf = _addcarry_u64(cf, lo3, hi2, p + 3);
            return hi3 + cf;
#else
#pragma message("Using basic function to compute 256bit * num_mul + num_add")
            return mul_add_x4_basic(in, num_mul, num_add, out);
#endif
        }

        template <typename UintTy>
        inline void mul_add_plus_self(UintTy in, UintTy num_mul, UintTy &out, UintTy &carry)
        {
            bool cf;
            UintTy hi, lo;
            mul_binary(in, num_mul, lo, hi);
            lo = add_half(lo, out, cf);
            hi += cf;
            out = add_half(lo, carry, cf);
            carry = hi + cf;
        }
        template <typename UintTy>
        inline UintTy add_eq_mul_add_x8_basic(const UintTy in[], size_t len, UintTy out[], UintTy num_mul, UintTy num_add)
        {
            len -= len % 8;
            for (size_t i = 0; i < len; i += 8)
            {
                mul_add_plus_self(in[i], num_mul, out[i], num_add);
                mul_add_plus_self(in[i + 1], num_mul, out[i + 1], num_add);
                mul_add_plus_self(in[i + 2], num_mul, out[i + 2], num_add);
                mul_add_plus_self(in[i + 3], num_mul, out[i + 3], num_add);
                mul_add_plus_self(in[i + 4], num_mul, out[i + 4], num_add);
                mul_add_plus_self(in[i + 5], num_mul, out[i + 5], num_add);
                mul_add_plus_self(in[i + 6], num_mul, out[i + 6], num_add);
                mul_add_plus_self(in[i + 7], num_mul, out[i + 7], num_add);
            }
            return num_add;
        }
        template <typename UintTy>
        inline UintTy add_eq_mul_add_x8(const UintTy in[], size_t len, UintTy out[], UintTy num_mul, UintTy num_add)
        {
            static_assert(std::is_integral<UintTy>::value, "Input type must be integral type");
            return add_eq_mul_add_x8_basic(in, len, out, num_mul, num_add);
        }

        inline uint64_t add_eq_mul_add_x8(const uint64_t in[], size_t len, uint64_t out[], uint64_t num_mul, uint64_t num_add)
        {
            len -= len % 8;
#if defined(HINT_ASM) && defined(__ADX__) && defined(__BMI2__)
#pragma message("Using ADX to compute 256bit += 256bit * num_mul + num_add")
            auto end = in + len;
            asm volatile(
                "MOVQ %2, %%R12\n\t"
                "MOVQ %4, %%R13\n\t"
                "LOOP_ADX%=:\n\t"
                "XORQ %%R11, %%R11\n\t"

                "MULX (%%R12), %%RAX, %%R10\n\t"
                "ADCX %0, %%RAX\n\t"
                "ADOX (%%R13), %%RAX\n\t"
                "MOVQ %%RAX, (%%R13)\n\t"

                "MULX 8(%%R12), %%RAX, %0\n\t"
                "ADCX %%R10, %%RAX\n\t"
                "ADOX 8(%%R13), %%RAX\n\t"
                "MOVQ %%RAX, 8(%%R13)\n\t"

                "MULX 16(%%R12), %%RAX, %%R10\n\t"
                "ADCX %0, %%RAX\n\t"
                "ADOX 16(%%R13), %%RAX\n\t"
                "MOVQ %%RAX, 16(%%R13)\n\t"

                "MULX 24(%%R12), %%RAX, %0\n\t"
                "ADCX %%R10, %%RAX\n\t"
                "ADOX 24(%%R13), %%RAX\n\t"
                "MOVQ %%RAX, 24(%%R13)\n\t"

                "MULX 32(%%R12), %%RAX, %%R10\n\t"
                "ADCX %0, %%RAX\n\t"
                "ADOX 32(%%R13), %%RAX\n\t"
                "MOVQ %%RAX, 32(%%R13)\n\t"

                "MULX 40(%%R12), %%RAX, %0\n\t"
                "ADCX %%R10, %%RAX\n\t"
                "ADOX 40(%%R13), %%RAX\n\t"
                "MOVQ %%RAX, 40(%%R13)\n\t"

                "MULX 48(%%R12), %%RAX, %%R10\n\t"
                "ADCX %0, %%RAX\n\t"
                "ADOX 48(%%R13), %%RAX\n\t"
                "MOVQ %%RAX, 48(%%R13)\n\t"

                "MULX 56(%%R12), %%RAX, %0\n\t"
                "ADCX %%R10, %%RAX\n\t"
                "ADOX 56(%%R13), %%RAX\n\t"
                "MOVQ %%RAX, 56(%%R13)\n\t"

                "ADCX %%R11, %0\n\t"
                "ADOX %%R11, %0\n\t"

                "ADDQ $64, %%R12\n\t"
                "ADDQ $64, %%R13\n\t"
                "CMPQ %3,  %%R12\n\t"
                "JB LOOP_ADX%=\n\t"
                : "+r"(num_add)
                : "d"(num_mul), "r"(in), "r"(end), "r"(out)
                : "cc", "memory", "rax", "r10", "r11", "r12", "r13");
            return num_add;
#else
#pragma message("Using basic function to compute 256bit += 256bit * num_mul + num_add")
            return add_eq_mul_add_x8_basic(in, len, out, num_mul, num_add);
#endif
        }

        // uint64_t to std::string
        inline std::string ui64to_string_base10(uint64_t input, uint8_t digits)
        {
            std::string result(digits, '0');
            for (uint8_t i = 0; i < digits; i++)
            {
                result[digits - i - 1] = static_cast<char>(input % 10 + '0');
                input /= 10;
            }
            return result;
        }

        namespace extend_int
        {
            class Uint128
            {
            public:
                Uint128() = default;
                constexpr Uint128(uint64_t l, uint64_t h = 0) : low(l), high(h) {}

                friend constexpr Uint128 operator+(Uint128 lhs, Uint128 rhs)
                {
                    rhs.low += lhs.low;
                    rhs.high += lhs.high + (rhs.low < lhs.low);
                    return rhs;
                }
                friend constexpr Uint128 operator-(Uint128 lhs, Uint128 rhs)
                {
                    rhs.low = lhs.low - rhs.low;
                    rhs.high = lhs.high - rhs.high - (rhs.low > lhs.low);
                    return rhs;
                }
                constexpr Uint128 operator+(uint64_t rhs)
                {
                    rhs = low + rhs;
                    return Uint128(rhs, high + (rhs < low));
                }
                constexpr Uint128 operator-(uint64_t rhs)
                {
                    rhs = low - rhs;
                    return Uint128(rhs, high - (rhs > low));
                }
                // Only compute the low * rhs.low
                friend Uint128 operator*(Uint128 lhs, Uint128 rhs)
                {
                    mul64x64to128(lhs.low, rhs.low, lhs.low, lhs.high);
                    return lhs;
                }
                // Only compute the low * rhs
                Uint128 operator*(uint64_t rhs) const
                {
                    Uint128 res;
                    mul64x64to128(low, rhs, res.low, res.high);
                    return res;
                }
                // Only compute the 128bit / 64 bit
                friend constexpr Uint128 operator/(Uint128 lhs, Uint128 rhs)
                {
                    return lhs / rhs.low;
                }
                // Only compute the 128bit % 64 bit
                friend constexpr Uint128 operator%(Uint128 lhs, Uint128 rhs)
                {
                    return lhs % rhs.low;
                }
                // Only compute the 128bit / 64 bit
                constexpr Uint128 operator/(uint64_t rhs) const
                {
                    Uint128 quot = *this;
                    quot.selfDivRem(rhs);
                    return quot;
                }
                // Only compute the 128bit % 64 bit
                constexpr Uint128 operator%(uint64_t rhs) const
                {
                    Uint128 quot = *this;
                    uint64_t rem = quot.selfDivRem(rhs);
                    return Uint128(rem);
                }
                constexpr Uint128 &operator+=(const Uint128 &rhs)
                {
                    return *this = *this + rhs;
                }
                constexpr Uint128 &operator-=(const Uint128 &rhs)
                {
                    return *this = *this - rhs;
                }
                constexpr Uint128 &operator+=(uint64_t rhs)
                {
                    return *this = *this + rhs;
                }
                constexpr Uint128 &operator-=(uint64_t rhs)
                {
                    return *this = *this - rhs;
                }
                // Only compute the low * rhs.low
                constexpr Uint128 &operator*=(const Uint128 &rhs)
                {
                    mul64x64to128_base(low, rhs.low, low, high);
                    return *this;
                }
                constexpr Uint128 &operator/=(const Uint128 &rhs)
                {
                    return *this = *this / rhs;
                }
                constexpr Uint128 &operator%=(const Uint128 &rhs)
                {
                    return *this = *this % rhs;
                }
                // Return *this % divisor, *this /= divisor
                constexpr uint64_t selfDivRem(uint64_t divisor)
                {
                    if ((divisor >> 32) == 0)
                    {
                        return div128by32(high, low, uint32_t(divisor));
                    }
                    uint64_t divid1 = high % divisor, divid0 = low;
                    high /= divisor;
                    low = div128by64to64(divid1, divid0, divisor);
                    return divid0;
                }
                static constexpr Uint128 mul64x64(uint64_t a, uint64_t b)
                {
                    Uint128 res{};
                    mul64x64to128_base(a, b, res.low, res.high);
                    return res;
                }
                friend constexpr bool operator<(Uint128 lhs, Uint128 rhs)
                {
                    if (lhs.high != rhs.high)
                    {
                        return lhs.high < rhs.high;
                    }
                    return lhs.low < rhs.low;
                }
                friend constexpr bool operator>(Uint128 lhs, Uint128 rhs)
                {
                    return rhs < lhs;
                }
                friend constexpr bool operator<=(Uint128 lhs, Uint128 rhs)
                {
                    return !(lhs > rhs);
                }
                friend constexpr bool operator>=(Uint128 lhs, Uint128 rhs)
                {
                    return !(lhs < rhs);
                }
                friend constexpr bool operator==(Uint128 lhs, Uint128 rhs)
                {
                    return lhs.high == rhs.high && lhs.low == rhs.low;
                }
                friend constexpr bool operator!=(Uint128 lhs, Uint128 rhs)
                {
                    return !(lhs == rhs);
                }
                constexpr Uint128 operator<<(int shift) const
                {
                    if (0 == shift)
                    {
                        return *this;
                    }
                    shift %= 128;
                    shift = shift < 0 ? shift + 128 : shift;
                    if (shift < 64)
                    {
                        return Uint128(low << shift, (high << shift) | (low >> (64 - shift)));
                    }
                    return Uint128(0, low << (shift - 64));
                }
                constexpr Uint128 operator>>(int shift) const
                {
                    if (0 == shift)
                    {
                        return *this;
                    }
                    shift %= 128;
                    shift = shift < 0 ? shift + 128 : shift;
                    if (shift < 64)
                    {
                        return Uint128((low >> shift) | (high << (64 - shift)), high >> shift);
                    }
                    return Uint128(high >> (shift - 64), 0);
                }
                constexpr Uint128 &operator<<=(int shift)
                {
                    return *this = *this << shift;
                }
                constexpr Uint128 &operator>>=(int shift)
                {
                    return *this = *this >> shift;
                }
                constexpr operator uint64_t() const
                {
                    return low;
                }
                std::string toStringBase10() const
                {
                    if (0 == high)
                    {
                        return std::to_string(low);
                    }
                    constexpr uint64_t BASE(10000'0000'0000'0000);
                    Uint128 copy(*this);
                    std::string s;
                    s = ui64to_string_base10(uint64_t(copy.selfDivRem(BASE)), 16) + s;
                    s = ui64to_string_base10(uint64_t(copy.selfDivRem(BASE)), 16) + s;
                    return std::to_string(uint64_t(copy.selfDivRem(BASE))) + s;
                }
                void printDec() const
                {
                    std::cout << std::dec << toStringBase10() << '\n';
                }
                void printHex() const
                {
                    std::cout << std::hex << "0x" << high << ' ' << low << std::dec << '\n';
                }
                constexpr uint64_t high64() const
                {
                    return high;
                }
                constexpr uint64_t low64() const
                {
                    return low;
                }

            private:
                uint64_t low, high;
            };

            class Uint192
            {
            public:
                using ConstRef = const Uint192 &;
                Uint192() = default;
                constexpr Uint192(uint64_t lo, uint64_t mi = 0, uint64_t hi = 0) : low(lo), mid(mi), high(hi) {}
                constexpr Uint192(Uint128 n) : low(n.low64()), mid(n.high64()), high(0) {}

                constexpr Uint192 &operator+=(ConstRef rhs)
                {
                    bool cf = false;
                    low = add_half(low, rhs.low, cf);
                    mid = add_carry(mid, rhs.mid, cf);
                    high = high + rhs.high + cf;
                    return *this;
                }
                constexpr Uint192 &operator-=(ConstRef rhs)
                {
                    bool bf = false;
                    low = sub_half(low, rhs.low, bf);
                    mid = sub_borrow(mid, rhs.mid, bf);
                    high = high - rhs.high - bf;
                    return *this;
                }
                constexpr Uint192 &operator/=(ConstRef rhs)
                {
                    return *this = *this / rhs;
                }
                constexpr Uint192 &operator%=(ConstRef rhs)
                {
                    return *this = *this % rhs;
                }
                friend constexpr Uint192 operator+(Uint192 lhs, ConstRef rhs)
                {
                    return lhs += rhs;
                }
                friend constexpr Uint192 operator-(Uint192 lhs, ConstRef rhs)
                {
                    return lhs += rhs;
                }
                constexpr Uint192 operator/(uint64_t rhs) const
                {
                    Uint192 result(*this);
                    result.selfDivRem(rhs);
                    return result;
                }
                constexpr Uint192 operator%(uint64_t rhs) const
                {
                    Uint192 result(*this);
                    return result.selfDivRem(rhs);
                }
                constexpr Uint192 subNorm(Uint192 mod) const
                {
                    bool bf = false;
                    mod.low = sub_half(low, mod.low, bf);
                    mod.mid = sub_borrow(mid, mod.mid, bf);
                    mod.high = high - mod.high - bf;
                    // mask = 0b111...1 if *this < mod
                    // res = mod if *this >= mod
                    auto mask = uint64_t(0) - uint64_t(mod.high > high), nmask = ~mask;
                    mod.low = (mod.low & nmask) | (low & mask);
                    mod.mid = (mod.mid & nmask) | (mid & mask);
                    mod.high = (mod.high & nmask) | (high & mask);
                    return mod;
                }
                constexpr Uint192 operator<<(int shift) const
                {
                    if (shift == 0)
                    {
                        return *this;
                    }
                    shift %= 192;
                    shift = shift < 0 ? shift + 192 : shift;
                    if (shift < 64)
                    {
                        return Uint192(low << shift, (mid << shift) | (low >> (64 - shift)), (high << shift) | (mid >> (64 - shift)));
                    }
                    else if (shift < 128)
                    {
                        shift -= 64;
                        return Uint192(0, low << shift, (mid << shift) | (low >> (64 - shift)));
                    }
                    return Uint192(0, 0, low << (shift - 128));
                }
                friend constexpr bool operator<(Uint192 lhs, ConstRef rhs)
                {
                    bool bf = false;
                    sub_half(lhs.low, rhs.low, bf);
                    sub_borrow(lhs.mid, rhs.mid, bf);
                    uint64_t high = lhs.high - rhs.high - bf;
                    return high > lhs.high;
                }
                friend constexpr bool operator>(ConstRef lhs, ConstRef rhs)
                {
                    return rhs < lhs;
                }
                friend constexpr bool operator<=(ConstRef lhs, ConstRef rhs)
                {
                    return !(lhs > rhs);
                }
                friend constexpr bool operator>=(ConstRef lhs, ConstRef rhs)
                {
                    return !(lhs < rhs);
                }
                friend constexpr bool operator==(ConstRef lhs, ConstRef rhs)
                {
                    return lhs.high == rhs.high && lhs.mid == rhs.mid && lhs.low == rhs.low;
                }
                friend constexpr bool operator!=(ConstRef lhs, ConstRef rhs)
                {
                    return !(lhs == rhs);
                }
                static Uint192 mul128x64(Uint128 a, uint64_t b)
                {
                    Uint192 result;
                    uint64_t lo, hi;
                    mul64x64to128(b, a.low64(), result.low, result.mid);
                    mul64x64to128(b, a.high64(), lo, result.high);
                    result.mid += lo;
                    result.high += (result.mid < lo);
                    return result;
                }
                static constexpr Uint192 mul64x64x64(uint64_t a, uint64_t b, uint64_t c)
                {
                    Uint192 result{};
                    uint64_t lo{}, hi{};
                    mul64x64to128_base(a, b, lo, hi);
                    mul64x64to128_base(c, lo, result.low, result.mid);
                    mul64x64to128_base(c, hi, lo, result.high);
                    result.mid += lo;
                    result.high += (result.mid < lo);
                    return result;
                }
                constexpr uint64_t selfDivRem(uint64_t divisor)
                {
                    uint64_t divid1 = high % divisor, divid0 = mid;
                    high /= divisor;
                    mid = div128by64to64(divid1, divid0, divisor);
                    divid1 = divid0, divid0 = low;
                    low = div128by64to64(divid1, divid0, divisor);
                    return divid0;
                }
                constexpr Uint192 rShift64() const
                {
                    return Uint192(mid, high, 0);
                }
                constexpr operator uint64_t() const
                {
                    return low;
                }
                constexpr uint64_t high64() const
                {
                    return high;
                }
                constexpr uint64_t mid64() const
                {
                    return mid;
                }
                constexpr uint64_t low64() const
                {
                    return low;
                }
                std::string toStringBase10() const
                {
                    if (high == 0)
                    {
                        return Uint128(mid, low).toStringBase10();
                    }
                    constexpr uint64_t BASE(10000'0000'0000'0000);
                    Uint192 copy(*this);
                    std::string s;
                    s = ui64to_string_base10(uint64_t(copy.selfDivRem(BASE)), 16) + s;
                    s = ui64to_string_base10(uint64_t(copy.selfDivRem(BASE)), 16) + s;
                    s = ui64to_string_base10(uint64_t(copy.selfDivRem(BASE)), 16) + s;
                    return std::to_string(uint64_t(copy.selfDivRem(BASE))) + s;
                }
                void printDec() const
                {
                    std::cout << std::dec << toStringBase10() << '\n';
                }
                void printHex() const
                {
                    std::cout << std::hex << "0x" << high << ' ' << mid << ' ' << low << std::dec << '\n';
                }

            private:
                uint64_t low, mid, high;
            };

            template <typename UInt128Ty>
            constexpr uint64_t high64(const UInt128Ty &n)
            {
                return n >> 64;
            }
            constexpr uint64_t high64(const Uint128 &n)
            {
                return n.high64();
            }

#ifdef HINT_INT128
            using Uint128Defaulf = __uint128_t;
#else
            using Uint128Defaulf = Uint128;
#endif // UINT128T
        }

        template <typename NumTy>
        class DivExecutor
        {
        public:
            constexpr DivExecutor(NumTy divisor_in) : divisor(divisor_in)
            {
                assert(divisor > 1);
                inv = getInv(divisor, shift);
                divisor <<= shift;
            }
            // Return dividend / divisor, dividend %= divisor
            NumTy divRem(NumTy dividend_hi, NumTy dividend_lo, NumTy &rem) const
            {
                if (shift > 0)
                {
                    dividend_hi = (dividend_hi << shift) | dividend_lo >> (NUM_BITS - shift);
                    dividend_lo <<= shift;
                }
                NumTy quot = divRemNorm(dividend_hi, dividend_lo, rem);
                rem >>= shift;
                return quot;
            }
            NumTy prodDivRem(NumTy a, NumTy b, NumTy &rem) const
            {
                mul_binary(a << shift, b, a, b);
                NumTy quot = divRemNorm(b, a, rem);
                rem >>= shift;
                return quot;
            }
            // Reference:N. Möller and T. Granlund, "Improved Division by Invariant Integers,"
            // in IEEE Transactions on Computers, vol. 60, no. 2, pp. 165-175, Feb. 2011, doi: 10.1109/TC.2010.143.
            // Best performance, optimized first branch.
            NumTy divRemNorm(NumTy dividend_hi, NumTy dividend_lo, NumTy &rem) const
            {
                NumTy lo, hi, quot, mask;
                mul_binary(dividend_hi, inv, lo, hi);
                lo += dividend_lo;
                hi += dividend_hi + (lo < dividend_lo);
                quot = hi + 1;
                rem = dividend_lo - quot * divisor;
                mask = -NumTy(rem > lo); // mask = -1 if rem > dividend
                rem += divisor & mask;   // rem += divisor if rem > dividend
                quot += mask;            // quot -= 1 if rem > dividend
                if (rem >= divisor)
                {
                    rem -= divisor;
                    quot++;
                }
                return quot;
            }

        private:
            static constexpr NumTy getInv(NumTy divisor, int &leading_zero)
            {
                constexpr NumTy MAX = hint::all_one<NumTy>(NUM_BITS);
                leading_zero = hint::hint_clz(divisor);
                divisor <<= leading_zero;
                NumTy rem{};
                return div_dword_word(MAX - divisor, MAX, divisor, rem);
            }

            NumTy divisor = 0;
            NumTy inv = 0;
            int shift = 0;
            static constexpr int NUM_BITS = sizeof(NumTy) * CHAR_BIT;
        };
        template <typename NumTy>
        constexpr int DivExecutor<NumTy>::NUM_BITS;

        template <typename NumTy>
        class BaseExecutor
        {
        public:
            static constexpr int NUM_BITS = sizeof(NumTy) * CHAR_BIT;
            static_assert(std::is_unsigned<NumTy>::value, "BaseExecutor only supports unsigned integer types");
            // using ProdTy = typename UintType<NUM_BITS>::NextType::Type;
            // using SignTy = typename UintType<NUM_BITS>::SignType;
            constexpr BaseExecutor(NumTy base_in) : base(base_in), div_exe(base_in)
            {
                assert(2 <= base);
                assert(base <= hint::all_one<NumTy>(NUM_BITS));
            }

            constexpr NumTy halfBase() const
            {
                return (base + 1) / 2;
            }

            constexpr bool checkDivisor(NumTy divisor) const
            {
                return divisor < base;
            }

            constexpr NumTy addCarry(NumTy a, NumTy b, bool &cf) const
            {
                return addHalf(a, b + cf, cf);
            }
            constexpr NumTy addHalf(NumTy a, NumTy b, bool &cf) const
            {
                if (a >= base - b)
                {
                    a -= base;
                    cf = true;
                }
                else
                {
                    cf = false;
                }
                return a + b;
            }
            // Return (a + cf) % base
            constexpr NumTy subBorrow(NumTy a, NumTy b, bool &bf) const
            {
                return subHalf(a, b + bf, bf);
            }
            constexpr NumTy subHalf(NumTy a, NumTy b, bool &bf) const
            {
                if (a < b)
                {
                    a += base;
                    bf = true;
                }
                else
                {
                    bf = false;
                }
                return a - b;
            }
            constexpr NumTy addCf(NumTy a, bool &cf) const
            {
                a += cf;
                if (a >= base)
                {
                    a -= base;
                    cf = true;
                }
                return a;
            }
            // Return (a - cf) % base
            constexpr NumTy subBf(NumTy a, bool &bf) const
            {
                auto flag = bf;
                if (a < flag)
                {
                    a += base;
                    bf = true;
                }
                return a - flag;
            }
            constexpr uint64_t divRemBase(extend_int::Uint192 &dividend) const
            {
                uint64_t rem = 0, hi = 0, mi = 0, lo = 0;
                hi = div_exe.divRem(rem, dividend.high64(), rem);
                mi = div_exe.divRem(rem, dividend.mid64(), rem);
                lo = div_exe.divRem(rem, dividend.low64(), rem);
                dividend = extend_int::Uint192{lo, mi, hi};
                return rem;
            }
            constexpr void add_carry_x4(const NumTy in1[4], const NumTy in2[4], NumTy sum[4], bool &carry) const
            {
                sum[0] = addCarry(in1[0], in2[0], carry);
                sum[1] = addCarry(in1[1], in2[1], carry);
                sum[2] = addCarry(in1[2], in2[2], carry);
                sum[3] = addCarry(in1[3], in2[3], carry);
            }
            constexpr void sub_borrow_x4(const NumTy in1[4], const NumTy in2[4], NumTy sum[4], bool &borrow) const
            {
                sum[0] = subBorrow(in1[0], in2[0], borrow);
                sum[1] = subBorrow(in1[1], in2[1], borrow);
                sum[2] = subBorrow(in1[2], in2[2], borrow);
                sum[3] = subBorrow(in1[3], in2[3], borrow);
            }
            constexpr void addCfX4(const NumTy in[4], NumTy sum[4], bool &carry) const
            {
                sum[0] = addCf(in[0], carry);
                sum[1] = addCf(in[1], carry);
                sum[2] = addCf(in[2], carry);
                sum[3] = addCf(in[3], carry);
            }
            constexpr void subBfX4(const NumTy in[4], NumTy sum[4], bool &borrow) const
            {
                sum[0] = subBf(in[0], borrow);
                sum[1] = subBf(in[1], borrow);
                sum[2] = subBf(in[2], borrow);
                sum[3] = subBf(in[3], borrow);
            }

            // a * b = lo + hi * base
            void mulInBase(NumTy a, NumTy b, NumTy &lo, NumTy &hi) const
            {
                hi = div_exe.prodDivRem(a, b, lo);
            }

            // lo + hi * base = [lo, hi]
            void dualBaseToBin(NumTy hi, NumTy lo, NumTy &lo_bin, NumTy &hi_bin) const
            {
                mul_binary(hi, base, lo_bin, hi_bin);
                lo_bin += lo;
                hi_bin += (lo_bin < lo);
            }

            // a * b + c = lo + hi * base
            void mulAdd(NumTy a, NumTy b, NumTy c, NumTy &lo, NumTy &hi) const
            {
                mulInBase(a, b, lo, hi);
                bool cf;
                lo = addHalf(lo, c, cf);
                hi += cf;
            }

            NumTy mulAddX4(const NumTy in[4], NumTy num_mul, NumTy num_add, NumTy out[4]) const
            {
                mulAdd(in[0], num_mul, num_add, out[0], num_add);
                mulAdd(in[1], num_mul, num_add, out[1], num_add);
                mulAdd(in[2], num_mul, num_add, out[2], num_add);
                mulAdd(in[3], num_mul, num_add, out[3], num_add);
                return num_add;
            }

            void mulAddPlusSelf(NumTy in, NumTy num_mul, NumTy &out, NumTy &carry) const
            {
                bool cf;
                NumTy hi, lo;
                mulInBase(in, num_mul, lo, hi);
                lo = addHalf(lo, out, cf);
                hi += cf;
                out = addHalf(lo, carry, cf);
                carry = hi + cf;
            }
            NumTy addEqMulAddX8(const NumTy in[], size_t len, NumTy out[], NumTy num_mul, NumTy num_add) const
            {
                len -= len % 8;
                for (size_t i = 0; i < len; i += 8)
                {
                    mulAddPlusSelf(in[i], num_mul, out[i], num_add);
                    mulAddPlusSelf(in[i + 1], num_mul, out[i + 1], num_add);
                    mulAddPlusSelf(in[i + 2], num_mul, out[i + 2], num_add);
                    mulAddPlusSelf(in[i + 3], num_mul, out[i + 3], num_add);
                    mulAddPlusSelf(in[i + 4], num_mul, out[i + 4], num_add);
                    mulAddPlusSelf(in[i + 5], num_mul, out[i + 5], num_add);
                    mulAddPlusSelf(in[i + 6], num_mul, out[i + 6], num_add);
                    mulAddPlusSelf(in[i + 7], num_mul, out[i + 7], num_add);
                }
                return num_add;
            }

            // NumTy divRemBase(ProdTy &dividend) const
            // {
            //     return div_exe.divRem(dividend);
            // }

        private:
            using DivExecutorTy = DivExecutor<NumTy>;
            NumTy base;
            DivExecutorTy div_exe;
        };
        template <typename NumTy>
        constexpr int BaseExecutor<NumTy>::NUM_BITS;

        template <typename UintTy>
        class BaseExecutorBinary
        {
        public:
            static_assert(std::is_integral<UintTy>::value, "UintTy must be integral type");
            static_assert(std::is_unsigned<UintTy>::value, "UintTy must be unsigned type");
            static constexpr int NUM_BITS = std::numeric_limits<UintTy>::digits;
            static constexpr UintTy MAX_NUM = all_one<UintTy>(NUM_BITS);
            constexpr BaseExecutorBinary() = default;

            constexpr BaseExecutorBinary(UintTy) { /*NOP*/ }

            static constexpr UintTy halfBase()
            {
                return UintTy(1) << (NUM_BITS - 1);
            }

            template <typename T>
            static constexpr bool checkDivisor(T divisor)
            {
                return divisor <= MAX_NUM;
            }

            static constexpr UintTy addCarry(UintTy a, UintTy b, bool &cf)
            {
                return add_carry(a, b, cf);
            }
            static constexpr UintTy addHalf(UintTy a, UintTy b, bool &cf)
            {
                return add_half(a, b, cf);
            }
            static constexpr UintTy subBorrow(UintTy a, UintTy b, bool &bf)
            {
                return sub_borrow(a, b, bf);
            }
            static constexpr UintTy subHalf(UintTy a, UintTy b, bool &bf)
            {
                return sub_half(a, b, bf);
            }
            static constexpr UintTy addCf(UintTy a, bool &cf)
            {
                a += cf;
                cf = a < cf;
                return a;
            }
            static constexpr UintTy subBf(UintTy a, bool &bf)
            {
                UintTy diff = a - bf;
                bf = diff > a;
                return diff;
            }

            static constexpr uint64_t divRemBase(extend_int::Uint192 &dividend)
            {
                uint64_t rem = dividend;
                dividend = dividend.rShift64();
                return rem;
            }
            static void addCarryX4(const UintTy in1[4], const UintTy in2[4], UintTy sum[4], bool &carry) noexcept
            {
                add_carry_x4(in1, in2, sum, carry);
            }
            static void subBorrowX4(const UintTy in1[4], const UintTy in2[4], UintTy diff[4], bool &borrow) noexcept
            {
                sub_borrow_x4(in1, in2, diff, borrow);
            }
            static void addCfX4(const UintTy in[4], UintTy sum[4], bool &carry) noexcept
            {
                add_cf_x4(in, sum, carry);
            }
            static void subBfX4(const UintTy in[4], UintTy diff[4], bool &borrow) noexcept
            {
                sub_bf_x4(in, diff, borrow);
            }

            static void mulInBase(UintTy a, UintTy b, UintTy &lo, UintTy &hi) noexcept
            {
                mul_binary(a, b, lo, hi);
            }

            static void dualBaseToBin(UintTy hi, UintTy lo, UintTy &lo_bin, UintTy &hi_bin) noexcept
            {
                lo_bin = lo;
                hi_bin = hi;
            }

            static void mulAdd(UintTy a, UintTy b, UintTy c, UintTy &lo, UintTy &hi) noexcept
            {
                mul_add(a, b, c, lo, hi);
            }

            static UintTy mulAddX4(const UintTy in[4], UintTy num_mul, UintTy num_add, UintTy out[4]) noexcept
            {
                return mul_add_x4(in, num_mul, num_add, out);
            }

            static void mulAddPlusSelf(UintTy in, UintTy num_mul, UintTy &out, UintTy &carry) noexcept
            {
                mul_add_plus_self(in, num_mul, out, carry);
            }
            static UintTy addEqMulAddX8(const UintTy in[], size_t len, UintTy out[], UintTy num_mul, UintTy num_add) noexcept
            {
                return add_eq_mul_add_x8(in, len, out, num_mul, num_add);
            }
        };
        template <uint64_t MOD1, uint64_t MOD2, uint64_t MOD3>
        struct CRTMod3T
        {
            using Uint128 = extend_int::Uint128;
            using Uint192 = extend_int::Uint192;

            static Uint192 crt3(uint64_t n1, uint64_t n2, uint64_t n3) noexcept
            {
                constexpr auto prod = Uint192::mul64x64x64(MOD1, MOD2, MOD3);
                constexpr auto mod12 = Uint128::mul64x64(MOD1, MOD2);
                constexpr auto mod13 = Uint128::mul64x64(MOD1, MOD3);
                constexpr auto mod23 = Uint128::mul64x64(MOD2, MOD3);

                Uint192 result = Uint192::mul128x64(mod23, n1);
                result += Uint192::mul128x64(mod13, n2);
                result += Uint192::mul128x64(mod12, n3);
                result = result.subNorm(prod);
                return result.subNorm(prod);
            }
        };

        // remove leading zeros, return the true length
        template <typename T>
        constexpr size_t count_ture_length(const T array[], size_t length)
        {
            if (nullptr == array)
            {
                return 0;
            }
            while (length > 0 && array[length - 1] == 0)
            {
                length--;
            }
            return length;
        }
        template <typename T>
        constexpr size_t count_trailing_zero(const T array[], size_t length)
        {
            if (nullptr == array)
            {
                return 0;
            }
            auto it = array, end = array + length;
            while (it < end && *it == 0)
            {
                it++;
            }
            return it - array;
        }

        template <typename T>
        constexpr void remove_leading_zeros(const T array[], size_t &length)
        {
            length = count_ture_length(array, length);
        }

        // Absolute compare, return 1 if a > b, -1 if a < b, 0 if a == b
        // Return the diffent length if a != b
        template <typename T>
        [[nodiscard]] constexpr auto abs_compare_with_length(const T in1[], size_t len1, const T in2[], size_t len2)
        {
            struct CompareResult
            {
                size_t diff_len = 0;
                int cmp = 0;
            };
            CompareResult result;
            if (len1 != len2)
            {
                result.diff_len = std::max(len1, len2);
                result.cmp = len1 > len2 ? 1 : -1;
                return result;
            }
            size_t i = len1;
            while (i > 0)
            {
                i--;
                if (in1[i] != in2[i])
                {
                    result.diff_len = i + 1;
                    result.cmp = in1[i] > in2[i] ? 1 : -1;
                    break;
                }
            }
            return result;
        }

        // Absolute compare, return 1 if a > b, -1 if a < b, 0 if a == b
        template <typename T>
        [[nodiscard]] constexpr int abs_compare(const T in1[], size_t len1, const T in2[], size_t len2)
        {
            return abs_compare_with_length(in1, len1, in2, len2).cmp;
        }
    }
    namespace transform
    {
        template <typename T>
        inline void transform2(T &sum, T &diff)
        {
            T temp0 = sum, temp1 = diff;
            sum = temp0 + temp1;
            diff = temp0 - temp1;
        }

        // 二进制逆序
        template <typename It>
        void binary_reverse_swap(It begin, It end)
        {
            const size_t len = end - begin;
            assert(hint::is_2pow(len));
            // 左下标小于右下标时交换,防止重复交换
            auto smaller_swap = [=](It it_left, It it_right)
            {
                if (it_left < it_right)
                {
                    std::swap(it_left[0], it_right[0]);
                }
            };
            // 若i的逆序数的迭代器为last,则返回i+1的逆序数的迭代器
            auto get_next_bitrev = [=](It last)
            {
                size_t k = len / 2, indx = last - begin;
                indx ^= k;
                while (k > indx)
                {
                    k >>= 1;
                    indx ^= k;
                };
                return begin + indx;
            };
            // 长度较短的普通逆序
            if (len <= 16)
            {
                for (auto i = begin + 1, j = begin + len / 2; i < end - 1; i++)
                {
                    smaller_swap(i, j);
                    j = get_next_bitrev(j);
                }
                return;
            }
            const size_t len_8 = len / 8;
            const auto last = begin + len_8;
            auto i0 = begin + 1, i1 = i0 + len / 2, i2 = i0 + len / 4, i3 = i1 + len / 4;
            for (auto j = begin + len / 2; i0 < last; i0++, i1++, i2++, i3++)
            {
                smaller_swap(i0, j);
                smaller_swap(i1, j + 1);
                smaller_swap(i2, j + 2);
                smaller_swap(i3, j + 3);
                smaller_swap(i0 + len_8, j + 4);
                smaller_swap(i1 + len_8, j + 5);
                smaller_swap(i2 + len_8, j + 6);
                smaller_swap(i3 + len_8, j + 7);
                j = get_next_bitrev(j);
            }
        }

        // 二进制逆序
        template <typename T>
        void binary_reverse_swap(T ary, const size_t len)
        {
            binary_reverse_swap(ary, ary + len);
        }
        // 自动类型，自检查快速数论变换
        namespace ntt
        {
            constexpr uint64_t MOD0 = 2485986994308513793, ROOT0 = 5; // 69 * 2^55 + 1
            constexpr uint64_t MOD1 = 1945555039024054273, ROOT1 = 5; // 27 * 2^56 + 1
            constexpr uint64_t MOD2 = 4179340454199820289, ROOT2 = 3; // 29 * 2^57 + 1
            constexpr uint32_t MOD3 = 998244353, ROOT3 = 3;           // 119 * 2^23 + 1
            constexpr uint32_t MOD4 = 754974721, ROOT4 = 11;          // 45 * 2^24 + 1
            constexpr uint32_t MOD5 = 469762049, ROOT5 = 3;           // 7 * 2^26 + 1

            //  Montgomery ModInt
            template <uint64_t MOD>
            class MontIntLazy
            {
            public:
                static constexpr int MOD_BITS = hint_log2(MOD) + 1;
                static_assert(MOD_BITS <= 30 || 32 < MOD_BITS && MOD_BITS <= 62, "MOD_BITS not in range [0, 30] or [33, 62]");

                using Uint128Fast = utility::extend_int::Uint128Defaulf;
                using Uint128 = utility::extend_int::Uint128;
                using IntType = typename std::conditional<MOD_BITS <= 30, uint32_t, uint64_t>::type;
                using ProdTypeFast = typename std::conditional<MOD_BITS <= 30, uint64_t, Uint128Fast>::type;
                using ProdType = typename std::conditional<MOD_BITS <= 30, uint64_t, Uint128>::type;
                static constexpr int R_BITS = sizeof(IntType) * CHAR_BIT;

                static constexpr IntType getR()
                {
                    constexpr IntType HALF = (IntType(1) << (R_BITS - 1)) % mod();
                    return HALF * 2 % mod();
                }
                static constexpr IntType R = getR();                               // R % MOD
                static constexpr IntType R2 = utility::mulMod(R, R, IntType(MOD)); // R^2 % MOD
                static constexpr IntType MOD_INV = inv_mod2pow(MOD, R_BITS);       // MOD^-1 % R
                static constexpr IntType MOD_INV_NEG = IntType(0) - MOD_INV;       // -MOD^-1 % R
                static constexpr IntType MOD2 = MOD * 2;                           // MOD * 2
                static_assert(IntType(MOD * MOD_INV) == 1, "MOD_INV not correct");

                constexpr MontIntLazy() = default;
                constexpr MontIntLazy(IntType n) : data(toMont(n)) {}

                constexpr IntType raw() const
                {
                    return data;
                }
                constexpr MontIntLazy &operator+=(MontIntLazy rhs)
                {
                    data += rhs.data;
                    data = norm2(data);
                    return *this;
                }
                constexpr MontIntLazy &operator-=(MontIntLazy rhs)
                {
                    const IntType mask = IntType(0) - (data < rhs.data);
                    data = data - rhs.data + (MOD2 & mask);
                    return *this;
                }
                constexpr MontIntLazy &operator*=(MontIntLazy rhs)
                {
                    data = redc(ProdType(data) * rhs.data);
                    return *this;
                }

                friend constexpr MontIntLazy operator+(MontIntLazy lhs, MontIntLazy rhs)
                {
                    return lhs += rhs;
                }
                friend constexpr MontIntLazy operator-(MontIntLazy lhs, MontIntLazy rhs)
                {
                    return lhs -= rhs;
                }
                friend MontIntLazy operator*(MontIntLazy lhs, MontIntLazy rhs)
                {
                    ProdTypeFast prod = ProdTypeFast(lhs.data) * rhs.data;
                    lhs.data = redcLazy(prod);
                    return lhs;
                }
                constexpr MontIntLazy operator-() const
                {
                    MontIntLazy res = this->norm1();
                    res.data = mod() - data;
                    return res;
                }
                constexpr MontIntLazy norm1() const
                {
                    MontIntLazy res{};
                    res.data = norm1(data);
                    return res;
                }
                constexpr MontIntLazy norm2() const
                {
                    MontIntLazy res{};
                    res.data = norm2(data);
                    return res;
                }
                template <IntType N = 1>
                constexpr MontIntLazy norm() const
                {
                    MontIntLazy res{};
                    res.data = norm<N>(data);
                    return res;
                }
                static constexpr IntType norm1(IntType n)
                {
                    const IntType mask = IntType(0) - (n >= mod());
                    return n - (MOD & mask);
                }
                static constexpr IntType norm2(IntType n)
                {
                    constexpr IntType MOD2 = mod<2>();
                    const IntType mask = IntType(0) - (n >= MOD2);
                    return n - (MOD2 & mask);
                }
                template <IntType N = 1>
                static constexpr IntType norm(IntType n)
                {
                    constexpr IntType MOD_N = mod<N>();
                    const IntType mask = IntType(0) - (n >= MOD_N);
                    return n - (MOD_N & mask);
                }
                constexpr MontIntLazy add(MontIntLazy rhs) const
                {
                    rhs.data = data + rhs.data;
                    return rhs;
                }
                constexpr MontIntLazy sub(MontIntLazy rhs) const
                {
                    rhs.data = data - rhs.data + MOD2;
                    return rhs;
                }
                constexpr operator IntType() const
                {
                    return toInt(data);
                }
                template <IntType N = 1>
                static constexpr IntType mod()
                {
                    constexpr IntType MOD_N = MOD * N;
                    return MOD_N;
                }
                static constexpr MontIntLazy montR()
                {
                    constexpr MontIntLazy res(R);
                    return res;
                }
                // R / 32 <= MOD
                constexpr MontIntLazy shrinkToMod4X(std::true_type) const
                {
                    MontIntLazy res{};
                    res.data = this->data % MOD;
                    return res;
                }
                // R / 16 <= MOD < R / 8
                constexpr MontIntLazy shrinkToMod4X(std::false_type) const
                {
                    return this->norm<8>().template norm<4>();
                }
                // R / 16 <= MOD < R / 8
                constexpr MontIntLazy shrinkToMod4(std::true_type) const
                {
                    constexpr IntType R_16 = IntType(1) << (R_BITS - 4);
                    using SMALL_MOD = std::integral_constant<bool, (MOD < R_16)>;
                    return shrinkToMod4X(SMALL_MOD{});
                }
                // R / 8 <= MOD < R / 4
                constexpr MontIntLazy shrinkToMod4(std::false_type) const
                {
                    return this->norm<4>();
                }
                // R <= MOD * 8 , R > MOD * 4
                constexpr MontIntLazy shrinkToMod4() const
                {
                    constexpr IntType R_8 = IntType(1) << (R_BITS - 3);
                    using SMALL_MOD = std::integral_constant<bool, (MOD < R_8)>;
                    return shrinkToMod4(SMALL_MOD{});
                }

                static constexpr uint32_t toMont(uint32_t n)
                {
                    return redc(uint64_t(n) * R2);
                }
                static constexpr uint32_t toInt(uint32_t n)
                {
                    return redc(uint64_t(n));
                }
                static constexpr uint64_t toMont(uint64_t n)
                {
                    uint64_t lo{}, hi{};
                    utility::mul64x64to128_base(n, R2, lo, hi);
                    return redc(Uint128{lo, hi});
                }
                static constexpr uint64_t toInt(uint64_t n)
                {
                    return redc(Uint128{n});
                }

                static constexpr uint32_t redcLazy(uint64_t n)
                {
                    uint32_t prod = uint32_t(n) * MOD_INV_NEG;
                    return (uint64_t(prod) * mod() + n) >> 32;
                }
                static constexpr uint32_t redc(uint64_t n)
                {
                    uint32_t res = redcLazy(n);
                    return norm1(res);
                }
                static uint64_t redcLazy(Uint128Fast n)
                {
                    uint64_t prod1 = uint64_t(n) * MOD_INV_NEG;
                    auto prod2 = Uint128Fast(prod1) * mod() + n;
                    return utility::extend_int::high64(prod2);
                }
                static constexpr uint64_t redc(Uint128 n)
                {
                    uint64_t prod1 = uint64_t(n) * MOD_INV_NEG, lo{}, hi{};
                    utility::mul64x64to128_base(prod1, mod(), lo, hi);
                    Uint128 prod2{lo, hi};
                    prod2 += n;
                    return norm<1>(prod2.high64());
                }

                static constexpr MontIntLazy one()
                {
                    constexpr MontIntLazy res(1);
                    return res;
                }
                static constexpr MontIntLazy negOne()
                {
                    constexpr MontIntLazy res(mod() - 1);
                    return res;
                }
                constexpr MontIntLazy inv() const
                {
                    return hint::qpow(*this, mod() - 2);
                }

            private:
                IntType data;
            };

            template <uint64_t MOD>
            constexpr int MontIntLazy<MOD>::MOD_BITS;
            template <uint64_t MOD>
            constexpr int MontIntLazy<MOD>::R_BITS;
            template <uint64_t MOD>
            constexpr typename MontIntLazy<MOD>::IntType MontIntLazy<MOD>::R;
            template <uint64_t MOD>
            constexpr typename MontIntLazy<MOD>::IntType MontIntLazy<MOD>::R2;
            template <uint64_t MOD>
            constexpr typename MontIntLazy<MOD>::IntType MontIntLazy<MOD>::MOD_INV;
            template <uint64_t MOD>
            constexpr typename MontIntLazy<MOD>::IntType MontIntLazy<MOD>::MOD_INV_NEG;
            template <uint64_t MOD>
            constexpr typename MontIntLazy<MOD>::IntType MontIntLazy<MOD>::MOD2;

            // in: in_out0<4p, in_ou1<4p
            // out: in_out0<4p, in_ou1<4p
            template <typename ModIntType>
            inline void dit_butterfly2(ModIntType &in_out0, ModIntType &in_out1, const ModIntType &omega)
            {
                auto x = in_out0.norm2();
                auto y = in_out1 * omega;
                in_out0 = x.add(y);
                in_out1 = x.sub(y);
            }

            // in: in_out0<2p, in_ou1<2p
            // out: in_out0<2p, in_ou1<2p
            template <typename ModIntType>
            inline void dif_butterfly2(ModIntType &in_out0, ModIntType &in_out1, const ModIntType &omega)
            {
                auto x = in_out0 + in_out1;
                auto y = in_out0.sub(in_out1);
                in_out0 = x;
                in_out1 = y * omega;
            }
            template <typename ModInt>
            class BinRevTable
            {
            public:
                static constexpr int MAX_LOG_LEN = CHAR_BIT * sizeof(ModInt);
                static constexpr size_t MAX_LEN = size_t(1) << MAX_LOG_LEN;

                using IntType = std::conditional_t<sizeof(ModInt) <= 4, uint32_t, uint64_t>;

                BinRevTable(IntType root_in, size_t factor = 1, size_t div = 1) : root(root_in)
                {
                    constexpr int LOG_LEN = hint::hint_ctz(ModInt::mod() - 1);
                    assert(hint::is_2pow(div));
                    table[0] = getOmega(2 * div, factor, false);
                    for (int i = 1; i <= LOG_LEN; i++)
                    {
                        const size_t rev_indx = 1;
                        const size_t last_indx = ((size_t(1) << i) - 1) << 1;
                        table[i] = getOmega((size_t(1) << (i + 1)) * div, (last_indx - rev_indx) * factor, true);
                    }
                }

                ModInt evenToOdd(ModInt last_even) const
                {
                    last_even = last_even * table[0];
                    return last_even.norm1();
                }

                ModInt omega1() const
                {
                    return table[0];
                }

                ModInt getNext(ModInt last, IntType indx) const
                {
                    if (0 == indx)
                    {
                        return ModInt::one();
                    }
                    indx = indx % 2 ? 0 : hint_ctz(indx);
                    last = last * table[indx];
                    return last.norm1();
                }

                ModInt getOmega(size_t n, size_t index, bool conj = false) const
                {
                    index = (ModInt::mod() - 1) / n * index;
                    ModInt omega = qpow(root, index);
                    return conj ? omega.inv() : omega;
                }

            private:
                ModInt table[MAX_LOG_LEN];
                ModInt root;
            };

            template <uint64_t MOD, uint64_t ROOT>
            struct NTT
            {
                static constexpr uint64_t mod()
                {
                    return MOD;
                }
                static constexpr uint64_t root()
                {
                    return ROOT;
                }
                static constexpr uint64_t rootInv()
                {
                    constexpr uint64_t IROOT = mod_inv<int64_t>(ROOT, MOD);
                    return IROOT;
                }

                static_assert(root() < mod(), "ROOT must be smaller than MOD");
                static constexpr int MOD_BITS = hint_log2(mod()) + 1;
                static constexpr int MAX_LOG_LEN = hint_ctz(mod() - 1);

                static constexpr size_t getMaxLen()
                {
                    if (MAX_LOG_LEN < sizeof(size_t) * CHAR_BIT)
                    {
                        return size_t(1) << MAX_LOG_LEN;
                    }
                    return size_t(1) << (sizeof(size_t) * CHAR_BIT - 1);
                }
                static constexpr size_t NTT_MAX_LEN = getMaxLen();

                using ModInt = MontIntLazy<MOD>;
                using IntType = typename ModInt::IntType;
                using TableType = BinRevTable<ModInt>;
                static_assert(std::is_trivial<ModInt>::value, "ModInt must be trivial");
                static const TableType table, itable;

                static constexpr size_t L1_BYTE = size_t(1) << 15;
                static constexpr size_t LONG_BYTE = size_t(1) << 29;
                static constexpr size_t ITER_THRESHOLD = L1_BYTE / sizeof(ModInt);
                static constexpr size_t LONG_THRESHOLD = LONG_BYTE / sizeof(ModInt);
                static constexpr size_t SHORT_THRESHOLD = 15;

                static size_t findFitLen(size_t conv_len) noexcept
                {
                    size_t result = findFitLen(conv_len, SHORT_THRESHOLD);
                    // 长度为2的幂次且过长时将造成L1缓存颠簸，需要将长度转为非2的幂次
                    constexpr size_t DIV_LEN = int_floor2(SHORT_THRESHOLD), MUL_LEN = DIV_LEN + 1;
                    static_assert(MUL_LEN <= SHORT_THRESHOLD, "MUL_LEN can't be larger than SHORT_THRESHOLD");
                    if (result >= LONG_THRESHOLD && is_2pow(result))
                    {
                        result /= DIV_LEN;
                        result *= MUL_LEN;
                    }
                    return result;
                }
                static void convolution(IntType in_out[], IntType in[], size_t conv_len,
                                        ModInt weight = ModInt::one()) noexcept
                {
                    assert(checkConvLen(conv_len)); // check if conv_len is too long
                    auto p1 = reinterpret_cast<ModInt *>(in_out), p2 = reinterpret_cast<ModInt *>(in);
                    if (conv_len <= SHORT_THRESHOLD)
                    {
                        std::copy(in_out, in_out + conv_len, p1);
                        std::copy(in, in + conv_len, p2);
                        convolutionCyclicShort(p1, p2, conv_len);
                        for (size_t i = 0; i < conv_len; i++)
                        {
                            in_out[i] = (p1[i] * weight).norm1();
                        }
                        return;
                    }
                    ModInt buffer[64], ibuffer[64];
                    convolutionCyclic(in_out, in, conv_len, conv_len, conv_len, buffer, ibuffer, weight);
                }

                static void convolution(IntType in_out[], IntType in[], size_t len1, size_t len2, size_t conv_len,
                                        ModInt weight = ModInt::one()) noexcept
                {
                    assert(len1 + len2 - 1 <= conv_len);
                    assert(checkConvLen(conv_len)); // check if conv_len is too long
                    auto p1 = reinterpret_cast<ModInt *>(in_out), p2 = reinterpret_cast<ModInt *>(in);
                    if (conv_len <= ITER_THRESHOLD)
                    {
                        std::fill(in_out + len1, in_out + conv_len, IntType{});
                        std::fill(in + len2, in + conv_len, IntType{});
                    }
                    if (conv_len <= SHORT_THRESHOLD)
                    {
                        std::copy(in_out, in_out + conv_len, p1);
                        std::copy(in, in + conv_len, p2);
                        convolutionCyclicShort(p1, p2, conv_len);
                        for (size_t i = 0; i < conv_len; i++)
                        {
                            in_out[i] = (p1[i] * weight).norm1();
                        }
                        return;
                    }
                    ModInt buffer[64], ibuffer[64];
                    convolutionCyclic(in_out, in, len1, len2, conv_len, buffer, ibuffer, weight);
                }

            private:
                static size_t forwardCyclic(ModInt in_out[], size_t len, ModInt buffer[], bool assign_buf) noexcept
                {
                    auto p_buf = buffer;
                    for (size_t rank = len; rank > SHORT_THRESHOLD * 2; rank /= 2, p_buf++)
                    {
                        size_t stride = rank / 2, i = 1;
                        auto it0 = in_out, it1 = it0 + stride;
                        for (const auto end = it1; it0 < end; it0++, it1++)
                        {
                            auto t0 = it0[0].norm2(), t1 = it1[0].norm2();
                            it0[0] = t0.add(t1), it1[0] = t0.sub(t1);
                        }
                        ModInt omega = ModInt::one();
                        for (auto begin = in_out + rank; begin < in_out + len; begin += rank, i++)
                        {
                            omega = table.getNext(omega, i);
                            it0 = begin, it1 = it0 + stride;
                            for (const auto end = it1; it0 < end; it0++, it1++)
                            {
                                dit_butterfly2(it0[0], it1[0], omega);
                            }
                        }
                        if (assign_buf)
                        {
                            p_buf[0] = omega;
                        }
                    }
                    return p_buf - buffer;
                }
                static size_t forwardIter(ModInt in_out[], size_t len, ModInt buffer[], size_t idx, bool assign_buf) noexcept
                {
                    auto p_buf = buffer;
                    for (size_t rank = len; rank > SHORT_THRESHOLD * 2; rank /= 2, idx *= 2, p_buf++)
                    {
                        size_t stride = rank / 2, i = idx;
                        ModInt omega = p_buf[0];
                        for (auto begin = in_out; begin < in_out + len; begin += rank, i++)
                        {
                            auto it0 = begin, it1 = it0 + stride;
                            omega = table.getNext(omega, i);
                            for (const auto end = it1; it0 < end; it0++, it1++)
                            {
                                dit_butterfly2(it0[0], it1[0], omega);
                            }
                        }
                        if (assign_buf)
                        {
                            p_buf[0] = omega;
                        }
                    }
                    return p_buf - buffer;
                }
                static void backwardCyclic(ModInt in_out[], size_t len, size_t rank, ModInt ibuffer[]) noexcept
                {
                    for (; rank <= len; rank *= 2, ibuffer++)
                    {
                        size_t stride = rank / 2, i = 1;
                        auto it0 = in_out, it1 = it0 + stride, end = it1;
                        for (; it0 < end; it0++, it1++)
                        {
                            transform2(it0[0], it1[0]);
                        }
                        ModInt omega = ModInt::one();
                        for (auto begin = in_out + rank; begin < in_out + len; begin += rank, i++)
                        {
                            omega = itable.getNext(omega, i);
                            it0 = begin, it1 = it0 + stride;
                            for (const auto end = it1; it0 < end; it0++, it1++)
                            {
                                dif_butterfly2(it0[0], it1[0], omega);
                            }
                        }
                        ibuffer[0] = omega;
                    }
                }
                static void backwardIter(ModInt in_out[], size_t len, size_t rank, ModInt ibuffer[], size_t idx) noexcept
                {
                    for (; rank <= len; rank *= 2, idx /= 2, ibuffer++)
                    {
                        size_t stride = rank / 2, i = idx;
                        ModInt omega = ibuffer[0];
                        for (auto begin = in_out; begin < in_out + len; begin += rank, i++)
                        {
                            auto it0 = begin, it1 = it0 + stride;
                            omega = itable.getNext(omega, i);
                            for (const auto end = it1; it0 < end; it0++, it1++)
                            {
                                dif_butterfly2(it0[0], it1[0], omega);
                            }
                        }
                        ibuffer[0] = omega;
                    }
                }

                static void forwardIter(ModInt inout0[], ModInt inout1[], size_t rank, ModInt omega) noexcept
                {
                    auto it0 = inout0, it1 = it0 + rank, end = it1;
                    auto it2 = inout1, it3 = it2 + rank;
                    for (; it0 < end; it0++, it1++, it2++, it3++)
                    {
                        dit_butterfly2(it0[0], it1[0], omega);
                        dit_butterfly2(it2[0], it3[0], omega);
                    }
                }

                static size_t convolutionIterCyclic(ModInt in_out[], ModInt in[], size_t conv_len,
                                                    ModInt buffer[], ModInt ibuffer[]) noexcept
                {
                    if (in_out != in)
                    {
                        forwardCyclic(in_out, conv_len, buffer, false);
                    }
                    size_t layers = forwardCyclic(in, conv_len, buffer, true) + 1;
                    size_t rank = conv_len >> layers;
                    ModInt omega = buffer[layers];
                    auto it0 = in_out, it1 = in, end = it0 + conv_len;
                    for (size_t i = 0; it0 < end; it0 += rank * 2, it1 += rank * 2, i++)
                    {
                        omega = table.getNext(omega, i);
                        forwardIter(it0, it1, rank, omega);
                        convolutionShort(it0, it1, rank, omega);
                        convolutionShort(it0 + rank, it1 + rank, rank, -omega);
                    }
                    buffer[layers] = omega;
                    size_t factor = size_t(1) << layers;
                    backwardCyclic(in_out, conv_len, rank * 2, ibuffer);
                    return factor;
                }
                static void convolutionIter(ModInt in_out[], ModInt in[], size_t conv_len,
                                            ModInt buffer[], ModInt ibuffer[], size_t idx) noexcept
                {
                    if (in_out != in)
                    {
                        forwardIter(in_out, conv_len, buffer, idx, false);
                    }
                    size_t layers = forwardIter(in, conv_len, buffer, idx, true) + 1;
                    size_t rank = conv_len >> layers;
                    ModInt omega = buffer[layers];
                    idx <<= (layers - 1);
                    auto it0 = in_out, it1 = in, end = it0 + conv_len;
                    for (size_t i = idx; it0 < end; it0 += rank * 2, it1 += rank * 2, i++)
                    {
                        omega = table.getNext(omega, i);
                        forwardIter(it0, it1, rank, omega);
                        convolutionShort(it0, it1, rank, omega);
                        convolutionShort(it0 + rank, it1 + rank, rank, -omega);
                    }
                    buffer[layers] = omega;
                    backwardIter(in_out, conv_len, rank * 2, ibuffer, idx);
                }

                static void convolutionCyclic(IntType in_out[], IntType in[], size_t len1, size_t len2, size_t conv_len,
                                              ModInt buffer[], ModInt ibuffer[], ModInt weight) noexcept
                {
                    auto p1 = reinterpret_cast<ModInt *>(in_out), p2 = reinterpret_cast<ModInt *>(in);
                    if (conv_len <= ITER_THRESHOLD)
                    {
                        for (size_t i = 0; i < conv_len; i++)
                        {
                            p1[i] = p1[i].shrinkToMod4();
                            p2[i] = p2[i].shrinkToMod4();
                        }
                        size_t factor = convolutionIterCyclic(p1, p2, conv_len, buffer, ibuffer);
                        ModInt inv = ModInt(factor).inv() * weight;
                        inv *= ModInt::montR();
                        for (size_t i = 0; i < conv_len; i++)
                        {
                            p1[i] = (inv * p1[i]).norm1();
                        }
                        return;
                    }
                    const size_t stride = conv_len / 4;
                    ModInt omega = table.omega1();
                    auto forwardIter = [stride, omega](ModInt inout[], size_t len)
                    {
                        auto it0 = inout, it1 = it0 + stride, it2 = it1 + stride, it3 = it2 + stride;
                        const auto end3 = inout + len, end2 = std::min(end3, it2 + stride);
                        const auto end1 = std::min(end3, it1 + stride), end0 = std::min(end3, it0 + stride);
                        for (; it3 < end3; it0++, it1++, it2++, it3++)
                        {
                            auto t0 = it0[0].shrinkToMod4().norm2(), t1 = it1[0].shrinkToMod4().norm2();
                            auto t2 = it2[0].shrinkToMod4().norm2(), t3 = it3[0].shrinkToMod4().norm2();
                            transform2(t0, t2);
                            auto diff = t1.sub(t3);
                            t1 = t1 + t3;
                            t3 = diff * omega;
                            it0[0] = t0.add(t1), it1[0] = t0.sub(t1), it2[0] = t2.add(t3), it3[0] = t2.sub(t3);
                        }
                        for (; it2 < end2; it0++, it1++, it2++, it3++)
                        {
                            auto t0 = it0[0].shrinkToMod4().norm2(), t1 = it1[0].shrinkToMod4().norm2();
                            auto t2 = it2[0].shrinkToMod4().norm2(), t3 = t1;
                            transform2(t0, t2);
                            t3 = t3 * omega;
                            it0[0] = t0.add(t1), it1[0] = t0.sub(t1), it2[0] = t2.add(t3), it3[0] = t2.sub(t3);
                        }
                        for (; it1 < end1; it0++, it1++, it2++, it3++)
                        {
                            auto t0 = it0[0].shrinkToMod4().norm2(), t1 = it1[0].shrinkToMod4().norm2();
                            auto t2 = t0, t3 = t1;
                            t3 = t3 * omega;
                            it0[0] = t0.add(t1), it1[0] = t0.sub(t1), it2[0] = t2.add(t3), it3[0] = t2.sub(t3);
                        }
                        for (; it0 < end0; it0++, it1++, it2++, it3++)
                        {
                            auto t0 = it0[0].shrinkToMod4();
                            it0[0] = t0, it1[0] = t0, it2[0] = t0, it3[0] = t0;
                        }
                        if (it0 < inout + stride)
                        {
                            const size_t size = (inout + stride - it0) * sizeof(ModInt);
                            std::memset(&it0[0], 0, size);
                            std::memset(&it1[0], 0, size);
                            std::memset(&it2[0], 0, size);
                            std::memset(&it3[0], 0, size);
                        }
                    };
                    forwardIter(p1, len1);
                    if (p2 != p1)
                    {
                        forwardIter(p2, len2);
                    }
                    size_t ret = convolutionCyclic(p1, p2, stride, buffer, ibuffer);
                    convolution(p1 + stride, p2 + stride, stride, buffer, ibuffer, 1);
                    convolution(p1 + stride * 2, p2 + stride * 2, stride, buffer, ibuffer, 2);
                    convolution(p1 + stride * 3, p2 + stride * 3, stride, buffer, ibuffer, 3);
                    auto it0 = p1, it1 = it0 + stride, it2 = it1 + stride, it3 = it2 + stride, end = it1;
                    ModInt omega_inv = itable.omega1();
                    ModInt inv = ModInt(ret * 4).inv() * weight;
                    inv *= ModInt::montR();
                    for (; it0 < end; it0++, it1++, it2++, it3++)
                    {
                        auto t0 = it0[0], t1 = it1[0], t2 = it2[0], t3 = it3[0];
                        transform2(t0, t1);
                        auto diff = t2.sub(t3);
                        t2 = t2 + t3;
                        t3 = diff * omega_inv;
                        it0[0] = (t0.add(t2) * inv).norm1(), it1[0] = (t1.add(t3) * inv).norm1();
                        it2[0] = (t0.sub(t2) * inv).norm1(), it3[0] = (t1.sub(t3) * inv).norm1();
                    }
                }
                static size_t convolutionCyclic(ModInt in_out[], ModInt in[], size_t conv_len,
                                                ModInt buffer[], ModInt ibuffer[]) noexcept
                {
                    if (conv_len <= ITER_THRESHOLD)
                    {
                        return convolutionIterCyclic(in_out, in, conv_len, buffer, ibuffer);
                    }
                    const size_t stride = conv_len / 4;
                    ModInt omega = buffer[0] = table.omega1();
                    auto forwardIter = [stride, omega](ModInt inout[])
                    {
                        auto it0 = inout, it1 = it0 + stride, it2 = it1 + stride, it3 = it2 + stride, end = it1;
                        for (; it0 < end; it0++, it1++, it2++, it3++)
                        {
                            auto t0 = it0[0].norm2(), t1 = it1[0].norm2(), t2 = it2[0].norm2(), t3 = it3[0].norm2();
                            transform2(t0, t2);
                            auto diff = t1.sub(t3);
                            t1 = t1 + t3;
                            t3 = diff * omega;
                            it0[0] = t0.add(t1), it1[0] = t0.sub(t1), it2[0] = t2.add(t3), it3[0] = t2.sub(t3);
                        }
                    };
                    forwardIter(in_out);
                    if (in != in_out)
                    {
                        forwardIter(in);
                    }
                    buffer++;
                    ibuffer++;
                    size_t ret = convolutionCyclic(in_out, in, stride, buffer, ibuffer);
                    convolution(in_out + stride, in + stride, stride, buffer, ibuffer, 1);
                    convolution(in_out + stride * 2, in + stride * 2, stride, buffer, ibuffer, 2);
                    convolution(in_out + stride * 3, in + stride * 3, stride, buffer, ibuffer, 3);
                    ibuffer--;
                    ModInt omega_inv = ibuffer[0] = itable.omega1();
                    auto it0 = in_out, it1 = it0 + stride, it2 = it1 + stride, it3 = it2 + stride, end = it1;
                    for (; it0 < end; it0++, it1++, it2++, it3++)
                    {
                        auto t0 = it0[0], t1 = it1[0], t2 = it2[0], t3 = it3[0];
                        transform2(t0, t1);
                        auto diff = t2.sub(t3);
                        t2 = t2 + t3;
                        t3 = diff * omega_inv;
                        it0[0] = t0 + t2, it1[0] = t1 + t3, it2[0] = t0 - t2, it3[0] = t1 - t3;
                    }
                    return ret * 4;
                }
                static void convolution(ModInt in_out[], ModInt in[], size_t conv_len,
                                        ModInt buffer[], ModInt ibuffer[], size_t idx) noexcept
                {
                    if (conv_len <= ITER_THRESHOLD)
                    {
                        convolutionIter(in_out, in, conv_len, buffer, ibuffer, idx);
                        return;
                    }
                    const size_t stride = conv_len / 4;
                    idx *= 2;
                    ModInt omega1 = table.getNext(buffer[0], idx), omega2 = buffer[0] = table.evenToOdd(omega1);
                    ModInt omega0 = (omega1 * omega1).norm1();
                    auto forwardIter = [stride, omega0, omega1, omega2](ModInt inout[])
                    {
                        auto it0 = inout, it1 = it0 + stride, it2 = it1 + stride, it3 = it2 + stride, end = it1;
                        for (; it0 < end; it0++, it1++, it2++, it3++)
                        {
                            auto t0 = it0[0], t1 = it1[0], t2 = it2[0], t3 = it3[0];
                            dit_butterfly2(t0, t2, omega0);
                            dit_butterfly2(t1, t3, omega0);
                            dit_butterfly2(t0, t1, omega1);
                            dit_butterfly2(t2, t3, omega2);
                            it0[0] = t0, it1[0] = t1, it2[0] = t2, it3[0] = t3;
                        }
                    };
                    forwardIter(in_out);
                    if (in != in_out)
                    {
                        forwardIter(in);
                    }
                    buffer++;
                    ibuffer++;
                    idx *= 2;
                    convolution(in_out, in, stride, buffer, ibuffer, idx);
                    convolution(in_out + stride, in + stride, stride, buffer, ibuffer, idx + 1);
                    convolution(in_out + stride * 2, in + stride * 2, stride, buffer, ibuffer, idx + 2);
                    convolution(in_out + stride * 3, in + stride * 3, stride, buffer, ibuffer, idx + 3);
                    ibuffer--;
                    idx /= 2;
                    ModInt omega_inv1 = itable.getNext(ibuffer[0], idx), omega_inv2 = ibuffer[0] = itable.evenToOdd(omega_inv1);
                    ModInt omega_inv0 = (omega_inv1 * omega_inv1).norm1();
                    auto it0 = in_out, it1 = it0 + stride, it2 = it1 + stride, it3 = it2 + stride, end = it1;
                    for (; it0 < end; it0++, it1++, it2++, it3++)
                    {
                        auto t0 = it0[0], t1 = it1[0], t2 = it2[0], t3 = it3[0];
                        dif_butterfly2(t0, t1, omega_inv1);
                        dif_butterfly2(t2, t3, omega_inv2);
                        dif_butterfly2(t0, t2, omega_inv0);
                        dif_butterfly2(t1, t3, omega_inv0);
                        it0[0] = t0, it1[0] = t1, it2[0] = t2, it3[0] = t3;
                    }
                }
                static void convolutionLinear(const ModInt in1[], const ModInt in2[], ModInt out[], size_t conv_len) noexcept
                {
                    if (0 == conv_len)
                    {
                        return;
                    }
                    ModInt x = in1[0].norm2().norm1();
                    for (size_t i = 0; i < conv_len; i++)
                    {
                        out[i] = in2[i] * x;
                    }
                    const size_t last_idx = conv_len - 1;
                    for (size_t i = 1; i < conv_len; i++)
                    {
                        x = in1[i].norm2().norm1();
                        auto out_it = out + i;
                        for (size_t j = 0; j < last_idx; j++)
                        {
                            out_it[j] += in2[j] * x;
                        }
                        out_it[last_idx] = in2[last_idx] * x;
                    }
                }
                // in_out = in_out * in % (x ^ conv_len - 1)
                static void convolutionCyclicShort(ModInt in_out[], ModInt in[], size_t conv_len) noexcept
                {
                    ModInt temp[SHORT_THRESHOLD * 2];
                    const size_t rem_len = conv_len - conv_len % 4;
                    convolutionLinear(in_out, in, temp, conv_len);
                    const size_t last_idx = conv_len - 1;
                    auto cyclic_it = temp + conv_len;
                    for (size_t i = 0; i < last_idx; i++)
                    {
                        in_out[i] = (temp[i] + cyclic_it[i]);
                    }
                    in_out[last_idx] = temp[last_idx];
                }
                // in_out = in_out * in % (x ^ conv_len - omega)
                static void convolutionShort(ModInt in_out[], ModInt in[], size_t conv_len, ModInt omega) noexcept
                {
                    ModInt temp[SHORT_THRESHOLD * 2];
                    const size_t rem_len = conv_len - conv_len % 4;
                    convolutionLinear(in_out, in, temp, conv_len);
                    const size_t last_idx = conv_len - 1;
                    auto cyclic_it = temp + conv_len;
                    for (size_t i = 0; i < last_idx; i++)
                    {
                        in_out[i] = temp[i] + cyclic_it[i] * omega;
                    }
                    in_out[last_idx] = temp[last_idx];
                }
                static constexpr bool checkConvLen(size_t conv_len) noexcept
                {
                    int tail_zeros = hint_ctz(conv_len);
                    size_t head = conv_len >> tail_zeros;
                    if (head > SHORT_THRESHOLD)
                    {
                        return false;
                    }
                    while (head * 2 <= SHORT_THRESHOLD && tail_zeros > 0)
                    {
                        head *= 2;
                        tail_zeros--;
                    }
                    return (size_t(1) << tail_zeros) <= NTT_MAX_LEN;
                }
                static constexpr size_t findFitLen(size_t len, size_t factor_range) noexcept
                {
                    const size_t conv_len = hint::int_ceil2(len);
                    size_t result = conv_len;
                    int shift = 2;
                    for (size_t i = 3; i <= factor_range; i += 2)
                    {
                        if ((size_t(1) << shift) <= i)
                        {
                            shift++;
                        }
                        size_t try_len = (conv_len >> shift) * i;
                        if (try_len >= len && try_len < result)
                        {
                            result = try_len;
                        }
                    }
                    return result;
                }
            };
            template <uint64_t MOD, uint64_t ROOT>
            const typename NTT<MOD, ROOT>::TableType NTT<MOD, ROOT>::table{NTT<MOD, ROOT>::root(), 1, 2};
            template <uint64_t MOD, uint64_t ROOT>
            const typename NTT<MOD, ROOT>::TableType NTT<MOD, ROOT>::itable{NTT<MOD, ROOT>::rootInv(), 1, 2};
            template <uint64_t MOD, uint64_t ROOT>
            constexpr int NTT<MOD, ROOT>::MOD_BITS;
            template <uint64_t MOD, uint64_t ROOT>
            constexpr int NTT<MOD, ROOT>::MAX_LOG_LEN;
            template <uint64_t MOD, uint64_t ROOT>
            constexpr size_t NTT<MOD, ROOT>::NTT_MAX_LEN;
            template <uint64_t MOD, uint64_t ROOT>
            constexpr size_t NTT<MOD, ROOT>::L1_BYTE;
            template <uint64_t MOD, uint64_t ROOT>
            constexpr size_t NTT<MOD, ROOT>::LONG_BYTE;
            template <uint64_t MOD, uint64_t ROOT>
            constexpr size_t NTT<MOD, ROOT>::ITER_THRESHOLD;
            template <uint64_t MOD, uint64_t ROOT>
            constexpr size_t NTT<MOD, ROOT>::LONG_THRESHOLD;
            template <uint64_t MOD, uint64_t ROOT>
            constexpr size_t NTT<MOD, ROOT>::SHORT_THRESHOLD;

            using NTT64_1 = NTT<MOD0, ROOT0>;
            using NTT64_2 = NTT<MOD1, ROOT1>;
            using NTT64_3 = NTT<MOD2, ROOT2>;
            using NTT32_1 = NTT<MOD3, ROOT3>;
            using NTT32_2 = NTT<MOD4, ROOT4>;
            using NTT32_3 = NTT<MOD5, ROOT5>;
        }

        namespace fft
        {
            // TODO: implement
        }
    }
    namespace arithm
    {
        namespace addition
        {
            // Addition, return carry
            template <typename NumTy, typename Executor>
            constexpr bool abs_add_long(const NumTy a[], const NumTy b[], size_t len, NumTy sum[], const Executor &exec)
            {
                bool carry = false;
                size_t i = 0;
                for (const size_t rem_len = len - len % 4; i < rem_len; i += 4)
                {
                    exec.addCarryX4(a + i, b + i, sum + i, carry);
                }
                for (; i < len; i++)
                {
                    sum[i] = exec.addCarry(a[i], b[i], carry);
                }
                return carry;
            }

            // BigInt a + num, return carry
            template <typename NumTy, typename Executor>
            constexpr bool abs_add_long_num(const NumTy a[], size_t len, NumTy num, NumTy sum[], const Executor &exec)
            {
                assert(len > 0);
                bool carry = false;
                sum[0] = exec.addHalf(a[0], num, carry);
                a++, sum++, len--;
                size_t i = 0;
                for (const size_t rem_len = len - len % 4; i < rem_len; i += 4)
                {
                    exec.addCfX4(a + i, sum + i, carry);
                }
                for (; i < len; i++)
                {
                    sum[i] = exec.addCf(a[i], carry);
                }
                return carry;
            }

            template <typename NumTy, typename Executor>
            constexpr bool abs_add(const NumTy a[], size_t len_a, const NumTy b[], size_t len_b, NumTy sum[],
                                   const Executor &exec, bool assign_carry = true)
            {
                if (len_a < len_b)
                {
                    std::swap(a, b);
                    std::swap(len_a, len_b);
                }
                bool carry = abs_add_long(a, b, len_b, sum, exec);
                if (len_a > len_b)
                {
                    carry = abs_add_long_num(a + len_b, len_a - len_b, NumTy(carry), sum + len_b, exec);
                }
                if (assign_carry)
                {
                    sum[len_a] = carry;
                }
                return carry;
            }

            // Subtraction, return borrow
            template <typename NumTy, typename Executor>
            constexpr bool abs_sub_long(const NumTy a[], const NumTy b[], size_t len, NumTy diff[], const Executor &exec)
            {
                bool borrow = false;
                size_t i = 0;
                for (const size_t rem_len = len - len % 4; i < rem_len; i += 4)
                {
                    exec.subBorrowX4(a + i, b + i, diff + i, borrow);
                }
                for (; i < len; i++)
                {
                    diff[i] = exec.subBorrow(a[i], b[i], borrow);
                }
                return borrow;
            }

            // BigInt a - num, return borrow
            template <typename NumTy, typename Executor>
            constexpr bool abs_sub_long_num(const NumTy a[], size_t len, NumTy num, NumTy diff[], const Executor &exec)
            {
                assert(len > 0);
                bool borrow = false;
                diff[0] = exec.subHalf(a[0], num, borrow);
                a++, diff++, len--;
                size_t i = 0;
                for (const size_t rem_len = len - len % 4; i < rem_len; i += 4)
                {
                    exec.subBfX4(a + i, diff + i, borrow);
                }
                for (; i < len; i++)
                {
                    diff[i] = exec.subBf(a[i], borrow);
                }
                return borrow;
            }

            template <typename NumTy, typename Executor>
            constexpr bool abs_sub(const NumTy a[], size_t len_a, const NumTy b[], size_t len_b, NumTy diff[],
                                   const Executor &exec)
            {
                assert(len_a >= len_b);
                bool borrow = abs_sub_long(a, b, len_b, diff, exec);
                if (len_a > len_b)
                {
                    borrow = abs_sub_long_num(a + len_b, len_a - len_b, NumTy(borrow), diff + len_b, exec);
                }
                return borrow;
            }

            template <typename NumTy, typename Executor>
            [[nodiscard]] constexpr int abs_difference(const NumTy a[], size_t len1, const NumTy b[], size_t len2, NumTy diff[],
                                                       const Executor &exec)
            {
                auto cmp_ab = utility::abs_compare_with_length(a, len1, b, len2);
                if (len1 == len2)
                {
                    std::fill(diff + cmp_ab.diff_len, diff + len1, NumTy(0));
                    len1 = len2 = cmp_ab.diff_len;
                }
                if (cmp_ab.cmp < 0)
                {
                    std::swap(a, b);
                    std::swap(len1, len2);
                }
                bool borrow = abs_sub(a, len1, b, len2, diff, exec);
                assert(!borrow);
                return cmp_ab.cmp;
            }

            // Binary absolute addtion a+b=sum, return the carry
            template <typename UintTy>
            constexpr void abs_add_binary(const UintTy a[], size_t len_a, const UintTy b[], size_t len_b, UintTy sum[])
            {
                constexpr utility::BaseExecutorBinary<UintTy> exec;
                abs_add(a, len_a, b, len_b, sum, exec);
            }
            template <typename UintTy>
            constexpr void abs_add_base(const UintTy a[], size_t len_a, const UintTy b[], size_t len_b, UintTy sum[], UintTy base_num)
            {
                const utility::BaseExecutor<UintTy> exec(base_num);
                abs_add(a, len_a, b, len_b, sum, exec);
            }
            template <typename UintTy>
            constexpr void abs_sub_binary(const UintTy a[], size_t len_a, const UintTy b[], size_t len_b, UintTy diff[])
            {
                constexpr utility::BaseExecutorBinary<UintTy> exec;
                abs_sub(a, len_a, b, len_b, diff, exec);
            }
            template <typename UintTy>
            constexpr void abs_sub_base(const UintTy a[], size_t len_a, const UintTy b[], size_t len_b, UintTy diff[], UintTy base_num)
            {
                const utility::BaseExecutor<UintTy> exec(base_num);
                abs_sub(a, len_a, b, len_b, diff, exec);
            }
        }
        namespace multiplication
        {
            using namespace utility;
            constexpr size_t STACK_MAX_LEN = 1024;
            constexpr size_t BASIC_THRESHOLD = 32;
            constexpr size_t KARATSUBA_THRESHOLD = 2048;
            // out = in * num_mul + num_add, return carry
            template <typename NumTy, typename Executor>
            [[nodiscard]] inline NumTy abs_mul_add_num(const NumTy in[], size_t len, NumTy out[], NumTy num_mul, NumTy num_add,
                                                       const Executor &exec)
            {
                if (0 == num_mul)
                {
                    std::fill(out, out + len, num_add);
                    return len > 0 ? 0 : num_add;
                }
                size_t i = 0;
                for (const size_t rem_len = len - len % 4; i < rem_len; i += 4)
                {
                    num_add = exec.mulAddX4(in + i, num_mul, num_add, out + i);
                }
                for (; i < len; i++)
                {
                    exec.mulAdd(in[i], num_mul, num_add, out[i], num_add);
                }
                return num_add;
            }

            // out = out + in * num_mul + num_add, return carry
            template <typename NumTy, typename Executor>
            [[nodiscard]] inline NumTy abs_add_equal_mul_add_num(const NumTy in[], size_t len, NumTy out[], NumTy num_mul, NumTy num_add,
                                                                 const Executor &exec)
            {
                if (0 == len)
                {
                    return num_add;
                }
                if (0 == num_mul)
                {
                    if (0 == num_add)
                    {
                        return 0;
                    }
                    return addition::abs_add_long_num(out, len, num_add, out, exec);
                }
                size_t rem = len % 8;
                num_add = exec.addEqMulAddX8(in, len - rem, out, num_mul, num_add);
                if (rem > 0)
                {
                    auto p_in = in + len - 7;
                    auto p_out = out + len - 7;
                    switch (rem)
                    {
                    case 7:
                        exec.mulAddPlusSelf(p_in[0], num_mul, p_out[0], num_add);
                    case 6:
                        exec.mulAddPlusSelf(p_in[1], num_mul, p_out[1], num_add);
                    case 5:
                        exec.mulAddPlusSelf(p_in[2], num_mul, p_out[2], num_add);
                    case 4:
                        exec.mulAddPlusSelf(p_in[3], num_mul, p_out[3], num_add);
                    case 3:
                        exec.mulAddPlusSelf(p_in[4], num_mul, p_out[4], num_add);
                    case 2:
                        exec.mulAddPlusSelf(p_in[5], num_mul, p_out[5], num_add);
                    case 1:
                        exec.mulAddPlusSelf(p_in[6], num_mul, p_out[6], num_add);
                    default:
                        break;
                    }
                }
                return num_add;
            }

            inline void abs_mul_basic_bin(const uint64_t in1[], size_t len1, const uint64_t in2[], size_t len2, uint64_t out[])
            {
                auto p1 = reinterpret_cast<const uint32_t *>(in1);
                auto p2 = reinterpret_cast<const uint32_t *>(in2);
                auto p_out = reinterpret_cast<uint32_t *>(out);
                len1 *= 2;
                len2 *= 2;
                std::vector<uint32_t> buf(len1 + len2);
                for (size_t i = 0; i < len1; i++)
                {
                    uint64_t num = p1[i], carry = 0;
                    for (size_t j = 0; j < len2; j++)
                    {
                        uint64_t prod = num * p2[j] + carry + buf[i + j];
                        buf[i + j] = prod;
                        carry = prod >> 32;
                    }
                    buf[i + len2] = carry;
                }
                std::copy(buf.begin(), buf.begin() + len1 + len2, p_out);
            }

            inline void mul_check(const uint64_t in1[], size_t len1, const uint64_t in2[], size_t len2, const uint64_t out[])
            {
                if (len1 * len2 == 0)
                {
                    return;
                }
                std::vector<uint64_t> buf(len1 + len2);
                abs_mul_basic_bin(in1, len1, in2, len2, buf.data());
                if (!std::equal(buf.begin(), buf.begin() + len1 + len2, out))
                {
                    ary_print(in1, len1);
                    ary_print(in2, len2);
                    for (size_t i = 0; i < len1 + len2; i++)
                    {
                        std::cout << "buf[" << i << "] = " << buf[i] << " out[" << i << "] = " << out[i] << std::endl;
                    }
                    exit(EXIT_FAILURE);
                }
            }

            // Basic multiplication algorithm, use_buf = true if the in and out overlap
            template <typename NumTy, typename WorkMem = std::vector<NumTy>, typename Executor>
            inline void abs_mul_basic(const NumTy in1[], size_t len1, const NumTy in2[], size_t len2, NumTy out[],
                                      const Executor &exec)
            {
                if (len1 < len2)
                {
                    std::swap(in1, in2);
                    std::swap(len1, len2); // Let in1 be the loonger one
                }
                if (0 == len2 || nullptr == in1 || nullptr == in2)
                {
                    return;
                }
                WorkMem in1_temp_dynamic;
                NumTy in1_temp[STACK_MAX_LEN], in2_temp[STACK_MAX_LEN];
                bool use_buf = mem_overlap(in1, in1 + len1, out, out + len1 + len2) ||
                               mem_overlap(in2, in2 + len2, out, out + len1 + len2);
                if (1 == len2)
                {
                    use_buf = (in1 < out) && (out < in1 + len1);
                }
                if (use_buf)
                {
                    assert(len2 <= STACK_MAX_LEN);
                    std::copy(in2, in2 + len2, in2_temp);
                    in2 = in2_temp;
                    if (len1 <= STACK_MAX_LEN)
                    {
                        std::copy(in1, in1 + len1, in1_temp);
                        in1 = in1_temp;
                    }
                    else
                    {
                        in1_temp_dynamic.assign(in1, in1 + len1);
                        in1 = in1_temp_dynamic.data();
                    }
                }
                if (1 == len2)
                {
                    out[len1] = abs_mul_add_num(in1, len1, out, in2[0], NumTy(0), exec);
                    return;
                }
                if (len2 <= STACK_MAX_LEN)
                {
                    std::swap(in1, in2);
                    std::swap(len1, len2); // Let in2 be the loonger one
                }
                out[len2] = abs_mul_add_num(in2, len2, out, in1[0], NumTy(0), exec);
                for (size_t i = 1; i < len1; i++)
                {
                    out[len2 + i] = abs_add_equal_mul_add_num(in2, len2, out + i, in1[i], NumTy(0), exec);
                }
            }
            // Karatsuba multiplication algorithm
            template <typename NumTy, typename WorkMem = std::vector<NumTy>, typename Executor>
            inline void abs_mul_karatusba(const NumTy in1[], size_t len1, const NumTy in2[], size_t len2, NumTy out[],
                                          const Executor &exec, NumTy *work_begin = nullptr, NumTy *work_end = nullptr)
            {
                if (0 == len1 || 0 == len2 || nullptr == in1 || nullptr == in2)
                {
                    return;
                }
                const size_t out_len = get_mul_len(len1, len2);
                {
                    const size_t true_len1 = count_ture_length(in1, len1), true_len2 = count_ture_length(in2, len2);
                    const size_t trail_len1 = count_trailing_zero(in1, true_len1), trail_len2 = count_trailing_zero(in2, true_len2);
                    if (true_len1 + true_len2 < out_len || trail_len1 + trail_len2 > 0)
                    {
                        auto in1_temp = in1 + trail_len1, in2_temp = in2 + trail_len2;
                        auto out_temp = out + trail_len1 + trail_len2;
                        abs_mul_karatusba(in1_temp, true_len1 - trail_len1, in2_temp, true_len2 - trail_len2, out_temp, exec, work_begin, work_end);
                        std::fill(out, out_temp, NumTy(0));
                        std::fill(out + get_mul_len(true_len1, true_len2), out + out_len, NumTy(0));
                        return;
                    }
                }
                if (len1 < len2)
                {
                    std::swap(in1, in2);
                    std::swap(len1, len2); // Let in1 be the loonger one
                }
                if (len2 <= BASIC_THRESHOLD)
                {
                    abs_mul_basic(in1, len1, in2, len2, out, exec);
                    return;
                }
                // Split A * B -> (AH * BASE + AL) * (BH * BASE + BL)
                // (AH * BASE + AL) * (BH * BASE + BL) = AH * BH * BASE^2 + (AH * BL + AL * BH) * BASE + AL * BL
                // Let M = AL * BL, N = AH * BH, K1 = (AH - AL), K2 = (BH - BL), K = K1 * K2 = AH * BH - (AH * BL + AL * BH) + AL * BL
                // A * B = N * BASE^2 + (M + N - K) * BASE + M
                const size_t base_len = (len1 + 1) / 2;
                size_t len1_low = base_len, len1_high = len1 - base_len;
                size_t len2_low = base_len, len2_high = len2 - base_len;
                if (len2 <= base_len)
                {
                    len2_low = len2;
                    len2_high = 0;
                }
                auto in1_high = in1 + len1_low, in2_high = in2 + len2_low;
                // Get length of every part
                size_t m_len = get_mul_len(len1_low, len2_low);
                size_t n_len = get_mul_len(len1_high, len2_high);

                // Get enough work_mem
                WorkMem work_mem;
                const size_t block_len = base_len * 2, work_size = block_len * 2 + 2;
                if (nullptr == work_begin || nullptr == work_end || work_end < work_begin + work_size)
                {
                    work_mem.resize(work_size * 2);
                    work_begin = work_mem.data();
                    work_end = work_begin + work_mem.size();
                }
                // Set pointer of every part
                auto m = work_begin, n = out + base_len * 2, k1 = m + block_len + 2, k2 = k1 + base_len, k = k1;
                int cmp1 = addition::abs_difference(in1, len1_low, in1_high, len1_high, k1, exec); // k1 = abs(AH - AL)
                int cmp2 = addition::abs_difference(in2, len2_low, in2_high, len2_high, k2, exec); // k2 = abs(BH - BL)
                size_t k1_len = get_sub_len(len1_low, len1_high);
                size_t k2_len = get_sub_len(len2_low, len2_high);
                // Compute M,N
                work_begin += work_size;
                abs_mul_karatusba(in1, len1_low, in2, len2_low, m, exec, work_begin, work_end);             // M = AL * BL
                abs_mul_karatusba(in1_high, len1_high, in2_high, len2_high, n, exec, work_begin, work_end); // N = AH * BH
                remove_leading_zeros(m, m_len);
                remove_leading_zeros(n, n_len);
                abs_mul_karatusba(k1, k1_len, k2, k2_len, k, exec, work_begin, work_end); // K = K1 * K2
                size_t k_len = count_ture_length(k, get_mul_len(k1_len, k2_len));
                // Compute K1,K2
                // Combine the result
                // out = M + N * BASE ^ 2
                {
                    std::copy(m, m + m_len, out);
                    std::fill(out + m_len, out + base_len * 2, NumTy(0));
                    std::fill(out + base_len * 2 + n_len, out + out_len, NumTy(0));
                }
                // out = N * BASE^2 + (M + N - K) * BASE + M
                addition::abs_add(m, m_len, n, n_len, m, exec);
                m_len = get_add_len(m_len, n_len);
                if (0 == k_len)
                {
                    // NOP
                }
                else if ((cmp1 > 0) == (cmp2 > 0))
                {
                    bool borrow = addition::abs_sub(m, m_len, k, k_len, m, exec);
                    assert(!borrow);
                }
                else
                {
                    addition::abs_add(m, m_len, k, k_len, m, exec);
                    m_len = get_add_len(m_len, k_len);
                }
                remove_leading_zeros(m, m_len);
                auto out_base = out + base_len;
                bool carry = addition::abs_add(m, m_len, out_base, out_len - base_len, out_base, exec, false);
                assert(!carry);
            }

            template <typename NumTy, typename Executor>
            void abs_mul_ntt(const NumTy in1[], size_t len1, const NumTy in2[], size_t len2, NumTy out[],
                             const Executor &exec) noexcept
            {
                if (0 == len1 || 0 == len2 || nullptr == in1 || nullptr == in2)
                {
                    return;
                }
                using NTT1 = transform::ntt::NTT64_1;
                using NTT2 = transform::ntt::NTT64_2;
                using NTT3 = transform::ntt::NTT64_3;
                using ModInt1 = NTT1::ModInt;
                using ModInt2 = NTT2::ModInt;
                using ModInt3 = NTT3::ModInt;
                constexpr NTT1::IntType MOD1 = NTT1::mod();
                constexpr NTT2::IntType MOD2 = NTT2::mod();
                constexpr NTT3::IntType MOD3 = NTT3::mod();
                constexpr auto mod23_mod1 = mulMod(MOD2 % MOD1, MOD3 % MOD1, MOD1);
                constexpr auto mod13_mod2 = mulMod(MOD1 % MOD2, MOD3 % MOD2, MOD2);
                constexpr auto mod12_mod3 = mulMod(MOD1 % MOD3, MOD2 % MOD3, MOD3);
                constexpr ModInt1 inv1 = hint::mod_inv<int64_t>(mod23_mod1, MOD1);
                constexpr ModInt2 inv2 = hint::mod_inv<int64_t>(mod13_mod2, MOD2);
                constexpr ModInt3 inv3 = hint::mod_inv<int64_t>(mod12_mod3, MOD3);

                constexpr ModInt1 mod23_1(NTT1::mod());
                size_t len = len1 + len2 - 1;
                size_t conv_len = NTT1::findFitLen(len);
                std::vector<uint64_t> conv1(conv_len), temp(conv_len);
                std::copy(in1, in1 + len1, conv1.data());
                std::copy(in2, in2 + len2, temp.data());
                NTT1::convolution(conv1.data(), temp.data(), len1, len2, conv_len, inv1);

                std::vector<uint64_t> conv2(conv_len);
                std::copy(in1, in1 + len1, conv2.data());
                std::copy(in2, in2 + len2, temp.data());
                NTT2::convolution(conv2.data(), temp.data(), len1, len2, conv_len, inv2);

                std::vector<uint64_t> conv3(conv_len);
                std::copy(in1, in1 + len1, conv3.data());
                std::copy(in2, in2 + len2, temp.data());
                NTT3::convolution(conv3.data(), temp.data(), len1, len2, conv_len, inv3);

                extend_int::Uint192 carry{};
                for (size_t i = 0; i < len; i++)
                {
                    carry += CRTMod3T<MOD1, MOD2, MOD3>::crt3(conv1[i], conv2[i], conv3[i]);
                    out[i] = exec.divRemBase(carry);
                }
                out[len] = uint64_t(carry);
            }
        }
        namespace division
        {
            // out = in / divisor, return remainder
            template <typename NumTy, typename Executor>
            inline NumTy abs_div_num(const NumTy in[], size_t len, NumTy out[], NumTy divisor, Executor &exec)
            {
                assert(divisor > 0);
                assert(exec.checkDivisor(divisor));
                const utility::DivExecutor<NumTy> div_exe{divisor};
                size_t i = len;
                NumTy remainder = 0;
                while (i > 0)
                {
                    i--;
                    NumTy hi, lo;
                    exec.dualBaseToBin(remainder, in[i], lo, hi);
                    out[i] = div_exe.divRem(hi, lo, remainder);
                }
                return remainder;
            }

            template <typename NumTy, typename Executor>
            inline void abs_div_basic(NumTy dividend[], size_t len1, const NumTy divisor[], size_t len2,
                                      NumTy quotient[], Executor &exec)
            {
                len2 = utility::count_ture_length(divisor, len2);
                assert(len2 > 0);
                assert(divisor[len2 - 1] >= exec.halfBase());
                if (len2 == 1)
                {
                    NumTy rem = abs_div_num_norm(dividend, len1, quotient, divisor[0], exec);
                    dividend[0] = rem;
                    return;
                }
                assert(len1 >= len2);
            }

            template <typename NumTy, typename Executor>
            inline void abs_div_rec(const NumTy in[], size_t len, NumTy out[], NumTy divisor, Executor &exec)
            {
            }

            template <typename NumTy, typename Executor>
            inline void abs_div(const NumTy in[], size_t len, NumTy out[], NumTy divisor, Executor &exec)
            {
            }
        }
        namespace base_conversion
        {

        }
    }
    using arithm::multiplication::BASIC_THRESHOLD;
    using arithm::multiplication::KARATSUBA_THRESHOLD;
    template <typename NumTy, typename Container, typename BaseExecutor, NumTy BASE>
    class HyperIntImpl
    {
    public:
        using Num = NumTy;
        using HyperUint = HyperIntImpl;

        HyperIntImpl() = default;
        HyperIntImpl(const NumTy &n)
        {
            data.reserve(2);
            data.push_back(n);
        }
        HyperIntImpl(const HyperUint &other) = default;
        HyperIntImpl(HyperUint &&other) = default;

        HyperUint &operator=(const HyperUint &other) = default;
        HyperUint &operator=(HyperUint &&other) = default;

        bool isZero() const
        {
            return data.size() == 0;
        }

        size_t limbSize() const
        {
            return data.size();
        }

        NumTy *limbPtr()
        {
            return data.data();
        }

        const NumTy *limbPtr() const
        {
            return data.data();
        }

        void fromString(const std::string &str)
        {
            for (auto c : str)
            {
                if (c >= '0' && c <= '9')
                {
                    this->mulAddNum(10, c - '0');
                }
            }
        }

        std::string toString() const
        {
            std::string str;
            HyperUint temp = *this;
            while (!temp.isZero())
            {
                NumTy rem = temp.divRemNum(10);
                str += char('0' + rem);
            }
            std::reverse(str.begin(), str.end());
            return str;
        }

        HyperUint &operator+=(const HyperUint &other)
        {
            size_t len1 = limbSize(), len2 = other.limbSize();
            size_t add_len = utility::get_add_len(len1, len2);
            data.resize(add_len);
            arithm::addition::abs_add(limbPtr(), len1, other.limbPtr(), len2, limbPtr(), exec);
            shrinkLeadingZeros();
            return *this;
        }

        // Crash when other is bigger than *this
        HyperUint &operator-=(const HyperUint &other)
        {
            size_t len1 = data.size(), len2 = other.data.size();
            assert(len1 >= len2);
            size_t sub_len = utility::get_sub_len(len1, len2);
            data.resize(sub_len);
            bool borrow = arithm::addition::abs_sub(limbPtr(), len1, other.limbPtr(), len2, limbPtr(), exec);
            assert(!borrow);
            shrinkLeadingZeros();
            return *this;
        }

        HyperUint &operator*=(const HyperUint &other)
        {
            size_t len1 = data.size(), len2 = other.data.size();
            size_t mul_len = utility::get_mul_len(len1, len2);
            data.resize(mul_len);
            if (mul_len <= KARATSUBA_THRESHOLD || len1 <= BASIC_THRESHOLD || len2 <= BASIC_THRESHOLD)
            {
                arithm::multiplication::abs_mul_karatusba(limbPtr(), len1, other.limbPtr(), len2, limbPtr(), exec);
            }
            else
            {
                arithm::multiplication::abs_mul_ntt(limbPtr(), len1, other.limbPtr(), len2, limbPtr(), exec);
            }
            shrinkLeadingZeros();
            return *this;
        }

        HyperUint &operator*=(Num other)
        {
            size_t len = data.size();
            size_t mul_len = len + 1;
            data.resize(mul_len);
            data[len] = arithm::multiplication::abs_mul_add_num(limbPtr(), len, limbPtr(), other, NumTy{}, exec);
            shrinkLeadingZeros();
            return *this;
        }

        void mulAddNum(Num num_mul, Num num_add)
        {
            size_t len = data.size();
            size_t res_len = len + 1;
            data.resize(res_len);
            data[len] = arithm::multiplication::abs_mul_add_num(limbPtr(), len, limbPtr(), num_mul, num_add, exec);
            shrinkLeadingZeros();
        }

        Num divRemNum(Num divisor)
        {
            size_t len = data.size();
            NumTy rem = arithm::division::abs_div_num(limbPtr(), len, limbPtr(), divisor, exec);
            shrinkLeadingZeros();
            return rem;
        }

        void shrinkLeadingZeros()
        {
            data.resize(utility::count_ture_length(limbPtr(), data.size()));
        }

        // HyperUint &operator/=(const HyperUint &other)
        // {
        // }

    private:
        using Vec = Container;
        using WorkMem = Container;
        Vec data;

        static constexpr BaseExecutor exec{BASE};
    };
    template <typename NumTy, typename Container, typename BaseExecutor, NumTy BASE>
    constexpr BaseExecutor HyperIntImpl<NumTy, Container, BaseExecutor, BASE>::exec;

    using HyperIntHex = HyperIntImpl<uint64_t, std::vector<uint64_t>, utility::BaseExecutorBinary<uint64_t>, 64>;
}
#endif