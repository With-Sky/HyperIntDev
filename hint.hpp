#ifndef HINT_HPP
#define HINT_HPP

#include <iostream>
#include <algorithm>
#include <random>
#include <string>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include <cassert>
#include <climits>
#include <cstdint>

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
        x &= (-x);
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
            if (b == 0)
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
        if (b == 0)
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
    template <typename UintTy>
    constexpr int hint_log2(UintTy n)
    {
        constexpr int bits = sizeof(UintTy) * CHAR_BIT;
        constexpr UintTy MASK = all_one<UintTy>(bits / 2) << (bits / 2);
        UintTy m = MASK;
        int res = 0, shift = bits / 2;
        while (shift > 0)
        {
            if ((n & m))
            {
                res += shift;
                n >>= shift;
            }
            shift /= 2;
            m >>= shift;
        }
        return res;
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
            __uint128_t x(a);
            x *= b;
            low = uint64_t(x), high = uint64_t(x >> 64);
#else
#if defined(HINT_WIN64) // Has _umul128
#pragma message("Using _umul128 to compute 64bit x 64bit to 128bit")
            HintULL lo, hi;
            lo = _umul128(a, b, &hi);
            low = lo, high = hi;
#else // No _umul128 or __uint128_t
#pragma message("Using basic function to compute 64bit x 64bit to 128bit")
            mul64x64to128_base(a, b, low, high);
#endif
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

        constexpr uint64_t div_dword_word(uint64_t dividend_hi64, uint64_t dividend_lo64, uint64_t divisor, uint64_t &remainder)
        {
            uint64_t quot = div128by64to64(dividend_hi64, dividend_lo64, divisor);
            remainder = dividend_lo64;
            return quot;
        }

        template <typename UintTy>
        inline UintTy div_dword_word(UintTy dividend_hi64, UintTy dividend_lo64, UintTy divisor, UintTy &remainder)
        {
            uint64_t dividend = (uint64_t(dividend_hi64) << 32) | dividend_lo64;
            remainder = dividend % divisor;
            return dividend / divisor;
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

        // template <typename NumTy, typename ProdTy>
        // class DivExecutor
        // {
        // public:
        //     constexpr DivExecutor(NumTy divisor_in) : divisor(divisor_in)
        //     {
        //         inv = getInv(divisor, shift);
        //         divisor <<= shift;
        //     }
        //     // Return dividend / divisor, dividend %= divisor
        //     NumTy divRem(ProdTy dividend, NumTy &rem) const
        //     {
        //         dividend <<= shift;
        //         NumTy quot = divRemNorm(dividend, rem);
        //         rem >>= shift;
        //         return quot;
        //     }
        //     // Return dividend / divisor, dividend %= divisor
        //     NumTy divRem(NumTy dividend_hi, NumTy dividend_lo, NumTy &rem) const
        //     {
        //         NumTy quot;
        //         if (shift > 0)
        //         {
        //             dividend_hi = (dividend_hi << shift) | dividend_lo >> (NUM_BITS - shift);
        //             dividend_lo <<= shift;
        //             quot = divRemNorm(dividend_hi, dividend_lo, rem);
        //             rem >>= shift;
        //         }
        //         else
        //         {
        //             quot = divRemNorm(dividend_hi, dividend_lo, rem);
        //         }
        //         return quot;
        //     }
        //     NumTy prodDivRem(NumTy a, NumTy b, NumTy &rem) const
        //     {
        //         ProdTy dividend = ProdTy(a << shift) * b;
        //         NumTy quot = this->divRemNorm(dividend, rem);
        //         rem >>= shift;
        //         return quot;
        //     }

        //     // NumTy div(ProdTy dividend) const
        //     // {
        //     //     return divRem(dividend);
        //     // }
        //     // NumTy mod(ProdTy dividend) const
        //     // {
        //     //     divRem(dividend);
        //     //     return dividend;
        //     // }
        //     NumTy divRemNorm(ProdTy dividend, NumTy &rem) const
        //     {
        //         return divRemNorm(dividend >> NUM_BITS, dividend, rem);
        //     }
        //     // Reference:N. Möller and T. Granlund, "Improved Division by Invariant Integers,"
        //     // in IEEE Transactions on Computers, vol. 60, no. 2, pp. 165-175, Feb. 2011, doi: 10.1109/TC.2010.143.
        //     // Best performance, optimized first branch.
        //     NumTy divRemNorm(NumTy dividend_hi, NumTy dividend_lo, NumTy &rem) const
        //     {
        //         NumTy lo, hi, quot, mask;
        //         mul_binary(dividend_hi, inv, lo, hi);
        //         lo += dividend_lo;
        //         hi += dividend_hi + (lo < dividend_lo);
        //         quot = hi + 1;
        //         rem = dividend_lo - quot * divisor;
        //         mask = -NumTy(rem > lo); // mask = -1 if rem > dividend
        //         rem += divisor & mask;   // rem += divisor if rem > dividend
        //         quot += mask;            // quot -= 1 if rem > dividend
        //         if (rem >= divisor)
        //         {
        //             rem -= divisor;
        //             quot++;
        //         }
        //         return quot;
        //     }

        // private:
        //     static constexpr NumTy getInv(NumTy divisor, int &leading_zero)
        //     {
        //         constexpr NumTy MAX = hint::all_one<NumTy>(NUM_BITS);
        //         leading_zero = hint::hint_clz(divisor);
        //         divisor <<= leading_zero;
        //         NumTy rem;
        //         return div_dword_word(MAX - divisor, MAX, divisor, rem);
        //     }

        //     NumTy divisor = 0;
        //     NumTy inv = 0;
        //     int shift = 0;
        //     static constexpr int NUM_BITS = sizeof(NumTy) * CHAR_BIT;
        // };
        // template <typename NumTy, typename ProdTy>
        // constexpr int DivExecutor<NumTy, ProdTy>::NUM_BITS;

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
                NumTy rem;
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
                assert(2 <= base && base <= hint::all_one<NumTy>(NUM_BITS));
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
            constexpr void addCarryX4(const NumTy in1[4], const NumTy in2[4], NumTy sum[4], bool &carry) const
            {
                sum[0] = addCarry(in1[0], in2[0], carry);
                sum[1] = addCarry(in1[1], in2[1], carry);
                sum[2] = addCarry(in1[2], in2[2], carry);
                sum[3] = addCarry(in1[3], in2[3], carry);
            }
            constexpr void subBorrowX4(const NumTy in1[4], const NumTy in2[4], NumTy sum[4], bool &borrow) const
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

            void mulInBase(NumTy a, NumTy b, NumTy &lo, NumTy &hi) const
            {
                hi = div_exe.prodDivRem(a, b, lo);
            }

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

        class BaseExecutorBinary
        {
        public:
            template <typename UintTy>
            static constexpr UintTy addCarry(UintTy a, UintTy b, bool &cf)
            {
                return add_carry(a, b, cf);
            }
            template <typename UintTy>
            static constexpr UintTy addHalf(UintTy a, UintTy b, bool &cf)
            {
                return add_half(a, b, cf);
            }
            template <typename UintTy>
            static constexpr UintTy subBorrow(UintTy a, UintTy b, bool &bf)
            {
                return sub_borrow(a, b, bf);
            }
            template <typename UintTy>
            static constexpr UintTy subHalf(UintTy a, UintTy b, bool &bf)
            {
                return sub_half(a, b, bf);
            }
            template <typename UintTy>
            static constexpr UintTy addCf(UintTy a, bool &cf)
            {
                a += cf;
                cf = a < cf;
                return a;
            }
            template <typename UintTy>
            static constexpr UintTy subBf(UintTy a, bool &bf)
            {
                UintTy diff = a - bf;
                bf = diff > a;
                return diff;
            }
            template <typename UintTy>
            static constexpr void addCarryX4Basic(UintTy in1[4], UintTy in2[4], UintTy sum[4], bool &carry)
            {
                sum[0] = add_carry(in1[0], in2[0], carry);
                sum[1] = add_carry(in1[1], in2[1], carry);
                sum[2] = add_carry(in1[2], in2[2], carry);
                sum[3] = add_carry(in1[3], in2[3], carry);
            }

            template <typename UintTy>
            static constexpr void addCarryX4(const UintTy in1[4], const UintTy in2[4], UintTy sum[4], bool &carry)
            {
                addCarryX4Basic(in1, in2, sum, carry);
            }

            static inline void addCarryX4(const uint64_t in1[4], const uint64_t in2[4], uint64_t sum[4], bool &carry)
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
                addCarryX4Basic(in1, in2, sum, carry);
#endif
            }

            template <typename UintTy>
            static constexpr void subBorrowX4Basic(const UintTy in1[4], const UintTy in2[4], UintTy diff[4], bool &borrow)
            {
                diff[0] = sub_borrow(in1[0], in2[0], borrow);
                diff[1] = sub_borrow(in1[1], in2[1], borrow);
                diff[2] = sub_borrow(in1[2], in2[2], borrow);
                diff[3] = sub_borrow(in1[3], in2[3], borrow);
            }
            template <typename UintTy>
            static constexpr void subBorrowX4(const UintTy in1[4], const UintTy in2[4], UintTy diff[4], bool &borrow)
            {
                subBorrowX4Basic(in1, in2, diff, borrow);
            }

            static inline void subBorrowX4(const uint64_t in1[4], const uint64_t in2[4], uint64_t diff[4], bool &borrow)
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
                subBorrowX4Basic(in1, in2, diff, borrow);
#endif
            }
            template <typename UintTy>
            static constexpr void addCfX4Basic(const UintTy in[4], UintTy sum[4], bool &carry)
            {
                sum[0] = add_carry(in[0], carry);
                sum[1] = add_carry(in[1], carry);
                sum[2] = add_carry(in[2], carry);
                sum[3] = add_carry(in[3], carry);
            }
            template <typename UintTy>
            static constexpr void addCfX4(const UintTy in[4], UintTy sum[4], bool &carry)
            {
                addCfX4Basic(in, sum, carry);
            }
            static inline void addCfX4(const uint64_t in[4], uint64_t sum[4], bool &carry)
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
                addCfX4Basic(in, sum, carry);
#endif
            }
            template <typename UintTy>
            static constexpr void subBfX4Basic(const UintTy in[4], UintTy diff[4], bool &borrow)
            {
                diff[0] = sub_borrow(in[0], borrow);
                diff[1] = sub_borrow(in[1], borrow);
                diff[2] = sub_borrow(in[2], borrow);
                diff[3] = sub_borrow(in[3], borrow);
            }
            template <typename UintTy>
            static constexpr void subBfX4(const UintTy in[4], UintTy diff[4], bool &borrow)
            {
                subBfX4Basic(in, diff, borrow);
            }
            static inline void subBfX4(const uint64_t in[4], uint64_t diff[4], bool &borrow)
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
                subBfX4Basic(in, diff, borrow);
#endif
            }

            template <typename UintTy>
            static void mulInBase(UintTy a, UintTy b, UintTy &lo, UintTy &hi)
            {
                constexpr int NUM_BITS = sizeof(UintTy) * CHAR_BIT;
                static_assert(NUM_BITS <= 32, "Only for 32bit or lower");
                using ProdTy = typename Uint<NUM_BITS>::NextType::Type;
                ProdTy prod = ProdTy(a) * b;
                hi = UintTy(prod >> NUM_BITS);
                lo = UintTy(prod);
            }

            static void mulInBase(uint64_t a, uint64_t b, uint64_t &lo, uint64_t &hi)
            {
                mul64x64to128(a, b, lo, hi);
            }

            template <typename UintTy>
            static void mulAdd(UintTy a, UintTy b, UintTy c, UintTy &lo, UintTy &hi)
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
            static void mulAdd(uint64_t a, uint64_t b, uint64_t c, uint64_t &lo, uint64_t &hi)
            {
                mul64x64to128(a, b, lo, hi);
                lo += c;
                hi += (lo < c);
            }

            template <typename UintTy>
            static UintTy mulAddX4Basic(const UintTy in[4], UintTy num_mul, UintTy num_add, UintTy out[4])
            {
                mulAdd(in[0], num_mul, num_add, out[0], num_add);
                mulAdd(in[1], num_mul, num_add, out[1], num_add);
                mulAdd(in[2], num_mul, num_add, out[2], num_add);
                mulAdd(in[3], num_mul, num_add, out[3], num_add);
                return num_add;
            }
            template <typename UintTy>
            static UintTy mulAddX4(const UintTy in[4], UintTy num_mul, UintTy num_add, UintTy out[4])
            {
                return mulAddX4Basic(in, num_mul, num_add, out);
            }

            static uint64_t mulAddX4(const uint64_t in[4], uint64_t num_mul, uint64_t num_add, uint64_t out[4])
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
                return mulAddX4Basic(in, num_mul, num_add, out);
#endif
            }

            template <typename UintTy>
            static void mulAddPlusSelf(UintTy in, UintTy num_mul, UintTy &out, UintTy &carry)
            {
                bool cf;
                UintTy hi, lo;
                mulInBase(in, num_mul, lo, hi);
                lo = addHalf(lo, out, cf);
                hi += cf;
                out = addHalf(lo, carry, cf);
                carry = hi + cf;
            }
            template <typename UintTy>
            static UintTy addEqMulAddX8Basic(const UintTy in[], size_t len, UintTy out[], UintTy num_mul, UintTy num_add)
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
            template <typename UintTy>
            static UintTy addEqMulAddX8(const UintTy in[], size_t len, UintTy out[], UintTy num_mul, UintTy num_add)
            {
                return addEqMulAddX8Basic(in, len, out, num_mul, num_add);
            }

            static uint64_t addEqMulAddX8(const uint64_t in[], size_t len, uint64_t out[], uint64_t num_mul, uint64_t num_add)
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
                return addEqMulAddX8Basic(in, len, out, num_mul, num_add);
#endif
            }

            // // return dividend % 2^64
            // static auto divRemBase(Uint192 &dividend)
            // {
            //     uint64_t lo = dividend;
            //     dividend = dividend.rShift64();
            // }
        };

        // class Uint128
        // {
        // private:
        //     uint64_t low, high;

        // public:
        //     constexpr Uint128() : low(0), high(0) {}
        //     constexpr Uint128(uint64_t l, uint64_t h = 0) : low(l), high(h) {}
        //     constexpr Uint128(std::pair<uint64_t, uint64_t> p) : low(p.first), high(p.second) {}

        //     constexpr Uint128 operator+(Uint128 rhs) const
        //     {
        //         rhs.low += low;
        //         rhs.high += high + (rhs.low < low);
        //         return rhs;
        //     }
        //     constexpr Uint128 operator-(Uint128 rhs) const
        //     {
        //         rhs.low = low - rhs.low;
        //         rhs.high = high - rhs.high - (rhs.low > low);
        //         return rhs;
        //     }
        //     constexpr Uint128 operator+(uint64_t rhs) const
        //     {
        //         rhs = low + rhs;
        //         return Uint128(rhs, high + (rhs < low));
        //     }
        //     constexpr Uint128 operator-(uint64_t rhs) const
        //     {
        //         rhs = low - rhs;
        //         return Uint128(rhs, high - (rhs > low));
        //     }
        //     // Only compute the low * rhs.low
        //     Uint128 operator*(const Uint128 &rhs) const
        //     {
        //         Uint128 res;
        //         mul64x64to128(low, rhs.low, res.low, res.high);
        //         return res;
        //     }
        //     // Only compute the low * rhs
        //     Uint128 operator*(uint64_t rhs) const
        //     {
        //         Uint128 res;
        //         mul64x64to128(low, rhs, res.low, res.high);
        //         return res;
        //     }
        //     // Only compute the 128bit / 64 bit
        //     constexpr Uint128 operator/(const Uint128 &rhs) const
        //     {
        //         return *this / rhs.low;
        //     }
        //     // Only compute the 128bit % 64 bit
        //     constexpr Uint128 operator%(const Uint128 &rhs) const
        //     {
        //         return *this % rhs.low;
        //     }
        //     // Only compute the 128bit / 64 bit
        //     constexpr Uint128 operator/(uint64_t rhs) const
        //     {
        //         Uint128 quot = *this;
        //         quot.selfDivRem(rhs);
        //         return quot;
        //     }
        //     // Only compute the 128bit % 64 bit
        //     constexpr Uint128 operator%(uint64_t rhs) const
        //     {
        //         Uint128 quot = *this;
        //         uint64_t rem = quot.selfDivRem(rhs);
        //         return Uint128(rem);
        //     }
        //     constexpr Uint128 &operator+=(const Uint128 &rhs)
        //     {
        //         return *this = *this + rhs;
        //     }
        //     constexpr Uint128 &operator-=(const Uint128 &rhs)
        //     {
        //         return *this = *this - rhs;
        //     }
        //     constexpr Uint128 &operator+=(uint64_t rhs)
        //     {
        //         return *this = *this + rhs;
        //     }
        //     constexpr Uint128 &operator-=(uint64_t rhs)
        //     {
        //         return *this = *this - rhs;
        //     }
        //     // Only compute the low * rhs.low
        //     constexpr Uint128 &operator*=(const Uint128 &rhs)
        //     {
        //         mul64x64to128_base(low, rhs.low, low, high);
        //         return *this;
        //     }
        //     constexpr Uint128 &operator/=(const Uint128 &rhs)
        //     {
        //         return *this = *this / rhs;
        //     }
        //     constexpr Uint128 &operator%=(const Uint128 &rhs)
        //     {
        //         return *this = *this % rhs;
        //     }
        //     // Return *this % divisor, *this /= divisor
        //     constexpr uint64_t selfDivRem(uint64_t divisor)
        //     {
        //         if ((divisor >> 32) == 0)
        //         {
        //             return div128by32(high, low, uint32_t(divisor));
        //         }
        //         uint64_t divid1 = high % divisor, divid0 = low;
        //         high /= divisor;
        //         low = div128by64to64(divid1, divid0, divisor);
        //         return divid0;
        //     }
        //     static constexpr Uint128 mul64x64(uint64_t a, uint64_t b)
        //     {
        //         Uint128 res;
        //         mul64x64to128_base(a, b, res.low, res.high);
        //         return res;
        //     }
        //     friend constexpr bool operator<(const Uint128 &lhs, const Uint128 &rhs)
        //     {
        //         if (lhs.high != rhs.high)
        //         {
        //             return lhs.high < rhs.high;
        //         }
        //         return lhs.low < rhs.low;
        //     }
        //     friend constexpr bool operator<=(const Uint128 &lhs, const Uint128 &rhs)
        //     {
        //         return !(rhs > lhs);
        //     }
        //     friend constexpr bool operator>(const Uint128 &lhs, const Uint128 &rhs)
        //     {
        //         return rhs < lhs;
        //     }
        //     friend constexpr bool operator>=(const Uint128 &lhs, const Uint128 &rhs)
        //     {
        //         return !(lhs < rhs);
        //     }
        //     friend constexpr bool operator==(const Uint128 &lhs, const Uint128 &rhs)
        //     {
        //         return lhs.low == rhs.low && lhs.high == rhs.high;
        //     }
        //     friend constexpr bool operator!=(const Uint128 &lhs, const Uint128 &rhs)
        //     {
        //         return !(lhs == rhs);
        //     }

        //     constexpr Uint128 operator<<(int shift) const
        //     {
        //         if (shift == 0)
        //         {
        //             return *this;
        //         }
        //         shift %= 128;
        //         shift = shift < 0 ? shift + 128 : shift;
        //         if (shift < 64)
        //         {
        //             return Uint128(low << shift, (high << shift) | (low >> (64 - shift)));
        //         }
        //         return Uint128(0, low << (shift - 64));
        //     }
        //     constexpr Uint128 operator>>(int shift) const
        //     {
        //         if (shift == 0)
        //         {
        //             return *this;
        //         }
        //         shift %= 128;
        //         shift = shift < 0 ? shift + 128 : shift;
        //         if (shift < 64)
        //         {
        //             return Uint128((low >> shift) | (high << (64 - shift)), high >> shift);
        //         }
        //         return Uint128(high >> (shift - 64), 0);
        //     }
        //     constexpr Uint128 &operator<<=(int shift)
        //     {
        //         return *this = *this << shift;
        //     }
        //     constexpr Uint128 &operator>>=(int shift)
        //     {
        //         return *this = *this >> shift;
        //     }
        //     constexpr uint64_t high64() const
        //     {
        //         return high;
        //     }
        //     constexpr uint64_t low64() const
        //     {
        //         return low;
        //     }
        //     constexpr explicit operator uint64_t() const
        //     {
        //         return low64();
        //     }
        //     std::string toStringBase10() const
        //     {
        //         if (high == 0)
        //         {
        //             return std::to_string(low);
        //         }
        //         constexpr uint64_t BASE(10000'0000'0000'0000);
        //         Uint128 copy(*this);
        //         std::string s;
        //         s = ui64to_string_base10(uint64_t(copy.selfDivRem(BASE)), 16) + s;
        //         s = ui64to_string_base10(uint64_t(copy.selfDivRem(BASE)), 16) + s;
        //         return std::to_string(uint64_t(copy.selfDivRem(BASE))) + s;
        //     }
        //     void printDec() const
        //     {
        //         std::cout << std::dec << toStringBase10() << '\n';
        //     }
        //     void printHex() const
        //     {
        //         std::cout << std::hex << "0x" << high << ' ' << low << std::dec << '\n';
        //     }
        // };

        // class Uint192
        // {
        //     friend Uint128;

        // private:
        //     uint64_t low, mid, high;

        // public:
        //     constexpr Uint192() : low(0), mid(0), high(0) {}
        //     constexpr Uint192(uint64_t low, uint64_t mi = 0, uint64_t high = 0) : low(low), mid(mi), high(high) {}
        //     constexpr Uint192(Uint128 n) : low(n.low64()), mid(n.high64()), high(0) {}
        //     constexpr Uint192 operator+(Uint192 rhs) const
        //     {
        //         bool cf = false;
        //         rhs.low = add_half(low, rhs.low, cf);
        //         rhs.mid = add_carry(mid, rhs.mid, cf);
        //         rhs.high = high + rhs.high + cf;
        //         return rhs;
        //     }
        //     constexpr Uint192 operator-(Uint192 rhs) const
        //     {
        //         bool bf = false;
        //         rhs.low = sub_half(low, rhs.low, bf);
        //         rhs.mid = sub_borrow(mid, rhs.mid, bf);
        //         rhs.high = high - rhs.high - bf;
        //         return rhs;
        //     }
        //     constexpr Uint192 operator/(uint64_t rhs) const
        //     {
        //         Uint192 result(*this);
        //         result.selfDivRem(rhs);
        //         return result;
        //     }
        //     constexpr Uint192 operator%(uint64_t rhs) const
        //     {
        //         Uint192 result(*this);
        //         return result.selfDivRem(rhs);
        //     }
        //     constexpr Uint192 &operator+=(const Uint192 &rhs)
        //     {
        //         return *this = *this + rhs;
        //     }
        //     constexpr Uint192 &operator-=(const Uint192 &rhs)
        //     {
        //         return *this = *this - rhs;
        //     }
        //     constexpr Uint192 &operator/=(const Uint192 &rhs)
        //     {
        //         return *this = *this / rhs;
        //     }
        //     constexpr Uint192 &operator%=(const Uint192 &rhs)
        //     {
        //         return *this = *this % rhs;
        //     }
        //     constexpr Uint192 operator<<(int shift) const
        //     {
        //         if (shift == 0)
        //         {
        //             return *this;
        //         }
        //         shift %= 192;
        //         shift = shift < 0 ? shift + 192 : shift;
        //         if (shift < 64)
        //         {
        //             return Uint192(low << shift, (mid << shift) | (low >> (64 - shift)), (high << shift) | (mid >> (64 - shift)));
        //         }
        //         else if (shift < 128)
        //         {
        //             shift -= 64;
        //             return Uint192(0, low << shift, (mid << shift) | (low >> (64 - shift)));
        //         }
        //         return Uint192(0, 0, low << (shift - 128));
        //     }
        //     friend constexpr bool operator<(const Uint192 &lhs, const Uint192 &rhs)
        //     {
        //         if (lhs.high != rhs.high)
        //         {
        //             return lhs.high < rhs.high;
        //         }
        //         if (lhs.mid != rhs.mid)
        //         {
        //             return lhs.mid < rhs.mid;
        //         }
        //         return lhs.low < rhs.low;
        //     }
        //     friend constexpr bool operator<=(const Uint192 &lhs, const Uint192 &rhs)
        //     {
        //         return !(rhs > lhs);
        //     }
        //     friend constexpr bool operator>(const Uint192 &lhs, const Uint128 &rhs)
        //     {
        //         return rhs < lhs;
        //     }
        //     friend constexpr bool operator>=(const Uint192 &lhs, const Uint192 &rhs)
        //     {
        //         return !(lhs < rhs);
        //     }
        //     friend constexpr bool operator==(const Uint192 &lhs, const Uint192 &rhs)
        //     {
        //         return lhs.low == rhs.low && lhs.mid == rhs.mid && lhs.high == rhs.high;
        //     }
        //     friend constexpr bool operator!=(const Uint192 &lhs, const Uint192 &rhs)
        //     {
        //         return !(lhs == rhs);
        //     }
        //     static constexpr Uint192 mul128x64(Uint128 a, uint64_t b)
        //     {
        //         auto prod1 = Uint128::mul64x64(b, a.low64());
        //         auto prod2 = Uint128::mul64x64(b, a.high64());
        //         Uint192 result;
        //         result.low = prod1.low64();
        //         result.mid = prod1.high64() + prod2.low64();
        //         result.high = prod2.high64() + (result.mid < prod1.high64());
        //         return result;
        //     }
        //     static constexpr Uint192 mul64x64x64(uint64_t a, uint64_t b, uint64_t c)
        //     {
        //         return mul128x64(Uint128::mul64x64(a, b), c);
        //     }
        //     constexpr uint64_t selfDivRem(uint64_t divisor)
        //     {
        //         uint64_t divid1 = high % divisor, divid0 = mid;
        //         high /= divisor;
        //         mid = div128by64to64(divid1, divid0, divisor);
        //         divid1 = divid0, divid0 = low;
        //         low = div128by64to64(divid1, divid0, divisor);
        //         return divid0;
        //     }
        //     constexpr Uint192 rShift64() const
        //     {
        //         return Uint192(mid, high, 0);
        //     }
        //     constexpr operator uint64_t() const
        //     {
        //         return low;
        //     }
        //     std::string toStringBase10() const
        //     {
        //         if (high == 0)
        //         {
        //             return Uint128(mid, low).toStringBase10();
        //         }
        //         constexpr uint64_t BASE(10000'0000'0000'0000);
        //         Uint192 copy(*this);
        //         std::string s;
        //         s = ui64to_string_base10(uint64_t(copy.selfDivRem(BASE)), 16) + s;
        //         s = ui64to_string_base10(uint64_t(copy.selfDivRem(BASE)), 16) + s;
        //         s = ui64to_string_base10(uint64_t(copy.selfDivRem(BASE)), 16) + s;
        //         return std::to_string(uint64_t(copy.selfDivRem(BASE))) + s;
        //     }
        //     void printDec() const
        //     {
        //         std::cout << std::dec << toStringBase10() << '\n';
        //     }
        //     void printHex() const
        //     {
        //         std::cout << std::hex << "0x" << high << ' ' << mid << ' ' << low << std::dec << '\n';
        //     }
        // };

        // template <typename Int128Type>
        // constexpr uint64_t high64(const Int128Type &n)
        // {
        //     return n >> 64;
        // }
        // constexpr uint64_t high64(const Uint128 &n)
        // {
        //     return n.high64();
        // }
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
        namespace ntt
        {

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
                constexpr utility::BaseExecutorBinary exec;
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
                constexpr utility::BaseExecutorBinary exec;
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
        }
        namespace division
        {

        }
        namespace base_conversion
        {

        }
    }
}
#endif