#pragma once
#ifndef HINT_BASIC_HPP
#define HINT_BASIC_HPP

#include <vector>
#include <complex>
#include <iostream>
#include <future>
#include <string>
#include <array>
#include <type_traits>
#include <random>
#include <utility>
#include <fstream>
#include <cassert>
#include <climits>
#include <cstdint>

// Windows 64bit fast multiply macro.
#if defined(_WIN64)
#include <intrin.h>
#define UMUL128
#endif //_WIN64

// GCC 64bit fast multiply macro.
#if defined(__SIZEOF_INT128__)
#define UINT128T
#endif //__SIZEOF_INT128__

namespace hint
{
    // bits of 1, equals to 2^bits - 1
    template <typename T>
    constexpr T all_one(int bits)
    {
        T temp = T(1) << (bits - 1);
        return temp - 1 + temp;
    }

    // return x == 2^k
    template <typename IntTy>
    constexpr bool is_2pow(IntTy x)
    {
        if (x <= 0)
        {
            return false;
        }
        return (x & (x - 1)) == 0;
    }

    // Leading zeros
    template <typename IntTy>
    constexpr int hint_clz(IntTy x)
    {
        constexpr uint32_t MASK32 = uint32_t(0xFFFF) << 16;
        int res = sizeof(IntTy) * CHAR_BIT;
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

    // Integer bit length
    template <typename IntTy>
    constexpr int hint_bit_length(IntTy x)
    {
        if (x == 0)
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

    // Fast power
    template <typename T, typename T1>
    constexpr T qpow(T m, T1 n)
    {
        T result = 1;
        while (n > 0)
        {
            if ((n & 1) != 0)
            {
                result *= m;
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
        while (n > 0)
        {
            if ((n & 1) != 0)
            {
                result *= m;
                result %= mod;
            }
            m *= m;
            m %= mod;
            n >>= 1;
        }
        return result;
    }

    // Get cloest power of 2 that not larger than n
    template <typename T>
    constexpr T int_floor2(T n)
    {
        constexpr int bits = sizeof(n) * CHAR_BIT;
        for (int i = 1; i < bits; i *= 2)
        {
            n |= (n >> i);
        }
        return (n >> 1) + 1;
    }

    // Get cloest power of 2 that not smaller than n
    template <typename T>
    constexpr T int_ceil2(T n)
    {
        constexpr int bits = sizeof(n) * CHAR_BIT;
        n--;
        for (int i = 1; i < bits; i *= 2)
        {
            n |= (n >> i);
        }
        return n + 1;
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

    // a * x + b * y = gcd(a,b)
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
        // Compute Integer multiplication, 64bit x 64bit to 128bit, basic algorithm
        // first is low 64bit, second is high 64bit
        constexpr void mul64x64to128_base(uint64_t a, uint64_t b, uint64_t &low, uint64_t &high)
        {
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

        inline void mul64x64to128(uint64_t a, uint64_t b, uint64_t &low, uint64_t &high)
        {
#if defined(UMUL128)
#pragma message("Using _umul128 to compute 64bit x 64bit to 128bit")
            unsigned long long lo, hi;
            lo = _umul128(a, b, &hi);
            low = lo, high = hi;
#else
#if defined(UINT128T) // No _umul128
#pragma message("Using __uint128_t to compute 64bit x 64bit to 128bit")
            __uint128_t x(a);
            x *= b;
            low = uint64_t(x), high = uint64_t(x >> 64);
#else // No __uint128_t
#pragma message("Using basic function to compute 64bit x 64bit to 128bit")
            mul64x64to128_base(a, b, low, high);
#endif // UINT128T
#endif // UMUL128
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

        class Uint128
        {
        private:
            uint64_t low, high;

        public:
            constexpr Uint128() : low(0), high(0) {}
            constexpr Uint128(uint64_t l, uint64_t h = 0) : low(l), high(h) {}
            constexpr Uint128(std::pair<uint64_t, uint64_t> p) : low(p.first), high(p.second) {}

            constexpr Uint128 operator+(Uint128 rhs) const
            {
                rhs.low += low;
                rhs.high += high + (rhs.low < low);
                return rhs;
            }
            constexpr Uint128 operator-(Uint128 rhs) const
            {
                rhs.low = low - rhs.low;
                rhs.high = high - rhs.high - (rhs.low > low);
                return rhs;
            }
            constexpr Uint128 operator+(uint64_t rhs) const
            {
                rhs = low + rhs;
                return Uint128(rhs, high + (rhs < low));
            }
            constexpr Uint128 operator-(uint64_t rhs) const
            {
                rhs = low - rhs;
                return Uint128(rhs, high - (rhs > low));
            }
            // Only compute the low * rhs.low
            Uint128 operator*(const Uint128 &rhs) const
            {
                Uint128 res;
                mul64x64to128(low, rhs.low, res.low, res.high);
                return res;
            }
            // Only compute the low * rhs
            Uint128 operator*(uint64_t rhs) const
            {
                Uint128 res;
                mul64x64to128(low, rhs, res.low, res.high);
                return res;
            }
            // Only compute the 128bit / 64 bit
            constexpr Uint128 operator/(const Uint128 &rhs) const
            {
                return *this / rhs.low;
            }
            // Only compute the 128bit % 64 bit
            constexpr Uint128 operator%(const Uint128 &rhs) const
            {
                return *this % rhs.low;
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
                Uint128 res;
                mul64x64to128_base(a, b, res.low, res.high);
                return res;
            }
            friend constexpr bool operator<(const Uint128 &lhs, const Uint128 &rhs)
            {
                if (lhs.high != rhs.high)
                {
                    return lhs.high < rhs.high;
                }
                return lhs.low < rhs.low;
            }
            friend constexpr bool operator<=(const Uint128 &lhs, const Uint128 &rhs)
            {
                return !(rhs > lhs);
            }
            friend constexpr bool operator>(const Uint128 &lhs, const Uint128 &rhs)
            {
                return rhs < lhs;
            }
            friend constexpr bool operator>=(const Uint128 &lhs, const Uint128 &rhs)
            {
                return !(lhs < rhs);
            }
            friend constexpr bool operator==(const Uint128 &lhs, const Uint128 &rhs)
            {
                return lhs.low == rhs.low && lhs.high == rhs.high;
            }
            friend constexpr bool operator!=(const Uint128 &lhs, const Uint128 &rhs)
            {
                return !(lhs == rhs);
            }

            constexpr Uint128 operator<<(int shift) const
            {
                if (shift == 0)
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
                if (shift == 0)
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
            constexpr uint64_t high64() const
            {
                return high;
            }
            constexpr uint64_t low64() const
            {
                return low;
            }
            constexpr explicit operator uint64_t() const
            {
                return low64();
            }
            std::string toStringBase10() const
            {
                if (high == 0)
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
        };

        class Uint192
        {
            friend Uint128;

        private:
            uint64_t low, mid, high;

        public:
            constexpr Uint192() : low(0), mid(0), high(0) {}
            constexpr Uint192(uint64_t low, uint64_t mi = 0, uint64_t high = 0) : low(low), mid(mi), high(high) {}
            constexpr Uint192(Uint128 n) : low(n.low64()), mid(n.high64()), high(0) {}
            constexpr Uint192 operator+(Uint192 rhs) const
            {
                bool cf = false;
                rhs.low = add_half(low, rhs.low, cf);
                rhs.mid = add_carry(mid, rhs.mid, cf);
                rhs.high = high + rhs.high + cf;
                return rhs;
            }
            constexpr Uint192 operator-(Uint192 rhs) const
            {
                bool bf = false;
                rhs.low = sub_half(low, rhs.low, bf);
                rhs.mid = sub_borrow(mid, rhs.mid, bf);
                rhs.high = high - rhs.high - bf;
                return rhs;
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
            constexpr Uint192 &operator+=(const Uint192 &rhs)
            {
                return *this = *this + rhs;
            }
            constexpr Uint192 &operator-=(const Uint192 &rhs)
            {
                return *this = *this - rhs;
            }
            constexpr Uint192 &operator/=(const Uint192 &rhs)
            {
                return *this = *this / rhs;
            }
            constexpr Uint192 &operator%=(const Uint192 &rhs)
            {
                return *this = *this % rhs;
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
            friend constexpr bool operator<(const Uint192 &lhs, const Uint192 &rhs)
            {
                if (lhs.high != rhs.high)
                {
                    return lhs.high < rhs.high;
                }
                if (lhs.mid != rhs.mid)
                {
                    return lhs.mid < rhs.mid;
                }
                return lhs.low < rhs.low;
            }
            friend constexpr bool operator<=(const Uint192 &lhs, const Uint192 &rhs)
            {
                return !(rhs > lhs);
            }
            friend constexpr bool operator>(const Uint192 &lhs, const Uint128 &rhs)
            {
                return rhs < lhs;
            }
            friend constexpr bool operator>=(const Uint192 &lhs, const Uint192 &rhs)
            {
                return !(lhs < rhs);
            }
            friend constexpr bool operator==(const Uint192 &lhs, const Uint192 &rhs)
            {
                return lhs.low == rhs.low && lhs.mid == rhs.mid && lhs.high == rhs.high;
            }
            friend constexpr bool operator!=(const Uint192 &lhs, const Uint192 &rhs)
            {
                return !(lhs == rhs);
            }
            static constexpr Uint192 mul128x64(Uint128 a, uint64_t b)
            {
                auto prod1 = Uint128::mul64x64(b, a.low64());
                auto prod2 = Uint128::mul64x64(b, a.high64());
                Uint192 result;
                result.low = prod1.low64();
                result.mid = prod1.high64() + prod2.low64();
                result.high = prod2.high64() + (result.mid < prod1.high64());
                return result;
            }
            static constexpr Uint192 mul64x64x64(uint64_t a, uint64_t b, uint64_t c)
            {
                return mul128x64(Uint128::mul64x64(a, b), c);
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
        };

        template <typename Int128Type>
        constexpr uint64_t high64(const Int128Type &n)
        {
            return n >> 64;
        }
        constexpr uint64_t high64(const Uint128 &n)
        {
            return n.high64();
        }
#ifdef UINT128T
        using Uint128Default = __uint128_t;
#else
        using Uint128Default = Uint128;
#endif // UINT128T

        template <int BITS>
        struct Uint
        {
            static constexpr int MAX_BITS = 128;
            static constexpr int MIN_BITS = 8;
            static constexpr int NORM_BITS = hint::int_ceil2(std::min(MAX_BITS, std::max(MIN_BITS, BITS)));
            using SignType = typename Uint<NORM_BITS>::SignType;
            using Type = typename Uint<NORM_BITS>::Type;
        };

        template <int BITS>
        constexpr int Uint<BITS>::MAX_BITS;
        template <int BITS>
        constexpr int Uint<BITS>::MIN_BITS;
        template <int BITS>
        constexpr int Uint<BITS>::NORM_BITS;

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

        template <>
        struct Uint<128>
        {
            using SignType = Uint128Default;
            using Type = Uint128Default;
        };

        template <>
        struct Uint<192>
        {
            using Type = Uint192;
        };

        template <int BITS>
        struct UintType
        {
            static constexpr int MAX_BITS = 128;
            static constexpr int MIN_BITS = 8;
            static constexpr int NORM_BITS = hint::int_ceil2(std::min(MAX_BITS, std::max(MIN_BITS, BITS)));
            static constexpr int NEXT_BITS = std::min(MAX_BITS, NORM_BITS * 2);
            static constexpr int LAST_BITS = std::max(MIN_BITS, NORM_BITS / 2);

            using SignType = typename Uint<NORM_BITS>::SignType;
            using Type = typename Uint<NORM_BITS>::Type;

            using NextType = UintType<NEXT_BITS>;
            using LastType = UintType<LAST_BITS>;
        };
        template <int BITS>
        constexpr int UintType<BITS>::MAX_BITS;
        template <int BITS>
        constexpr int UintType<BITS>::MIN_BITS;
        template <int BITS>
        constexpr int UintType<BITS>::NORM_BITS;
        template <int BITS>
        constexpr int UintType<BITS>::NEXT_BITS;
        template <int BITS>
        constexpr int UintType<BITS>::LAST_BITS;
    }

    namespace mod_int
    {
        using extend_int::Uint128;

        //  Montgomery for mod > 2^32
        //  default R = 2^64
        template <uint64_t MOD, typename Int128Type = Uint128>
        class MontInt64Lazy
        {
        private:
            static_assert(MOD > UINT32_MAX, "Montgomery64 modulus must be greater than 2^32");
            static_assert(hint_log2(MOD) < 62, "MOD can't be larger than 62 bits");
            uint64_t data;

        public:
            using IntType = uint64_t;

            MontInt64Lazy() : data(0) {}
            constexpr MontInt64Lazy(uint64_t n) : data(mulMontCompileTime(n, rSquare())) {}

            constexpr MontInt64Lazy operator+(MontInt64Lazy rhs) const
            {
                rhs.data = data + rhs.data;
                rhs.data = rhs.data < mod2() ? rhs.data : rhs.data - mod2();
                return rhs;
            }
            constexpr MontInt64Lazy operator-(MontInt64Lazy rhs) const
            {
                rhs.data = data - rhs.data;
                rhs.data = rhs.data > data ? rhs.data + mod2() : rhs.data;
                return rhs;
            }
            MontInt64Lazy operator*(MontInt64Lazy rhs) const
            {
                rhs.data = mulMontRunTimeLazy(data, rhs.data);
                return rhs;
            }
            constexpr MontInt64Lazy &operator+=(const MontInt64Lazy &rhs)
            {
                return *this = *this + rhs;
            }
            constexpr MontInt64Lazy &operator-=(const MontInt64Lazy &rhs)
            {
                return *this = *this - rhs;
            }
            constexpr MontInt64Lazy &operator*=(const MontInt64Lazy &rhs)
            {
                data = mulMontCompileTime(data, rhs.data);
                return *this;
            }
            constexpr MontInt64Lazy largeNorm() const
            {
                MontInt64Lazy res;
                res.data = data >= mod2() ? data - mod2() : data;
                return res;
            }
            constexpr MontInt64Lazy add(MontInt64Lazy rhs) const
            {
                rhs.data = data + rhs.data;
                return rhs;
            }
            constexpr MontInt64Lazy sub(MontInt64Lazy rhs) const
            {
                rhs.data = data - rhs.data + mod2();
                return rhs;
            }
            constexpr operator uint64_t() const
            {
                return toInt(data);
            }

            static constexpr uint64_t mod()
            {
                return MOD;
            }
            static constexpr uint64_t mod2()
            {
                return MOD * 2;
            }
            static constexpr uint64_t modInv()
            {
                constexpr uint64_t mod_inv = inv_mod2pow(mod(), 64); //(mod_inv * mod)%(2^64) = 1
                return mod_inv;
            }
            static constexpr uint64_t modInvNeg()
            {
                constexpr uint64_t mod_inv_neg = uint64_t(0 - modInv()); //(mod_inv_neg + mod_inv)%(2^64) = 0
                return mod_inv_neg;
            }
            static constexpr uint64_t rSquare()
            {
                constexpr Int128Type r = (Int128Type(1) << 64) % Int128Type(mod()); // R % mod
                constexpr uint64_t r2 = uint64_t(qpow(r, 2, Int128Type(mod())));    // R^2 % mod
                return r2;
            }
            static_assert((mod() * modInv()) == 1, "mod_inv not correct");

            static constexpr uint64_t toMont(uint64_t n)
            {
                return mulMontCompileTime(n, rSquare());
            }
            static constexpr uint64_t toInt(uint64_t n)
            {
                return redc(Int128Type(n));
            }

            static uint64_t redcFastLazy(const Int128Type &input)
            {
                Int128Type n = uint64_t(input) * modInvNeg();
                n = n * mod();
                n += input;
                return hint::extend_int::high64(n);
            }
            static uint64_t redcFast(const Int128Type &input)
            {
                uint64_t n = redcFastLazy(input);
                return n < mod() ? n : n - mod();
            }
            static constexpr uint64_t redc(const Int128Type &input)
            {
                Int128Type n = uint64_t(input) * modInvNeg();
                n *= Int128Type(mod());
                n += input;
                uint64_t m = hint::extend_int::high64(n);
                return m < mod() ? m : m - mod();
            }
            static uint64_t mulMontRunTime(uint64_t a, uint64_t b)
            {
                return redcFast(Int128Type(a) * b);
            }
            static uint64_t mulMontRunTimeLazy(uint64_t a, uint64_t b)
            {
                return redcFastLazy(Int128Type(a) * b);
            }
            static constexpr uint64_t mulMontCompileTime(uint64_t a, uint64_t b)
            {
                Int128Type prod(a);
                prod *= Int128Type(b);
                return redc(prod);
            }
        };

        //  Montgomery for mod < 2^30
        //  default R = 2^32
        template <uint32_t MOD>
        class MontInt32Lazy
        {
        private:
            static_assert(hint_log2(MOD) < 30, "MOD can't be larger than 30 bits");
            uint32_t data;

        public:
            using IntType = uint32_t;

            MontInt32Lazy() : data(0) {}
            constexpr MontInt32Lazy(uint32_t n) : data(toMont(n)) {}

            constexpr MontInt32Lazy operator+(MontInt32Lazy rhs) const
            {
                rhs.data = data + rhs.data;
                rhs.data = rhs.data < mod2() ? rhs.data : rhs.data - mod2();
                return rhs;
            }
            constexpr MontInt32Lazy operator-(MontInt32Lazy rhs) const
            {
                rhs.data = data - rhs.data;
                rhs.data = rhs.data > data ? rhs.data + mod2() : rhs.data;
                return rhs;
            }
            constexpr MontInt32Lazy operator*(MontInt32Lazy rhs) const
            {
                rhs.data = redcLazy(uint64_t(data) * rhs.data);
                return rhs;
            }
            constexpr MontInt32Lazy &operator+=(const MontInt32Lazy &rhs)
            {
                return *this = *this + rhs;
            }
            constexpr MontInt32Lazy &operator-=(const MontInt32Lazy &rhs)
            {
                return *this = *this - rhs;
            }
            constexpr MontInt32Lazy &operator*=(const MontInt32Lazy &rhs)
            {
                data = redc(uint64_t(data) * rhs.data);
                return *this;
            }
            constexpr MontInt32Lazy largeNorm() const
            {
                MontInt32Lazy res;
                res.data = data >= mod2() ? data - mod2() : data;
                return res;
            }
            constexpr MontInt32Lazy add(MontInt32Lazy rhs) const
            {
                rhs.data = data + rhs.data;
                return rhs;
            }
            constexpr MontInt32Lazy sub(MontInt32Lazy rhs) const
            {
                rhs.data = data - rhs.data + mod2();
                return rhs;
            }
            constexpr operator uint32_t() const
            {
                return toInt(data);
            }

            static constexpr uint32_t mod()
            {
                return MOD;
            }
            static constexpr uint32_t mod2()
            {
                return MOD * 2;
            }
            static constexpr uint32_t modInv()
            {
                constexpr uint32_t mod_inv = uint32_t(inv_mod2pow(mod(), 32));
                return mod_inv;
            }
            static constexpr uint32_t modNegInv()
            {
                constexpr uint32_t mod_neg_inv = uint32_t(0 - modInv());
                return mod_neg_inv;
            }
            static_assert((mod() * modInv()) == 1, "mod_inv not correct");

            static constexpr uint32_t toMont(uint32_t n)
            {
                return (uint64_t(n) << 32) % MOD;
            }
            static constexpr uint32_t toInt(uint32_t n)
            {
                return redc(n);
            }

            static constexpr uint32_t redcLazy(uint64_t n)
            {
                uint32_t prod = uint32_t(n) * modNegInv();
                return (uint64_t(prod) * mod() + n) >> 32;
            }
            static constexpr uint32_t redc(uint64_t n)
            {
                uint32_t res = redcLazy(n);
                return res < mod() ? res : res - mod();
            }
        };

        //  ModInt for mod < 2^32
        template <uint32_t MOD>
        class ModInt32
        {
        private:
            uint32_t data;

        public:
            using IntType = uint32_t;

            constexpr ModInt32() : data(0) {}
            constexpr ModInt32(uint32_t in) : data(in) {}

            constexpr ModInt32 operator+(ModInt32 in) const
            {
                uint32_t diff = MOD - data;
                return in.data > diff ? in.data - diff : in.data + data;
            }
            constexpr ModInt32 operator-(ModInt32 in) const
            {
                in.data = data - in.data;
                return in.data > data ? in.data + MOD : in.data;
            }
            constexpr ModInt32 operator*(ModInt32 in) const
            {
                return mul64(in) % MOD;
            }
            constexpr ModInt32 &operator+=(ModInt32 in)
            {
                return *this = *this + in;
            }
            constexpr ModInt32 &operator-=(ModInt32 in)
            {
                return *this = *this - in;
            }
            constexpr ModInt32 &operator*=(ModInt32 in)
            {
                return *this = *this * in;
            }
            constexpr ModInt32 largeNorm() const
            {
                return data;
            }
            constexpr ModInt32 add(ModInt32 n) const
            {
                return *this + n;
            }
            constexpr ModInt32 sub(ModInt32 n) const
            {
                return *this - n;
            }
            constexpr uint64_t mul64(ModInt32 in) const
            {
                return uint64_t(data) * in.data;
            }
            constexpr operator uint32_t() const
            {
                return data;
            }
            static constexpr uint32_t mod()
            {
                return MOD;
            }
        };

    };

    namespace transform
    {
        template <typename T>
        inline void transform2(T &sum, T &diff)
        {
            T temp0 = sum, temp1 = diff;
            sum = temp0 + temp1;
            diff = temp0 - temp1;
        }

        // Multi mode, self checking, fast number theoretic transform.
        namespace ntt
        {
            constexpr uint64_t MOD0 = 2485986994308513793, ROOT0 = 5;
            constexpr uint64_t MOD1 = 1945555039024054273, ROOT1 = 5;
            constexpr uint64_t MOD2 = 4179340454199820289, ROOT2 = 3;
            constexpr uint64_t MOD3 = 754974721, ROOT3 = 11;
            constexpr uint64_t MOD4 = 469762049, ROOT4 = 3;
            constexpr uint64_t MOD5 = 3489660929, ROOT5 = 3;
            constexpr uint64_t MOD6 = 3221225473, ROOT6 = 5;

            using namespace extend_int;
            using namespace mod_int;

            template <typename IntType>
            constexpr bool check_inv(uint64_t n, uint64_t n_inv, uint64_t mod)
            {
                n %= mod;
                n_inv %= mod;
                IntType m(n);
                m *= IntType(n_inv);
                m %= IntType(mod);
                return m == IntType(1);
            }

            // 快速计算两模数的中国剩余定理
            template <uint32_t MOD1, uint32_t MOD2>
            inline uint64_t crt2(uint32_t num1, uint32_t num2)
            {
                constexpr uint64_t inv1 = mod_inv<int64_t>(MOD1, MOD2);
                constexpr uint64_t inv2 = mod_inv<int64_t>(MOD2, MOD1);
                static_assert(check_inv<uint64_t>(inv1, MOD1, MOD2), "Inv1 error");
                static_assert(check_inv<uint64_t>(inv2, MOD2, MOD1), "Inv2 error");
                if (num1 > num2)
                {
                    return (uint64_t(num1 - num2) * uint64_t(inv2) % MOD1) * MOD2 + num2;
                }
                else
                {
                    return (uint64_t(num2 - num1) * uint64_t(inv1) % MOD2) * MOD1 + num1;
                }
            }

            // 快速计算两模数的中国剩余定理
            template <uint64_t MOD1, uint64_t MOD2, typename Int128Type = Uint128>
            inline Int128Type crt2(uint64_t num1, uint64_t num2)
            {
                constexpr uint64_t inv1 = mod_inv<int64_t>(MOD1, MOD2);
                constexpr uint64_t inv2 = mod_inv<int64_t>(MOD2, MOD1);
                static_assert(check_inv<Int128Type>(inv1, MOD1, MOD2), "Inv1 error");
                static_assert(check_inv<Int128Type>(inv2, MOD2, MOD1), "Inv2 error");
                if (num1 > num2)
                {
                    return (Int128Type(num1 - num2) * Int128Type(inv2) % Int128Type(MOD1)) * Int128Type(MOD2) + num2;
                }
                else
                {
                    return (Int128Type(num2 - num1) * Int128Type(inv1) % Int128Type(MOD2)) * Int128Type(MOD1) + num1;
                }
            }

            // 快速计算两模数的中国剩余定理
            template <typename ModInt1, typename ModInt2>
            inline Uint128 crt2(ModInt1 num1, ModInt2 num2)
            {
                constexpr uint64_t MOD1 = ModInt1::mod();
                constexpr uint64_t MOD2 = ModInt2::mod();
                constexpr ModInt1 MOD2_INV1 = mod_inv<int64_t>(MOD2, MOD1);
                constexpr ModInt2 MOD1_INV2 = mod_inv<int64_t>(MOD1, MOD2);
                constexpr Uint128 MOD12 = Uint128::mul64x64(MOD1, MOD2);
                static_assert(check_inv<Uint128>(MOD1, MOD1_INV2, MOD2), "INV1 error");
                static_assert(check_inv<Uint128>(MOD2, MOD2_INV1, MOD1), "INV2 error");
                num1 = num1 * MOD2_INV1;
                num2 = num2 * MOD1_INV2;
                Uint128 result = Uint128(MOD2) * uint64_t(num1);
                result += Uint128(MOD1) * uint64_t(num2);
                return result < MOD12 ? result : result - MOD12;
            }

            // 快速计算三模数的中国剩余定理
            template <typename ModInt1, typename ModInt2, typename ModInt3>
            inline Uint192 crt3(ModInt1 n1, ModInt2 n2, ModInt3 n3)
            {
                constexpr uint64_t MOD1 = ModInt1::mod(), MOD2 = ModInt2::mod(), MOD3 = ModInt3::mod();
                constexpr Uint192 MOD123 = Uint192::mul64x64x64(MOD1, MOD2, MOD3);                                   // MOD1*MOD2*MOD3
                constexpr Uint128 MOD12 = Uint128::mul64x64(MOD1, MOD2);                                             // MOD1*MOD2
                constexpr Uint128 MOD23 = Uint128::mul64x64(MOD2, MOD3);                                             // MOD2*MOD3
                constexpr Uint128 MOD13 = Uint128::mul64x64(MOD1, MOD3);                                             // MOD1*MOD3
                constexpr uint64_t MOD23_M1 = uint64_t(Uint128::mul64x64(MOD2 % MOD1, MOD3 % MOD1) % Uint128(MOD1)); // (MOD2*MOD3)  mod MOD1
                constexpr uint64_t MOD13_M2 = uint64_t(Uint128::mul64x64(MOD1 % MOD2, MOD3 % MOD2) % Uint128(MOD2)); // (MOD1*MOD3)  mod MOD2
                constexpr uint64_t MOD12_M3 = uint64_t(Uint128::mul64x64(MOD1 % MOD3, MOD2 % MOD3) % Uint128(MOD3)); // (MOD1*MOD2)  mod MOD3
                constexpr ModInt1 MOD23_INV1 = mod_inv<int64_t>(MOD23_M1, MOD1);                                     // (MOD2*MOD3)^-1 mod MOD1
                constexpr ModInt2 MOD13_INV2 = mod_inv<int64_t>(MOD13_M2, MOD2);                                     // (MOD1*MOD3)^-1 mod MOD2
                constexpr ModInt3 MOD12_INV3 = mod_inv<int64_t>(MOD12_M3, MOD3);                                     // (MOD1*MOD2)^-1 mod MOD3
                static_assert(check_inv<Uint128>(MOD23_INV1, MOD23_M1, MOD1), "INV1 error");
                static_assert(check_inv<Uint128>(MOD13_INV2, MOD13_M2, MOD2), "INV2 error");
                static_assert(check_inv<Uint128>(MOD12_INV3, MOD12_M3, MOD3), "INV3 error");
                n1 = n1 * MOD23_INV1;
                n2 = n2 * MOD13_INV2;
                n3 = n3 * MOD12_INV3;
                Uint192 result = Uint192::mul128x64(MOD23, uint64_t(n1));
                result += Uint192::mul128x64(MOD13, uint64_t(n2));
                result += Uint192::mul128x64(MOD12, uint64_t(n3));
                result = result < MOD123 ? result : result - MOD123;
                return result < MOD123 ? result : result - MOD123;
            }
            namespace split_radix
            {
                template <uint64_t ROOT, typename ModIntType>
                inline ModIntType mul_w41(ModIntType n)
                {
                    constexpr ModIntType W_4_1 = qpow(ModIntType(ROOT), (ModIntType::mod() - 1) / 4);
                    return n * W_4_1;
                }

                // in: in_out0<4p, in_ou1<4p; in_out2<2p, in_ou3<2p
                // out: in_out0<4p, in_ou1<4p; in_out2<4p, in_ou3<4p
                template <uint64_t ROOT, typename ModIntType>
                inline void dit_butterfly244(ModIntType &in_out0, ModIntType &in_out1, ModIntType &in_out2, ModIntType &in_out3)
                {
                    ModIntType temp0, temp1, temp2, temp3;
                    temp0 = in_out0.largeNorm();
                    temp1 = in_out1.largeNorm();
                    temp2 = in_out2 + in_out3;
                    temp3 = in_out2.sub(in_out3);
                    temp3 = mul_w41<ROOT>(temp3);
                    in_out0 = temp0.add(temp2);
                    in_out2 = temp0.sub(temp2);
                    in_out1 = temp1.add(temp3);
                    in_out3 = temp1.sub(temp3);
                }

                // in: in_out0<2p, in_ou1<2p; in_out2<2p, in_ou3<2p
                // out: in_out0<2p, in_ou1<2p; in_out2<4p, in_ou3<4p
                template <uint64_t ROOT, typename ModIntType>
                inline void dif_butterfly244(ModIntType &in_out0, ModIntType &in_out1, ModIntType &in_out2, ModIntType &in_out3)
                {
                    ModIntType temp0, temp1, temp2, temp3;
                    temp0 = in_out0.add(in_out2);
                    temp2 = in_out0 - in_out2;
                    temp1 = in_out1.add(in_out3);
                    temp3 = in_out1.sub(in_out3);
                    temp3 = mul_w41<ROOT>(temp3);
                    in_out0 = temp0.largeNorm();
                    in_out1 = temp1.largeNorm();
                    in_out2 = temp2.add(temp3);
                    in_out3 = temp2.sub(temp3);
                }

                // in: in_out0<4p, in_ou1<4p
                // out: in_out0<4p, in_ou1<4p
                template <typename ModIntType>
                inline void dit_butterfly2(ModIntType &in_out0, ModIntType &in_out1, const ModIntType &omega)
                {
                    auto x = in_out0.largeNorm();
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

                template <size_t MAX_LEN, uint64_t ROOT, typename ModIntType>
                struct NTTShort
                {
                    static constexpr size_t NTT_LEN = MAX_LEN;
                    static constexpr int LOG_LEN = hint_log2(NTT_LEN);
                    struct TableType
                    {
                        std::array<ModIntType, NTT_LEN> omega_table;
                        // Compute in compile time if need.
                        /*constexpr*/ TableType()
                        {
                            for (int omega_log_len = 0; omega_log_len <= LOG_LEN; omega_log_len++)
                            {
                                size_t omega_len = size_t(1) << omega_log_len, omega_count = omega_len / 2;
                                auto it = &omega_table[omega_len / 2];
                                ModIntType root = qpow(ModIntType(ROOT), (ModIntType::mod() - 1) / omega_len);
                                ModIntType omega(1);
                                for (size_t i = 0; i < omega_count; i++)
                                {
                                    it[i] = omega;
                                    omega *= root;
                                }
                            }
                        }
                        constexpr ModIntType &operator[](size_t i)
                        {
                            return omega_table[i];
                        }
                        constexpr const ModIntType &operator[](size_t i) const
                        {
                            return omega_table[i];
                        }
                        constexpr const ModIntType *getOmegaIt(size_t len) const
                        {
                            return &omega_table[len / 2];
                        }
                    };

                    static TableType table;

                    static void dit(ModIntType in_out[], size_t len)
                    {
                        len = std::min(NTT_LEN, len);
                        size_t rank = len;
                        if (hint_log2(len) % 2 == 0)
                        {
                            NTTShort<4, ROOT, ModIntType>::dit(in_out, len);
                            for (size_t i = 4; i < len; i += 4)
                            {
                                NTTShort<4, ROOT, ModIntType>::dit(in_out + i);
                            }
                            rank = 16;
                        }
                        else
                        {
                            NTTShort<8, ROOT, ModIntType>::dit(in_out, len);
                            for (size_t i = 8; i < len; i += 8)
                            {
                                NTTShort<8, ROOT, ModIntType>::dit(in_out + i);
                            }
                            rank = 32;
                        }
                        for (; rank <= len; rank *= 4)
                        {
                            size_t gap = rank / 4;
                            auto omega_it = table.getOmegaIt(rank), last_omega_it = table.getOmegaIt(rank / 2);
                            auto it0 = in_out, it1 = in_out + gap, it2 = in_out + gap * 2, it3 = in_out + gap * 3;
                            for (size_t j = 0; j < len; j += rank)
                            {
                                for (size_t i = 0; i < gap; i++)
                                {
                                    auto temp0 = it0[j + i], temp1 = it1[j + i], temp2 = it2[j + i], temp3 = it3[j + i], omega = last_omega_it[i];
                                    dit_butterfly2(temp0, temp1, omega);
                                    dit_butterfly2(temp2, temp3, omega);
                                    dit_butterfly2(temp0, temp2, omega_it[i]);
                                    dit_butterfly2(temp1, temp3, omega_it[gap + i]);
                                    it0[j + i] = temp0, it1[j + i] = temp1, it2[j + i] = temp2, it3[j + i] = temp3;
                                }
                            }
                        }
                    }
                    static void dif(ModIntType in_out[], size_t len)
                    {
                        len = std::min(NTT_LEN, len);
                        size_t rank = len;
                        for (; rank >= 16; rank /= 4)
                        {
                            size_t gap = rank / 4;
                            auto omega_it = table.getOmegaIt(rank), last_omega_it = table.getOmegaIt(rank / 2);
                            auto it0 = in_out, it1 = in_out + gap, it2 = in_out + gap * 2, it3 = in_out + gap * 3;
                            for (size_t j = 0; j < len; j += rank)
                            {
                                for (size_t i = 0; i < gap; i++)
                                {
                                    auto temp0 = it0[j + i], temp1 = it1[j + i], temp2 = it2[j + i], temp3 = it3[j + i], omega = last_omega_it[i];
                                    dif_butterfly2(temp0, temp2, omega_it[i]);
                                    dif_butterfly2(temp1, temp3, omega_it[gap + i]);
                                    dif_butterfly2(temp0, temp1, omega);
                                    dif_butterfly2(temp2, temp3, omega);
                                    it0[j + i] = temp0, it1[j + i] = temp1, it2[j + i] = temp2, it3[j + i] = temp3;
                                }
                            }
                        }
                        if (hint_log2(rank) % 2 == 0)
                        {
                            NTTShort<4, ROOT, ModIntType>::dif(in_out, len);
                            for (size_t i = 4; i < len; i += 4)
                            {
                                NTTShort<4, ROOT, ModIntType>::dif(in_out + i);
                            }
                        }
                        else
                        {
                            NTTShort<8, ROOT, ModIntType>::dif(in_out, len);
                            for (size_t i = 8; i < len; i += 8)
                            {
                                NTTShort<8, ROOT, ModIntType>::dif(in_out + i);
                            }
                        }
                    }
                };
                template <size_t LEN, uint64_t ROOT, typename ModIntType>
                typename NTTShort<LEN, ROOT, ModIntType>::TableType NTTShort<LEN, ROOT, ModIntType>::table;
                template <size_t LEN, uint64_t ROOT, typename ModIntType>
                constexpr size_t NTTShort<LEN, ROOT, ModIntType>::NTT_LEN;
                template <size_t LEN, uint64_t ROOT, typename ModIntType>
                constexpr int NTTShort<LEN, ROOT, ModIntType>::LOG_LEN;

                template <uint64_t ROOT, typename ModIntType>
                struct NTTShort<0, ROOT, ModIntType>
                {
                    static void dit(ModIntType in_out[]) {}
                    static void dif(ModIntType in_out[]) {}
                    static void dit(ModIntType in_out[], size_t len) {}
                    static void dif(ModIntType in_out[], size_t len) {}
                };

                template <uint64_t ROOT, typename ModIntType>
                struct NTTShort<1, ROOT, ModIntType>
                {
                    static void dit(ModIntType in_out[]) {}
                    static void dif(ModIntType in_out[]) {}
                    static void dit(ModIntType in_out[], size_t len) {}
                    static void dif(ModIntType in_out[], size_t len) {}
                };

                template <uint64_t ROOT, typename ModIntType>
                struct NTTShort<2, ROOT, ModIntType>
                {
                    static void dit(ModIntType in_out[])
                    {
                        transform2(in_out[0], in_out[1]);
                    }
                    static void dif(ModIntType in_out[])
                    {
                        transform2(in_out[0], in_out[1]);
                    }
                    static void dit(ModIntType in_out[], size_t len)
                    {
                        if (len < 2)
                        {
                            return;
                        }
                        dit(in_out);
                    }
                    static void dif(ModIntType in_out[], size_t len)
                    {
                        if (len < 2)
                        {
                            return;
                        }
                        dif(in_out);
                    }
                };

                template <uint64_t ROOT, typename ModIntType>
                struct NTTShort<4, ROOT, ModIntType>
                {
                    static void dit(ModIntType in_out[])
                    {
                        auto temp0 = in_out[0].largeNorm();
                        auto temp1 = in_out[1].largeNorm();
                        auto temp2 = in_out[2].largeNorm();
                        auto temp3 = in_out[3].largeNorm();

                        transform2(temp0, temp1);
                        auto sum = temp2.add(temp3);
                        auto dif = temp2.sub(temp3);
                        temp2 = sum.largeNorm();
                        temp3 = mul_w41<ROOT>(dif);

                        in_out[0] = temp0.add(temp2);
                        in_out[1] = temp1.add(temp3);
                        in_out[2] = temp0.sub(temp2);
                        in_out[3] = temp1.sub(temp3);
                    }
                    static void dif(ModIntType in_out[])
                    {
                        auto temp0 = in_out[0];
                        auto temp1 = in_out[1];
                        auto temp2 = in_out[2];
                        auto temp3 = in_out[3];

                        transform2(temp0, temp2);
                        auto sum = temp1.add(temp3);
                        auto dif = temp1.sub(temp3);
                        temp1 = sum.largeNorm();
                        temp3 = mul_w41<ROOT>(dif);

                        in_out[0] = temp0 + temp1;
                        in_out[1] = temp0 - temp1;
                        in_out[2] = temp2 + temp3;
                        in_out[3] = temp2 - temp3;
                    }
                    static void dit(ModIntType in_out[], size_t len)
                    {
                        if (len < 4)
                        {
                            NTTShort<2, ROOT, ModIntType>::dit(in_out, len);
                            return;
                        }
                        dit(in_out);
                    }
                    static void dif(ModIntType in_out[], size_t len)
                    {
                        if (len < 4)
                        {
                            NTTShort<2, ROOT, ModIntType>::dif(in_out, len);
                            return;
                        }
                        dif(in_out);
                    }
                };

                template <uint64_t ROOT, typename ModIntType>
                struct NTTShort<8, ROOT, ModIntType>
                {
                    static void dit(ModIntType in_out[])
                    {
                        static constexpr ModIntType w1 = qpow(ModIntType(ROOT), (ModIntType::mod() - 1) / 8);
                        static constexpr ModIntType w2 = qpow(w1, 2);
                        static constexpr ModIntType w3 = qpow(w1, 3);
                        auto temp0 = in_out[0].largeNorm();
                        auto temp1 = in_out[1].largeNorm();
                        auto temp2 = in_out[2].largeNorm();
                        auto temp3 = in_out[3].largeNorm();
                        auto temp4 = in_out[4].largeNorm();
                        auto temp5 = in_out[5].largeNorm();
                        auto temp6 = in_out[6].largeNorm();
                        auto temp7 = in_out[7].largeNorm();

                        transform2(temp0, temp1);
                        transform2(temp4, temp5);
                        auto sum = temp2.add(temp3);
                        auto dif = temp2.sub(temp3);
                        temp2 = sum.largeNorm();
                        temp3 = mul_w41<ROOT>(dif);
                        sum = temp6.add(temp7);
                        dif = temp6.sub(temp7);
                        temp6 = sum.largeNorm();
                        temp7 = mul_w41<ROOT>(dif);

                        transform2(temp0, temp2);
                        transform2(temp1, temp3);
                        sum = temp4.add(temp6);
                        dif = temp4.sub(temp6);
                        temp4 = sum.largeNorm();
                        temp6 = dif * w2;
                        sum = temp5.add(temp7);
                        dif = temp5.sub(temp7);
                        temp5 = sum * w1;
                        temp7 = dif * w3;

                        in_out[0] = temp0.add(temp4);
                        in_out[1] = temp1.add(temp5);
                        in_out[2] = temp2.add(temp6);
                        in_out[3] = temp3.add(temp7);
                        in_out[4] = temp0.sub(temp4);
                        in_out[5] = temp1.sub(temp5);
                        in_out[6] = temp2.sub(temp6);
                        in_out[7] = temp3.sub(temp7);
                    }
                    static void dif(ModIntType in_out[])
                    {
                        static constexpr ModIntType w1 = qpow(ModIntType(ROOT), (ModIntType::mod() - 1) / 8);
                        static constexpr ModIntType w2 = qpow(w1, 2);
                        static constexpr ModIntType w3 = qpow(w1, 3);
                        auto temp0 = in_out[0];
                        auto temp1 = in_out[1];
                        auto temp2 = in_out[2];
                        auto temp3 = in_out[3];
                        auto temp4 = in_out[4];
                        auto temp5 = in_out[5];
                        auto temp6 = in_out[6];
                        auto temp7 = in_out[7];

                        transform2(temp0, temp4);
                        auto sum = temp1.add(temp5);
                        auto dif = temp1.sub(temp5);
                        temp1 = sum.largeNorm();
                        temp5 = dif * w1;
                        sum = temp2.add(temp6);
                        dif = temp2.sub(temp6);
                        temp2 = sum.largeNorm();
                        temp6 = dif * w2;
                        sum = temp3.add(temp7);
                        dif = temp3.sub(temp7);
                        temp3 = sum.largeNorm();
                        temp7 = dif * w3;

                        transform2(temp0, temp2);
                        transform2(temp4, temp6);
                        sum = temp1.add(temp3);
                        dif = temp1.sub(temp3);
                        temp1 = sum.largeNorm();
                        temp3 = mul_w41<ROOT>(dif);
                        sum = temp5.add(temp7);
                        dif = temp5.sub(temp7);
                        temp5 = sum.largeNorm();
                        temp7 = mul_w41<ROOT>(dif);

                        in_out[0] = temp0 + temp1;
                        in_out[1] = temp0 - temp1;
                        in_out[2] = temp2 + temp3;
                        in_out[3] = temp2 - temp3;
                        in_out[4] = temp4 + temp5;
                        in_out[5] = temp4 - temp5;
                        in_out[6] = temp6 + temp7;
                        in_out[7] = temp6 - temp7;
                    }
                    static void dit(ModIntType in_out[], size_t len)
                    {
                        if (len < 8)
                        {
                            NTTShort<4, ROOT, ModIntType>::dit(in_out, len);
                            return;
                        }
                        dit(in_out);
                    }
                    static void dif(ModIntType in_out[], size_t len)
                    {
                        if (len < 8)
                        {
                            NTTShort<4, ROOT, ModIntType>::dif(in_out, len);
                            return;
                        }
                        dif(in_out);
                    }
                };

                template <uint64_t MOD, uint64_t ROOT, typename Int128Type = Uint128Default>
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
                    static_assert(check_inv<Int128Type>(root(), rootInv(), mod()), "IROOT * ROOT % MOD must be 1");
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

                    using INTT = NTT<mod(), rootInv(), Int128Type>;
                    using ModInt32Type = typename std::conditional<(MOD_BITS > 30), ModInt32<uint32_t(MOD)>, MontInt32Lazy<uint32_t(MOD)>>::type;
                    using ModInt64Type = MontInt64Lazy<MOD, Int128Type>;
                    using ModIntType = typename std::conditional<(MOD_BITS > 32), ModInt64Type, ModInt32Type>::type;
                    using IntType = typename ModIntType::IntType;

                    static constexpr size_t L2_BYTE = size_t(1) << 20; // 1MB L2 cache size, change this if you know your cache size.
                    static constexpr size_t LONG_THRESHOLD = std::min(L2_BYTE / sizeof(ModIntType), NTT_MAX_LEN);
                    using NTTTemplate = NTTShort<LONG_THRESHOLD, root(), ModIntType>;

                    static void dit244(ModIntType in_out[], size_t ntt_len)
                    {
                        ntt_len = std::min(int_floor2(ntt_len), NTT_MAX_LEN);
                        if (ntt_len <= LONG_THRESHOLD)
                        {
                            NTTTemplate::dit(in_out, ntt_len);
                            return;
                        }
                        size_t quarter_len = ntt_len / 4;
                        dit244(in_out + quarter_len * 3, ntt_len / 4);
                        dit244(in_out + quarter_len * 2, ntt_len / 4);
                        dit244(in_out, ntt_len / 2);
                        const ModIntType unit_omega1 = qpow(ModIntType(root()), (mod() - 1) / ntt_len);
                        const ModIntType unit_omega3 = qpow(unit_omega1, 3);
                        ModIntType omega1(1), omega3(1);
                        auto it0 = in_out, it1 = in_out + quarter_len, it2 = in_out + quarter_len * 2, it3 = in_out + quarter_len * 3;
                        for (size_t i = 0; i < quarter_len; i++)
                        {
                            ModIntType temp0 = it0[i], temp1 = it1[i], temp2 = it2[i] * omega1, temp3 = it3[i] * omega3;
                            dit_butterfly244<ROOT>(temp0, temp1, temp2, temp3);
                            it0[i] = temp0, it1[i] = temp1, it2[i] = temp2, it3[i] = temp3;
                            omega1 = omega1 * unit_omega1;
                            omega3 = omega3 * unit_omega3;
                        }
                    }
                    static void dif244(ModIntType in_out[], size_t ntt_len)
                    {
                        ntt_len = std::min(int_floor2(ntt_len), NTT_MAX_LEN);
                        if (ntt_len <= LONG_THRESHOLD)
                        {
                            NTTTemplate::dif(in_out, ntt_len);
                            return;
                        }
                        size_t quarter_len = ntt_len / 4;
                        const ModIntType unit_omega1 = qpow(ModIntType(root()), (mod() - 1) / ntt_len);
                        const ModIntType unit_omega3 = qpow(unit_omega1, 3);
                        ModIntType omega1(1), omega3(1);
                        auto it0 = in_out, it1 = in_out + quarter_len, it2 = in_out + quarter_len * 2, it3 = in_out + quarter_len * 3;
                        for (size_t i = 0; i < quarter_len; i++)
                        {
                            ModIntType temp0 = it0[i], temp1 = it1[i], temp2 = it2[i], temp3 = it3[i];
                            dif_butterfly244<ROOT>(temp0, temp1, temp2, temp3);
                            it0[i] = temp0, it1[i] = temp1, it2[i] = temp2 * omega1, it3[i] = temp3 * omega3;
                            omega1 = omega1 * unit_omega1;
                            omega3 = omega3 * unit_omega3;
                        }
                        dif244(in_out, ntt_len / 2);
                        dif244(in_out + quarter_len * 3, ntt_len / 4);
                        dif244(in_out + quarter_len * 2, ntt_len / 4);
                    }
                    static void convolution(ModIntType in1[], ModIntType in2[], ModIntType out[], size_t ntt_len, bool normlize = true)
                    {
                        dif244(in1, ntt_len);
                        dif244(in2, ntt_len);
                        if (normlize)
                        {
                            const ModIntType inv_len(qpow(ModIntType(ntt_len), mod() - 2));
                            for (size_t i = 0; i < ntt_len; i++)
                            {
                                out[i] = in1[i] * in2[i] * inv_len;
                            }
                        }
                        else
                        {
                            for (size_t i = 0; i < ntt_len; i++)
                            {
                                out[i] = in1[i] * in2[i];
                            }
                        }
                        INTT::dit244(out, ntt_len);
                    }
                    static void convolutionRecursion(ModIntType in1[], ModIntType in2[], ModIntType out[], size_t ntt_len, bool normlize = true)
                    {
                        if (ntt_len <= LONG_THRESHOLD)
                        {
                            NTTTemplate::dif(in1, ntt_len);
                            if (in1 != in2)
                            {
                                NTTTemplate::dif(in2, ntt_len);
                            }
                            if (normlize)
                            {
                                const ModIntType inv_len(qpow(ModIntType(ntt_len), mod() - 2));
                                for (size_t i = 0; i < ntt_len; i++)
                                {
                                    out[i] = in1[i] * in2[i] * inv_len;
                                }
                            }
                            else
                            {
                                for (size_t i = 0; i < ntt_len; i++)
                                {
                                    out[i] = in1[i] * in2[i];
                                }
                            }
                            INTT::NTTTemplate::dit(out, ntt_len);
                            return;
                        }
                        const size_t quarter_len = ntt_len / 4;
                        ModIntType unit_omega1 = qpow(ModIntType(root()), (mod() - 1) / ntt_len);
                        ModIntType unit_omega3 = qpow(unit_omega1, 3);
                        ModIntType omega1(1), omega3(1);
                        if (in1 != in2)
                        {
                            for (size_t i = 0; i < quarter_len; i++)
                            {
                                ModIntType temp0 = in1[i], temp1 = in1[quarter_len + i], temp2 = in1[quarter_len * 2 + i], temp3 = in1[quarter_len * 3 + i];
                                dif_butterfly244<ROOT>(temp0, temp1, temp2, temp3);
                                in1[i] = temp0, in1[quarter_len + i] = temp1, in1[quarter_len * 2 + i] = temp2 * omega1, in1[quarter_len * 3 + i] = temp3 * omega3;

                                temp0 = in2[i], temp1 = in2[quarter_len + i], temp2 = in2[quarter_len * 2 + i], temp3 = in2[quarter_len * 3 + i];
                                dif_butterfly244<ROOT>(temp0, temp1, temp2, temp3);
                                in2[i] = temp0, in2[quarter_len + i] = temp1, in2[quarter_len * 2 + i] = temp2 * omega1, in2[quarter_len * 3 + i] = temp3 * omega3;

                                omega1 = omega1 * unit_omega1;
                                omega3 = omega3 * unit_omega3;
                            }
                        }
                        else
                        {
                            for (size_t i = 0; i < quarter_len; i++)
                            {
                                ModIntType temp0 = in1[i], temp1 = in1[quarter_len + i], temp2 = in1[quarter_len * 2 + i], temp3 = in1[quarter_len * 3 + i];
                                dif_butterfly244<ROOT>(temp0, temp1, temp2, temp3);
                                in1[i] = temp0, in1[quarter_len + i] = temp1, in1[quarter_len * 2 + i] = temp2 * omega1, in1[quarter_len * 3 + i] = temp3 * omega3;

                                omega1 = omega1 * unit_omega1;
                                omega3 = omega3 * unit_omega3;
                            }
                        }

                        convolutionRecursion(in1, in2, out, ntt_len / 2, false);
                        convolutionRecursion(in1 + quarter_len * 2, in2 + quarter_len * 2, out + quarter_len * 2, ntt_len / 4, false);
                        convolutionRecursion(in1 + quarter_len * 3, in2 + quarter_len * 3, out + quarter_len * 3, ntt_len / 4, false);

                        unit_omega1 = qpow(ModIntType(rootInv()), (mod() - 1) / ntt_len);
                        unit_omega3 = qpow(unit_omega1, 3);
                        if (normlize)
                        {
                            const ModIntType inv_len(qpow(ModIntType(ntt_len), mod() - 2));
                            omega1 = inv_len, omega3 = inv_len;
                            for (size_t i = 0; i < quarter_len; i++)
                            {
                                ModIntType temp0 = out[i] * inv_len, temp1 = out[quarter_len + i] * inv_len, temp2 = out[quarter_len * 2 + i] * omega1, temp3 = out[quarter_len * 3 + i] * omega3;
                                dit_butterfly244<rootInv()>(temp0, temp1, temp2, temp3);
                                out[i] = temp0, out[quarter_len + i] = temp1, out[quarter_len * 2 + i] = temp2, out[quarter_len * 3 + i] = temp3;

                                omega1 = omega1 * unit_omega1;
                                omega3 = omega3 * unit_omega3;
                            }
                        }
                        else
                        {
                            omega1 = 1, omega3 = 1;
                            for (size_t i = 0; i < quarter_len; i++)
                            {
                                ModIntType temp0 = out[i], temp1 = out[quarter_len + i], temp2 = out[quarter_len * 2 + i] * omega1, temp3 = out[quarter_len * 3 + i] * omega3;
                                dit_butterfly244<rootInv()>(temp0, temp1, temp2, temp3);
                                out[i] = temp0, out[quarter_len + i] = temp1, out[quarter_len * 2 + i] = temp2, out[quarter_len * 3 + i] = temp3;

                                omega1 = omega1 * unit_omega1;
                                omega3 = omega3 * unit_omega3;
                            }
                        }
                    }
                };
                template <uint64_t MOD, uint64_t ROOT, typename Int128Type>
                constexpr int NTT<MOD, ROOT, Int128Type>::MOD_BITS;
                template <uint64_t MOD, uint64_t ROOT, typename Int128Type>
                constexpr int NTT<MOD, ROOT, Int128Type>::MAX_LOG_LEN;
                template <uint64_t MOD, uint64_t ROOT, typename Int128Type>
                constexpr size_t NTT<MOD, ROOT, Int128Type>::NTT_MAX_LEN;
            } // namespace split_radix
            using NTT0 = split_radix::NTT<MOD0, ROOT0>; // using 64bit integer, Montgomery speed up
            using NTT1 = split_radix::NTT<MOD1, ROOT1>; // using 64bit integer, Montgomery speed up
            using NTT2 = split_radix::NTT<MOD2, ROOT2>; // using 64bit integer, Montgomery speed up
            using NTT3 = split_radix::NTT<MOD3, ROOT3>; // using 32bit integer, Montgomery speed up
            using NTT4 = split_radix::NTT<MOD4, ROOT4>; // using 32bit integer, Montgomery speed up
            using NTT5 = split_radix::NTT<MOD5, ROOT5>; // using 32bit integer
            using NTT6 = split_radix::NTT<MOD6, ROOT6>; // using 32bit integer
        }
    } // namespace transform
} // namespace hint
#endif