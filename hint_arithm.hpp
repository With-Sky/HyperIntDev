#include "hint_basic.hpp"

#ifndef HINT_ARITHM_HPP
#define HINT_ARITHM_HPP
namespace hint
{
    namespace arithm
    {
        using namespace hint::extend_int;
        namespace support
        {
            template <typename T>
            void ary_print(const T a[], size_t len, bool is_rev = false)
            {
                if (len == 0)
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
                std::cout << "]\n";
            }

            template <typename T>
            void store_arr(const T arr[], size_t len, const std::string &path)
            {
                if (len == 0)
                {
                    return;
                }
                std::ofstream file(path);
                for (size_t i = 0; i < len - 1; i++)
                {
                    file << arr[i] << "\n";
                }
                file << arr[len - 1];
            }

            template <typename T>
            void load_arr(std::vector<T> &arr, const std::string &path)
            {
                std::ifstream file(path);
                arr.clear();
                T num;
                while (file >> num)
                {
                    arr.push_back(num);
                }
            }

            // 去除前导零，返回实际长度，如果数组为空，则返回0
            template <typename T>
            constexpr size_t count_ture_length(const T array[], size_t length)
            {
                if (array == nullptr)
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
            constexpr void remove_leading_zeros(const T array[], size_t &length)
            {
                length = count_ture_length(array, length);
            }

            // Absolute compare, return 1 if a > b, -1 if a < b, 0 if a == b
            // Return the diffence length if a != b
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

            // 二进制移位,shift小于字的位数,输入输出可以为相同指针
            template <typename WordTy>
            constexpr WordTy lshift_in_word_half(const WordTy in[], size_t len, WordTy out[], int shift)
            {
                constexpr int WORD_BITS = sizeof(WordTy) * CHAR_BIT;
                assert(shift >= 0 && shift < WORD_BITS);
                if (0 == len)
                {
                    return 0;
                }
                if (0 == shift)
                {
                    std::copy(in, in + len, out);
                    return 0;
                }
                // [n,last] -> [?,n >> shift_rem | last << shift]
                WordTy last = in[len - 1], ret = last;
                const int shift_rem = WORD_BITS - shift;
                size_t i = len - 1;
                while (i > 0)
                {
                    i--;
                    WordTy n = in[i];
                    out[i + 1] = (last << shift) | (n >> shift_rem);
                    last = n;
                }
                out[0] = last << shift;
                return ret >> shift_rem;
            }
            template <typename WordTy>
            constexpr void lshift_in_word(const WordTy in[], size_t len, WordTy out[], int shift)
            {
                if (0 == len)
                {
                    return;
                }
                assert(shift >= 0 && size_t(shift) < sizeof(WordTy) * CHAR_BIT);
                uint64_t last = lshift_in_word_half(in, len, out, shift);
                out[len] = last;
            }

            template <typename WordTy>
            constexpr void lshift_binary(const WordTy in[], size_t len, WordTy out[], size_t shift)
            {
                if (0 == len)
                {
                    return;
                }
                constexpr int WORD_BIT = sizeof(WordTy) * CHAR_BIT;
                size_t block_shift = shift / WORD_BIT, in_word_shift = shift % WORD_BIT;
                if (0 == in_word_shift)
                {
                    std::copy(in, in + len, out + block_shift);
                    return;
                }
                lshift_in_word(in, len, out + block_shift, in_word_shift);
                std::fill_n(out, block_shift, WordTy{});
            }

            template <typename WordTy>
            constexpr void rshift_in_word(const WordTy in[], size_t len, WordTy out[], int shift)
            {
                constexpr int WORD_BITS = sizeof(WordTy) * CHAR_BIT;
                if (0 == len)
                {
                    return;
                }
                if (0 == shift)
                {
                    std::copy(in, in + len, out);
                    return;
                }
                assert(shift >= 0 && size_t(shift) < sizeof(WordTy) * CHAR_BIT);
                WordTy last = in[0];
                const int shift_rem = WORD_BITS - shift;
                for (size_t i = 1; i < len; i++)
                {
                    WordTy n = in[i];
                    out[i - 1] = (last >> shift) | (n << shift_rem);
                    last = n;
                }
                out[len - 1] = last >> shift;
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
                if (l_len == 0 || r_len == 0)
                {
                    return 0;
                }
                return l_len + r_len;
            }

            size_t get_div_len(size_t l_len, size_t r_len)
            {
                if (l_len < r_len)
                {
                    return 0;
                }
                return l_len - r_len + 1;
            }

            template <typename NumTy, typename ProdTy>
            class DivSupporter
            {
            private:
                NumTy divisor = 0;
                NumTy inv = 0;
                int shift = 0;
                enum : int
                {
                    NUM_BITS = sizeof(NumTy) * CHAR_BIT
                };

            public:
                constexpr DivSupporter(NumTy divisor_in) : divisor(divisor_in)
                {
                    inv = getInv(divisor, shift);
                    divisor <<= shift;
                }
                // Return dividend / divisor, dividend %= divisor
                NumTy divMod(ProdTy &dividend) const
                {
                    dividend <<= shift;
                    NumTy r = NumTy(dividend);
                    dividend = (dividend >> NUM_BITS) * inv + dividend;
                    NumTy q1 = NumTy(dividend >> NUM_BITS) + 1;
                    r -= q1 * divisor;
                    if (r > NumTy(dividend))
                    {
                        q1--;
                        r += divisor;
                    }
                    if (r >= divisor)
                    {
                        q1++;
                        r -= divisor;
                    }
                    dividend = r >> shift;
                    return q1;
                }
                NumTy div(ProdTy dividend) const
                {
                    return divMod(dividend);
                }
                NumTy mod(ProdTy dividend) const
                {
                    divMod(dividend);
                    return dividend;
                }

                static constexpr NumTy getInv(NumTy divisor, int &leading_zero)
                {
                    constexpr NumTy MAX = hint::all_one<NumTy>(NUM_BITS);
                    leading_zero = hint::hint_clz(divisor);
                    divisor <<= leading_zero;
                    ProdTy x = ProdTy(MAX - divisor) << NUM_BITS;
                    return NumTy((x + MAX) / divisor);
                }
            };

            template <typename NumTy>
            class BaseExecutor
            {
            private:
                static constexpr int NUM_BIT = sizeof(NumTy) * CHAR_BIT;
                using ProdTy = typename UintType<NUM_BIT>::NextType::Type;
                using SignTy = typename UintType<NUM_BIT>::SignType;
                using DivSupporterTy = DivSupporter<NumTy, ProdTy>;
                NumTy base;
                DivSupporterTy div_supporter;

            public:
                constexpr BaseExecutor(NumTy base_in) : base(base_in), div_supporter(base_in) {}

                constexpr NumTy addCarry(NumTy a, NumTy b, bool &cf) const
                {
                    a = a + b + cf;
                    cf = (a >= base);
                    if (cf)
                    {
                        a -= base;
                    }
                    return a;
                }
                constexpr NumTy addHalf(NumTy a, NumTy b, bool &cf) const
                {
                    a = a + b;
                    cf = (a >= base);
                    if (cf)
                    {
                        a -= base;
                    }
                    return a;
                }

                constexpr NumTy subBorrow(SignTy a, SignTy b, bool &bf) const
                {
                    a = a - b - bf;
                    bf = (a < 0);
                    if (bf)
                    {
                        a += base;
                    }
                    return a;
                }
                constexpr NumTy subHalf(SignTy a, SignTy b, bool &bf) const
                {
                    a = a - b;
                    bf = (a < 0);
                    if (bf)
                    {
                        a += base;
                    }
                    return a;
                }

                NumTy divMod(ProdTy &dividend) const
                {
                    return div_supporter.divMod(dividend);
                }
            };
            template <typename NumTy>
            constexpr int BaseExecutor<NumTy>::NUM_BIT;
        }

        namespace addition_binary
        {
            using namespace support;

            template <typename NumTy, typename Executor>
            constexpr bool abs_add_long(const NumTy a[], const NumTy b[], size_t len, NumTy sum[], const Executor &exec)
            {
                bool carry = false;
                size_t i = 0;
                for (const size_t rem_len = len - len % 4; i < rem_len; i += 4)
                {
                    sum[i] = exec.addCarry(a[i], b[i], carry);
                    sum[i + 1] = exec.addCarry(a[i + 1], b[i + 1], carry);
                    sum[i + 2] = exec.addCarry(a[i + 2], b[i + 2], carry);
                    sum[i + 3] = exec.addCarry(a[i + 3], b[i + 3], carry);
                }
                for (; i < len; i++)
                {
                    sum[i] = exec.addCarry(a[i], b[i], carry);
                }
                return carry;
            }

            template <typename NumTy, typename Executor>
            constexpr bool abs_add_long_num(const NumTy a[], size_t len, NumTy num, NumTy sum[], const Executor &exec)
            {
                assert(len > 0);
                bool carry = false;
                sum[0] = exec.addHalf(a[0], num, carry);
                for (size_t i = 1; i < len; i++)
                {
                    sum[i] = exec.addHalf(a[i], NumTy(carry), carry);
                }
                return carry;
            }

            template <typename NumTy, typename Executor>
            constexpr bool abs_add_long_long(const NumTy a[], size_t len_a, const NumTy b[], size_t len_b, NumTy sum[], const Executor &exec)
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
                return carry;
            }

            template <typename NumTy, typename Executor>
            constexpr bool abs_sub_long(const NumTy a[], const NumTy b[], size_t len, NumTy diff[], const Executor &exec)
            {
                bool borrow = false;
                for (size_t i = 0; i < len; i++)
                {
                    diff[i] = exec.subBorrow(a[i], b[i], borrow);
                }
                return borrow;
            }

            template <typename NumTy, typename Executor>
            constexpr bool abs_sub_long_num(const NumTy a[], size_t len, NumTy num, NumTy diff[], const Executor &exec)
            {
                assert(len > 0);
                bool borrow = false;
                diff[0] = exec.subHalf(a[0], num, borrow);
                for (size_t i = 1; i < len; i++)
                {
                    diff[i] = exec.subHalf(a[i], NumTy(borrow), borrow);
                }
                return borrow;
            }

            template <typename NumTy, typename Executor>
            constexpr bool abs_sub_long_long(const NumTy a[], size_t len_a, const NumTy b[], size_t len_b, NumTy diff[], const Executor &exec)
            {
                assert(len_a >= len_b);
                bool borrow = abs_sub_long(a, b, len_b, diff, exec);
                if (len_a > len_b)
                {
                    borrow = abs_sub_long_num(a + len_b, len_a - len_b, NumTy(borrow), diff + len_b, exec);
                }
                return borrow;
            }

            // Binary absolute addtion a+b=sum, return the carry
            template <typename UintTy>
            constexpr bool abs_add_binary_half(const UintTy a[], size_t len_a, const UintTy b[], size_t len_b, UintTy sum[])
            {
                bool carry = false;
                size_t i = 0, min_len = std::min(len_a, len_b);
                for (; i < min_len; i++)
                {
                    sum[i] = add_carry(a[i], b[i], carry);
                }
                for (; i < len_a; i++)
                {
                    sum[i] = add_half(a[i], UintTy(carry), carry);
                }
                for (; i < len_b; i++)
                {
                    sum[i] = add_half(b[i], UintTy(carry), carry);
                }
                return carry;
            }
            // Binary absolute addtion a+b=sum, return the carry
            template <typename UintTy>
            constexpr void abs_add_binary(const UintTy a[], size_t len_a, const UintTy b[], size_t len_b, UintTy sum[])
            {
                bool carry = abs_add_binary_half(a, len_a, b, len_b, sum);
                sum[std::max(len_a, len_b)] = carry;
            }

            // Binary absolute subtraction a-b=diff, return the borrow
            template <typename UintTy>
            constexpr bool abs_sub_binary(const UintTy a[], size_t len_a, const UintTy b[], size_t len_b, UintTy diff[])
            {
                bool borrow = false;
                size_t i = 0, min_len = std::min(len_a, len_b), rem_len = min_len - min_len % 4;
                for (; i < rem_len; i += 4)
                {
                    diff[i] = sub_borrow(a[i], b[i], borrow);
                    diff[i + 1] = sub_borrow(a[i + 1], b[i + 1], borrow);
                    diff[i + 2] = sub_borrow(a[i + 2], b[i + 2], borrow);
                    diff[i + 3] = sub_borrow(a[i + 3], b[i + 3], borrow);
                }
                // borrow = bf != 0;
                for (; i < min_len; i++)
                {
                    diff[i] = sub_borrow(a[i], b[i], borrow);
                }
                for (; i < len_a; i++)
                {
                    diff[i] = sub_half(a[i], UintTy(borrow), borrow);
                }
                for (; i < len_b; i++)
                {
                    UintTy b0 = b[i];
                    diff[i] = UintTy(0) - borrow - b0;
                    borrow = borrow || (b0 > 0);
                }
                return borrow;
            }

            template <typename UintTy>
            [[nodiscard]] constexpr int abs_difference_binary(const UintTy a[], size_t len1, const UintTy b[], size_t len2, UintTy diff[])
            {
                auto cmp_ab = abs_compare_with_length(a, len1, b, len2);
                if (len1 == len2)
                {
                    std::fill(diff + cmp_ab.diff_len, diff + len1, UintTy(0));
                    len1 = len2 = cmp_ab.diff_len;
                }
                if (cmp_ab.cmp < 0)
                {
                    std::swap(a, b);
                    std::swap(len1, len2);
                }
                abs_sub_binary(a, len1, b, len2, diff);
                return cmp_ab.cmp;
            }
            template <typename UintTy>
            [[nodiscard]] bool addWithSign(const UintTy in1[], size_t len1, bool sign1, const UintTy in2[], size_t len2, bool sign2, UintTy out[])
            {
                bool res_sign = false;
                if (sign1 == sign2)
                {
                    bool carry = arithm::addition_binary::abs_add_binary(in1, len1, in2, len2, out);
                    return sign1;
                }
                int cmp = arithm::addition_binary::abs_difference_binary(in1, len1, in2, len2, out);
                if (0 == cmp)
                {
                    return false;
                }
                return cmp > 0 ? sign1 : sign2;
            }

            // a - num
            template <typename UintTy>
            constexpr bool abs_sub_num_binary(const UintTy a[], size_t len_a, UintTy num, UintTy diff[])
            {
                assert(len_a > 0);
                bool borrow = false;
                size_t i = 1;
                diff[0] = sub_half(a[0], num, borrow);
                for (; i < len_a; i++)
                {
                    diff[i] = sub_half(a[i], UintTy(borrow), borrow);
                }
                return borrow;
            }

            // a + num
            template <typename UintTy>
            constexpr bool abs_add_num_binary_half(const UintTy a[], size_t len_a, UintTy num, UintTy sum[])
            {
                assert(len_a > 0);
                bool carry = false;
                size_t i = 1;
                sum[0] = add_half(a[0], num, carry);
                for (; i < len_a; i++)
                {
                    sum[i] = add_half(a[i], UintTy(carry), carry);
                }
                return carry;
            }

            // a + num
            template <typename UintTy>
            constexpr void abs_add_num_binary(const UintTy a[], size_t len_a, UintTy num, UintTy sum[])
            {
                assert(len_a > 0);
                sum[len_a] = abs_add_num_binary_half(a, len_a, num, sum);
            }

            // Binary absolute addtion a+b=sum, return the carry
            template <typename UintTy>
            constexpr bool binary_complement_code(const UintTy a[], size_t len_a, UintTy complement[])
            {
                assert(len_a > 0);
                bool borrow = false;
                for (size_t i = 0; i < len_a; i++)
                {
                    UintTy n = a[i];
                    complement[i] = UintTy(0) - n - borrow;
                    borrow = n > 0 || borrow;
                }
                return borrow;
            }

        }

        namespace multiplication
        {
            using namespace addition_binary;

            inline uint64_t abs_mul_add_num64_half(const uint64_t in[], size_t len, uint64_t out[], uint64_t num_add, uint64_t num_mul)
            {
                size_t i = 0;
                uint64_t prod_lo, prod_hi;
                for (const size_t rem_len = len - len % 4; i < rem_len; i += 4)
                {
                    mul64x64to128(in[i], num_mul, prod_lo, prod_hi);
                    prod_lo += num_add;
                    out[i] = prod_lo;
                    num_add = prod_hi + (prod_lo < num_add);

                    mul64x64to128(in[i + 1], num_mul, prod_lo, prod_hi);
                    prod_lo += num_add;
                    out[i + 1] = prod_lo;
                    num_add = prod_hi + (prod_lo < num_add);

                    mul64x64to128(in[i + 2], num_mul, prod_lo, prod_hi);
                    prod_lo += num_add;
                    out[i + 2] = prod_lo;
                    num_add = prod_hi + (prod_lo < num_add);

                    mul64x64to128(in[i + 3], num_mul, prod_lo, prod_hi);
                    prod_lo += num_add;
                    out[i + 3] = prod_lo;
                    num_add = prod_hi + (prod_lo < num_add);
                }
                for (; i < len; i++)
                {
                    mul64x64to128(in[i], num_mul, prod_lo, prod_hi);
                    prod_lo += num_add;
                    out[i] = prod_lo;
                    num_add = prod_hi + (prod_lo < num_add);
                }
                return num_add;
            }

            /// @brief 2^64 base long integer multiply 64bit number, add another 64bit number to product.
            /// @param in Input long integer.
            /// @param len Length of input long integer.
            /// @param out Output long integer, equals to input * num_mul + num_add
            /// @param num_add The 64 bit number to add.
            /// @param num_mul The 64 bit number to multiply.
            inline void abs_mul_add_num64(const uint64_t in[], size_t len, uint64_t out[], uint64_t num_add, uint64_t num_mul)
            {
                out[len] = abs_mul_add_num64_half(in, len, out, num_add, num_mul);
            }

            // in * num_mul + in_out -> in_out
            inline void mul64_sub_proc(const uint64_t in[], size_t len, uint64_t in_out[], uint64_t num_mul)
            {
                uint64_t carry = 0;
                size_t i = 0;
                for (const size_t rem_len = len - len % 2; i < rem_len; i += 2)
                {
                    bool cf;
                    uint64_t prod_lo, prod_hi;
                    mul64x64to128(in[i], num_mul, prod_lo, prod_hi);
                    prod_lo = add_half(prod_lo, in_out[i], cf);
                    prod_hi += cf;
                    in_out[i] = add_half(prod_lo, carry, cf);
                    carry = prod_hi + cf;

                    mul64x64to128(in[i + 1], num_mul, prod_lo, prod_hi);
                    prod_lo = add_half(prod_lo, in_out[i + 1], cf);
                    prod_hi += cf;
                    in_out[i + 1] = add_half(prod_lo, carry, cf);
                    carry = prod_hi + cf;

                    mul64x64to128(in[i + 2], num_mul, prod_lo, prod_hi);
                    prod_lo = add_half(prod_lo, in_out[i + 2], cf);
                    prod_hi += cf;
                    in_out[i + 2] = add_half(prod_lo, carry, cf);
                    carry = prod_hi + cf;

                    mul64x64to128(in[i + 3], num_mul, prod_lo, prod_hi);
                    prod_lo = add_half(prod_lo, in_out[i + 3], cf);
                    prod_hi += cf;
                    in_out[i + 3] = add_half(prod_lo, carry, cf);
                    carry = prod_hi + cf;
                }
                for (; i < len; i++)
                {
                    bool cf;
                    uint64_t prod_lo, prod_hi;
                    mul64x64to128(in[i], num_mul, prod_lo, prod_hi);
                    prod_lo = add_half(prod_lo, in_out[i], cf);
                    prod_hi += cf;
                    in_out[i] = add_half(prod_lo, carry, cf);
                    carry = prod_hi + cf;
                }
                in_out[len] = carry;
            }

            // 小学乘法
            inline void abs_mul64_classic(const uint64_t in1[], size_t len1, const uint64_t in2[], size_t len2, uint64_t out[], uint64_t *work_begin = nullptr, uint64_t *work_end = nullptr)
            {
                const size_t out_len = get_mul_len(len1, len2);
                remove_leading_zeros(in1, len1);
                remove_leading_zeros(in2, len2);
                if (len1 < len2)
                {
                    std::swap(in1, in2);
                    std::swap(len1, len2); // Let in1 be the loonger one
                }
                if (0 == len2 || nullptr == in1 || nullptr == in2)
                {
                    std::fill_n(out, out_len, uint64_t(0));
                    return;
                }
                if (1 == len2)
                {
                    abs_mul_add_num64(in1, len1, out, 0, in2[0]);
                    return;
                }
                // Get enough work memory
                std::vector<uint64_t> work_mem;
                const size_t work_size = get_mul_len(len1, len2);
                if (work_begin + work_size > work_end)
                {
                    work_mem.resize(work_size);
                    work_begin = work_mem.data();
                    work_end = work_begin + work_mem.size();
                }
                else
                {
                    // Clear work_mem that may used
                    std::fill_n(work_begin, work_size, uint64_t(0));
                }
                auto out_temp = work_begin;
                for (size_t i = 0; i < len1; i++)
                {
                    mul64_sub_proc(in2, len2, out_temp + i, in1[i]);
                }
                std::copy(out_temp, out_temp + work_size, out);
                std::fill(out + work_size, out + out_len, uint64_t(0));
            }

            inline void mul_check(const uint64_t in1[], size_t len1, const uint64_t in2[], size_t len2, const uint64_t out[])
            {
                std::vector<uint64_t> out2(len1 + len2);
                abs_mul64_classic(in1, len1, in2, len2, out2.data());
                for (size_t i = 0; i < len1 + len2; i++)
                {
                    if (out[i] != out2[i])
                    {
                        std::cout << "mul check fail" << std::endl;
                        store_arr(in1, len1, "in1.txt");
                        store_arr(in2, len2, "in2.txt");
                        exit(1);
                    }
                }
                std::cout << "mul check success" << std::endl;
            }

            // Karatsuba 乘法
            inline void abs_mul64_karatsuba(const uint64_t in1[], size_t len1, const uint64_t in2[], size_t len2, uint64_t out[], uint64_t *work_begin = nullptr, uint64_t *work_end = nullptr)
            {
                const size_t out_len = get_mul_len(len1, len2);
                remove_leading_zeros(in1, len1);
                remove_leading_zeros(in2, len2);
                if (len1 < len2)
                {
                    std::swap(in1, in2);
                    std::swap(len1, len2); // Let in1 be the loonger one
                }
                if (0 == len2 || nullptr == in1 || nullptr == in2)
                {
                    std::fill_n(out, out_len, uint64_t(0));
                    return;
                }
                constexpr size_t KARATSUBA_THRESHOLD = 24;
                if (len2 < KARATSUBA_THRESHOLD)
                {
                    abs_mul64_classic(in1, len1, in2, len2, out, work_begin, work_end);
                    std::fill(out + get_mul_len(len1, len2), out + out_len, uint64_t(0));
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
                // Get length of every part
                size_t m_len = get_mul_len(len1_low, len2_low);
                size_t n_len = get_mul_len(len1_high, len2_high);

                // Get enough work_mem
                std::vector<uint64_t> work_mem;
                const size_t work_size = m_len + n_len + get_mul_len(len1_low, len2_low);
                if (work_begin + work_size > work_end)
                {
                    work_mem.resize(work_size * 2);
                    work_begin = work_mem.data();
                    work_end = work_begin + work_mem.size();
                }
                // Set pointer of every part
                auto m = work_begin, n = m + m_len, k1 = n + n_len, k2 = k1 + len1_low, k = k1;

                // Compute M,N
                abs_mul64_karatsuba(in1, len1_low, in2, len2_low, m, work_begin + work_size, work_end); // M = AL * BL
                if (len2_high > 0)
                {
                    abs_mul64_karatsuba(in1 + base_len, len1_high, in2 + base_len, len2_high, n, work_begin + work_size, work_end); // N = AH * BH
                }

                // Compute K1,K2
                remove_leading_zeros(in1, len1_low);
                remove_leading_zeros(in2, len2_low);
                int cmp1 = abs_difference_binary(in1, len1_low, in1 + base_len, len1_high, k1); // k1 = abs(AH - AL)
                int cmp2 = abs_difference_binary(in2, len2_low, in2 + base_len, len2_high, k2); // k2 = abs(BH - BL)
                size_t k1_len = get_sub_len(len1_low, len1_high);
                size_t k2_len = get_sub_len(len2_low, len2_high);
                remove_leading_zeros(k1, k1_len);
                remove_leading_zeros(k2, k2_len);

                // Compute K1*K2 = K
                abs_mul64_karatsuba(k1, k1_len, k2, k2_len, k, work_begin + work_size, work_end);
                size_t k_len = get_mul_len(k1_len, k2_len);
                remove_leading_zeros(k, k_len);

                // Combine the result
                // out = M + N * BASE ^ 2
                {
                    std::copy(m, m + m_len, out);
                    std::fill(out + m_len, out + base_len * 2, uint64_t(0));
                    std::copy(n, n + n_len, out + base_len * 2);
                    std::fill(out + base_len * 2 + n_len, out + out_len, uint64_t(0));
                }
                auto out_base = out + base_len;
                // out = M + N * BASE ^ 2 + (M + N) ^ BASE
                {
                    m_len = std::min(m_len, out_len - base_len);
                    n_len = std::min(n_len, out_len - base_len);
                    if (m_len < n_len)
                    {
                        std::swap(m_len, n_len);
                        std::swap(m, n);
                    }
                    int carry = 0;
                    bool cf;
                    size_t i = 0;
                    for (; i < n_len; i++)
                    {
                        uint64_t sum = add_half(m[i], uint64_t(carry), cf);
                        carry = cf;
                        sum = add_half(n[i], sum, cf);
                        carry += cf;
                        out_base[i] = add_half(out_base[i], sum, cf);
                        carry += cf;
                    }
                    for (; i < m_len; i++)
                    {
                        uint64_t sum = add_half(m[i], uint64_t(carry), cf);
                        carry = cf;
                        out_base[i] = add_half(out_base[i], sum, cf);
                        carry += cf;
                    }
                    for (; i < out_len - base_len; i++)
                    {
                        out_base[i] = add_half(out_base[i], uint64_t(carry), cf);
                        carry = cf;
                    }
                }
                // out = M + N * BASE ^ 2 + (M + N - K) ^ BASE
                {
                    k_len = std::min(k_len, out_len - base_len);
                    if ((cmp1 > 0) == (cmp2 > 0))
                    {

                        abs_sub_binary(out_base, out_len - base_len, k, k_len, out_base);
                    }
                    else
                    {
                        abs_add_binary_half(out_base, out_len - base_len, k, k_len, out_base);
                    }
                }
            }
            // NTT square
            inline void abs_sqr64_ntt(const uint64_t in[], size_t len, uint64_t out[])
            {
                using namespace hint::transform::ntt;
                if (0 == len || in == nullptr)
                {
                    return;
                }
                size_t out_len = len * 2, conv_len = out_len - 1;
                size_t ntt_len = hint::int_ceil2(conv_len);
                std::vector<NTT0::ModIntType> buffer1(ntt_len);
                {
                    std::copy(in, in + len, buffer1.begin());
                    NTT0::convolutionRecursion(buffer1.data(), buffer1.data(), buffer1.data(), ntt_len);
                };
                std::vector<NTT1::ModIntType> buffer2(ntt_len);
                {
                    std::copy(in, in + len, buffer2.begin());
                    NTT1::convolutionRecursion(buffer2.data(), buffer2.data(), buffer2.data(), ntt_len);
                };
                std::vector<NTT2::ModIntType> buffer3(ntt_len);
                {
                    std::copy(in, in + len, buffer3.begin());
                    NTT2::convolutionRecursion(buffer3.data(), buffer3.data(), buffer3.data(), ntt_len);
                };
                Uint192 carry = 0;
                for (size_t i = 0; i < conv_len; i++)
                {
                    carry += crt3(buffer1[i], buffer2[i], buffer3[i]);
                    out[i] = uint64_t(carry);
                    carry = carry.rShift64();
                }
                out[conv_len] = uint64_t(carry);
            }
            // NTT multiplication
            inline void abs_mul64_ntt(const uint64_t in1[], size_t len1, const uint64_t in2[], size_t len2, uint64_t out[])
            {
                if (0 == len1 || 0 == len2 || in1 == nullptr || in2 == nullptr)
                {
                    return;
                }
                if (in1 == in2)
                {
                    abs_sqr64_ntt(in1, len1, out); // Square
                    return;
                }
                using namespace hint::transform::ntt;
                size_t out_len = len1 + len2, conv_len = out_len - 1;
                size_t ntt_len = hint::int_ceil2(conv_len);
                std::vector<NTT0::ModIntType> buffer1(ntt_len);
                {
                    std::vector<NTT0::ModIntType> buffer2(ntt_len);
                    std::copy(in2, in2 + len2, buffer2.begin());
                    std::copy(in1, in1 + len1, buffer1.begin());
                    NTT0::convolutionRecursion(buffer1.data(), buffer2.data(), buffer1.data(), ntt_len);
                };
                std::vector<NTT1::ModIntType> buffer3(ntt_len);
                {
                    std::vector<NTT1::ModIntType> buffer4(ntt_len);
                    std::copy(in2, in2 + len2, buffer4.begin());
                    std::copy(in1, in1 + len1, buffer3.begin());
                    NTT1::convolutionRecursion(buffer3.data(), buffer4.data(), buffer3.data(), ntt_len);
                };
                std::vector<NTT2::ModIntType> buffer5(ntt_len);
                {
                    std::vector<NTT2::ModIntType> buffer6(ntt_len);
                    std::copy(in2, in2 + len2, buffer6.begin());
                    std::copy(in1, in1 + len1, buffer5.begin());
                    NTT2::convolutionRecursion(buffer5.data(), buffer6.data(), buffer5.data(), ntt_len);
                };
                Uint192 carry = 0;
                for (size_t i = 0; i < conv_len; i++)
                {
                    carry += crt3(buffer1[i], buffer3[i], buffer5[i]);
                    out[i] = uint64_t(carry);
                    carry = carry.rShift64();
                }
                out[conv_len] = uint64_t(carry);
            }

            inline void abs_mul64_balanced(const uint64_t in1[], size_t len1, const uint64_t in2[], size_t len2, uint64_t out[], uint64_t *work_begin = nullptr, uint64_t *work_end = nullptr)
            {
                if (len1 < len2)
                {
                    std::swap(in1, in2);
                    std::swap(len1, len2);
                }
                if (len2 <= 24)
                {
                    abs_mul64_classic(in1, len1, in2, len2, out, work_begin, work_end);
                }
                else if (len2 <= 1536)
                {
                    abs_mul64_karatsuba(in1, len1, in2, len2, out, work_begin, work_end);
                }
                else
                {
                    abs_mul64_ntt(in1, len1, in2, len2, out);
                }
            }

            inline void abs_mul64(const uint64_t in1[], size_t len1, const uint64_t in2[], size_t len2, uint64_t out[], uint64_t *work_begin = nullptr, uint64_t *work_end = nullptr)
            {
                if (len1 < len2)
                {
                    std::swap(in1, in2);
                    std::swap(len1, len2);
                }
                if (len2 <= 24)
                {
                    abs_mul64_balanced(in1, len1, in2, len2, out, work_begin, work_end);
                    return;
                }
                // Get enough work memory
                std::vector<uint64_t> work_mem;
                const size_t work_size = len2 * 3 + len1; // len1 + len2 + len2 * 2,存放结果以及平衡乘积
                if (work_begin + work_size > work_end)
                {
                    work_mem.resize(work_size + len2 * 6); // 为karatsuba做准备
                    work_begin = work_mem.data();
                    work_end = work_begin + work_mem.size();
                }
                else
                {
                    std::fill_n(work_begin, work_size, uint64_t(0));
                }
                auto balance_prod = work_begin, total_prod = balance_prod + len2 * 2;
                size_t rem = len1 % len2, i = len2;
                abs_mul64_balanced(in2, len2, in1, len2, total_prod, work_begin + work_size, work_end);
                for (; i < len1 - rem; i += len2)
                {
                    abs_mul64_balanced(in2, len2, in1 + i, len2, balance_prod, work_begin + work_size, work_end);
                    abs_add_binary_half(balance_prod, len2 * 2, total_prod + i, len2, total_prod + i);
                }
                if (rem > 0)
                {
                    abs_mul64(in2, len2, in1 + i, rem, balance_prod, work_begin + work_size, work_end);
                    abs_add_binary_half(total_prod + i, len2, balance_prod, len2 + rem, total_prod + i);
                }
                std::copy(total_prod, total_prod + len1 + len2, out);
            }
            namespace ssa
            {

                // MOD = 2^(WORD_BIT*block_word)+1
                template <typename WordTy = uint64_t>
                class FermatFFT
                {
                public:
                    using Word = WordTy;
                    static constexpr int WORD_BIT = sizeof(Word) * CHAR_BIT;

                    FermatFFT(size_t new_block_word) : block_word(new_block_word)
                    {
                        assert(block_word >= 1);
                    }

                    size_t blockBit() const
                    {
                        return block_word * WORD_BIT;
                    }

                    size_t getNumLen() const
                    {
                        return block_word + 1;
                    }
                    // len([[0,1,..,words_per_block-1],[...]]) = fft_len
                    void difForward(Word in_out[], size_t fft_len, Word *work_begin = nullptr, Word *work_end = nullptr)
                    {
                        assert(blockBit() % fft_len == 0);
                        assert(is_2pow(fft_len));
                        if (fft_len <= 1)
                        {
                            return;
                        }
                        const size_t num_len = getNumLen();
                        if (2 == fft_len)
                        {
                            transform2(in_out, in_out + num_len);
                            return;
                        }
                        const size_t unit_bit = getUnitBit(fft_len) % (blockBit() * 2);
                        const size_t word_len = fft_len * num_len, half_len = word_len / 2;
                        size_t omega_bit = 0;
                        for (auto x = in_out; x < in_out + half_len; x += num_len)
                        {
                            auto y = x + half_len;
                            transform2(x, y);
                            mul2Pow(y, omega_bit, work_begin, work_end);
                            omega_bit = omegaAddBit(omega_bit, unit_bit);
                        }
                        difForward(in_out, fft_len / 2, work_begin, work_end);
                        difForward(in_out + half_len, fft_len / 2, work_begin, work_end);
                    }

                    void ditBackward(Word in_out[], size_t fft_len, Word *work_begin = nullptr, Word *work_end = nullptr)
                    {
                        assert(blockBit() % fft_len == 0);
                        assert(is_2pow(fft_len));
                        if (fft_len <= 1)
                        {
                            return;
                        }
                        const size_t num_len = getNumLen();
                        if (2 == fft_len)
                        {
                            transform2(in_out, in_out + num_len);
                            return;
                        }
                        const size_t word_len = fft_len * num_len, half_len = word_len / 2;
                        ditBackward(in_out, fft_len / 2, work_begin, work_end);
                        ditBackward(in_out + half_len, fft_len / 2, work_begin, work_end);
                        const size_t unit_bit = getInvUnitBit(fft_len) % (blockBit() * 2);
                        size_t omega_bit = 0;
                        for (auto x = in_out; x < in_out + half_len; x += num_len)
                        {
                            auto y = x + half_len;
                            mul2Pow(y, omega_bit, work_begin, work_end);
                            transform2(x, y);
                            omega_bit = omegaAddBit(omega_bit, unit_bit);
                        }
                    }

                    // unit = 2 ** (word_bit * block_word * 2 / rank)
                    size_t getUnitBit(size_t rank) const
                    {
                        return 2 * blockBit() / rank;
                    }
                    //  inv_unit = 2 ** (word_bit * 2 (block_word - block_word / rank))
                    size_t getInvUnitBit(size_t rank) const
                    {
                        return 2 * (blockBit() - blockBit() / rank);
                    }

                    size_t omegaAddBit(size_t omega_bit, size_t inc) const
                    {
                        omega_bit += inc;
                        return omega_bit < blockBit() * 2 ? omega_bit : omega_bit - blockBit() * 2;
                    }

                    // x * 2 ^ k % MOD
                    void mul2Pow(Word x[], size_t k, Word *work_begin = nullptr, Word *work_end = nullptr) const
                    {
                        // k %= (blockBit() * 2);
                        assert(k < blockBit() * 2);
                        bool sign = k >= blockBit();
                        k = sign ? k - blockBit() : k;
                        if (k == 0)
                        {
                            if (sign)
                            {
                                modSub(x, x);
                            }
                            return;
                        }
                        std::vector<Word> work_mem;
                        const size_t work_size = getNumLen() * 2;
                        if (work_begin + work_size > work_end)
                        {
                            work_mem.resize(work_size);
                            work_begin = work_mem.data();
                            work_end = work_begin + work_mem.size();
                        }
                        else
                        {
                            std::fill(work_begin + getNumLen() + k / WORD_BIT, work_begin + getNumLen() * 2, Word{});
                        }
                        auto add = work_begin, sub = add + block_word;
                        support::lshift_binary(x, getNumLen(), add, k);
                        size_t add_len = block_word, sub_len = getNumLen();
                        if (sign)
                        {
                            std::swap(add, sub);
                            std::swap(add_len, sub_len);
                        }
                        bool borrow = hint::arithm::addition_binary::abs_sub_binary(add, add_len, sub, sub_len, x);
                        if (borrow)
                        {
                            addMod(x);
                        }
                    }

                    void subInMod(const Word x[], const Word y[], Word diff[]) const
                    {
                        bool borrow = arithm::addition_binary::abs_sub_binary(x, getNumLen(), y, getNumLen(), diff);
                        if (borrow)
                        {
                            addMod(diff);
                        }
                    }
                    void addInMod(const Word x[], const Word y[], Word sum[]) const
                    {
                        arithm::addition_binary::abs_add_binary_half(x, getNumLen(), y, getNumLen(), sum);
                        if (compareWithMod(sum) >= 0)
                        {
                            subMod(sum);
                        }
                    }

                    void transform2(Word x[], Word y[]) const
                    {
                        bool carry = false, borrow = false;
                        Word x0, y0, n = 0;
                        {
                            x0 = x[0], y0 = y[0];
                            x[0] = hint::add_half(x0, y0, carry);
                            y[0] = hint::sub_half(x0, y0, borrow);
                        }
                        for (size_t i = 1; i < block_word; i++)
                        {
                            x0 = x[i], y0 = y[i];
                            y[i] = hint::sub_borrow(x0, y0, borrow);
                            x[i] = hint::add_carry(x0, y0, carry);
                            n = std::max(n, x[i]);
                        }
                        {
                            x0 = x[block_word], y0 = y[block_word];
                            x[block_word] = hint::add_carry(x0, y0, carry);
                            y[block_word] = hint::sub_borrow(x0, y0, borrow);
                        }
                        if (x[block_word] > 1)
                        {
                            subMod(x);
                        }
                        else if (x[block_word] == 1)
                        {
                            if (n > 0)
                            {
                                subMod(x);
                            }
                            else if (x[0] > 1)
                            {
                                subMod(x);
                            }
                        }
                        if (borrow)
                        {
                            addMod(y);
                        }
                    }

                    void addMod(Word n[]) const
                    {
                        bool carry = false;
                        n[0] = hint::add_half<Word>(n[0], 1, carry);
                        for (size_t i = 1; i < block_word; i++)
                        {
                            n[i] = hint::add_half<Word>(n[i], carry, carry);
                        }
                        n[block_word] = n[block_word] + 1 + carry;
                    }
                    bool subMod(Word n[]) const
                    {
                        bool borrow = false;
                        n[0] = hint::sub_half<Word>(n[0], 1, borrow);
                        for (size_t i = 1; i < block_word; i++)
                        {
                            n[i] = hint::sub_half<Word>(n[i], Word(borrow), borrow);
                        }
                        n[block_word] = hint::sub_borrow<Word>(n[block_word], 1, borrow);
                        return borrow;
                    }

                    void modSub(const Word in[], Word out[]) const
                    {
                        if (support::count_ture_length(in, getNumLen()) == 0 && in != out)
                        {
                            std::fill_n(out, getNumLen(), Word{});
                        }
                        bool borrow = false;
                        out[0] = hint::sub_half<Word>(1, in[0], borrow);
                        for (size_t i = 1; i < block_word; i++)
                        {
                            Word n0 = in[i];
                            out[i] = Word(0) - borrow - n0;
                            borrow = borrow || n0 > 0;
                        }
                        out[block_word] = hint::sub_borrow<Word>(1, in[block_word], borrow);
                        assert(!borrow);
                    }

                    // n cmp 2^(WORD_BIT * block_word) + 1
                    int compareWithMod(const Word n[]) const
                    {
                        size_t len = getNumLen();
                        assert(len >= 2);
                        if (n[len - 1] != 1)
                        {
                            return n[len - 1] > 1 ? 1 : -1;
                        }
                        for (size_t i = len - 2; i >= 1; i--)
                        {
                            if (n[i] > 0)
                            {
                                return 1;
                            }
                        }
                        if (n[0] != 1)
                        {
                            return n[0] > 1 ? 1 : -1;
                        }
                        return 0;
                    }

                    // n cmp 2^(WORD_BIT * block_word - 1)
                    int compareWithHalfMod(const Word n[]) const
                    {
                        constexpr Word HALF = Word(1) << (WORD_BIT - 1);
                        size_t len = getNumLen();
                        assert(len >= 2);
                        if (n[len - 1] > 0)
                        {
                            return 1;
                        }
                        if (n[len - 2] != HALF)
                        {
                            return n[len - 2] > HALF ? 1 : -1;
                        }
                        size_t i = len - 2;
                        while (i > 0)
                        {
                            i--;
                            if (n[i] > 0)
                            {
                                return 1;
                            }
                        }
                        return 0;
                    }

                    bool addWithSign(const Word in1[], bool sign1, const Word in2[], bool sign2, Word out[])
                    {
                        size_t num_len = getNumLen();
                        if (sign1 == sign2)
                        {
                            bool carry = arithm::addition_binary::abs_add_binary_half(in1, num_len, in2, num_len, out);
                            assert(!carry);
                            return sign1;
                        }
                        int cmp = arithm::addition_binary::abs_difference_binary(in1, num_len, in2, num_len, out);
                        if (0 == cmp)
                        {
                            return false;
                        }
                        return cmp > 0 ? sign1 : sign2;
                    }

                    // 0<=in<MOD, rem = in % 2^(WORD_BIT*block_len), in /= 2^(WORD_BIT*block_len)
                    bool divRemBlockPow(Word in[], bool is_neg, Word rem[], size_t block_len) const
                    {
                        size_t num_len = getNumLen();
                        assert(num_len > block_len);
                        std::copy(in, in + block_len, rem);
                        std::copy(in + block_len, in + num_len, in);
                        std::fill(in + num_len - block_len, in + num_len, Word{});
                        if (is_neg && support::count_ture_length(rem, block_len) > 0)
                        {
                            arithm::addition_binary::binary_complement_code(rem, block_len, rem);
                            arithm::addition_binary::abs_add_num_binary(in, num_len - block_len, Word(1), in);
                        }
                        return is_neg;
                    }

                    size_t maxFFTLen() const
                    {
                        return size_t(1) << hint_ctz(block_word);
                    }

                    size_t nextMinBlockWord(size_t fft_len) const
                    {
                        assert(block_word % fft_len == 0);
                        size_t split_word = block_word / fft_len;
                        return split_word * 2 + 1;
                    }

                    size_t nextFFTLength() const
                    {
                        size_t next_block_word = nextBlockWord();
                        size_t next_fft_len = next_block_word * WORD_BIT;
                        next_fft_len = std::min(maxFFTLen(), next_fft_len);
                        while (next_block_word >= nextMinBlockWord(next_fft_len / 2))
                        {
                            next_fft_len /= 2;
                        }
                        return next_fft_len;
                    }

                    size_t nextBlockWord() const
                    {
                        constexpr int WORD_BIT_LOG = hint_log2(WORD_BIT);
                        size_t fft_len_max = maxFFTLen();
                        int block_len_log = 1;
                        while ((size_t(1) << block_len_log) < nextMinBlockWord(std::min(fft_len_max, size_t(1) << (block_len_log + WORD_BIT_LOG))))
                        {
                            block_len_log++;
                        }
                        return size_t(1) << block_len_log;
                    }

                    size_t nextSplitWord() const
                    {
                        size_t next_fft_len = nextFFTLength();
                        assert(block_word % next_fft_len == 0);
                        return block_word / next_fft_len;
                    }

                    static void splitNumToArr(const Word in[], size_t len, Word out[], size_t num_len, size_t split_word)
                    {
                        assert(num_len >= split_word);
                        size_t rem = len - len % split_word;
                        size_t i = 0, j = 0;
                        for (; i < rem; i += split_word, j += num_len)
                        {
                            std::copy(in + i, in + i + split_word, out + j);
                        }
                        if (len > rem)
                        {
                            std::copy(in + rem, in + len, out + j);
                        }
                    }

                    void mulMod(const Word in1[], const Word in2[], Word out[], Word *work_begin = nullptr, Word *work_end = nullptr)
                    {
                        size_t len1 = support::count_ture_length(in1, getNumLen());
                        size_t len2 = support::count_ture_length(in2, getNumLen());
                        if (0 == len1 || 0 == len2)
                        {
                            return;
                        }
                        if (len1 + len2 <= block_word)
                        {
                            mul(in1, len1, in2, len2, out);
                            std::fill(out + len1 + len2, out + getNumLen(), Word{});
                            return;
                        }
                        if (len1 == getNumLen())
                        {
                            modSub(in2, out);
                            return;
                        }
                        if (len2 == getNumLen())
                        {
                            modSub(in1, out);
                            return;
                        }
                        if (block_word <= 128)
                        {
                            // Get enough work_mem
                            std::vector<Word> work_mem;
                            const size_t work_size = len1 + len2;
                            if (work_begin + work_size > work_end)
                            {
                                work_mem.resize(work_size * 4);
                                work_begin = work_mem.data();
                                work_end = work_begin + work_mem.size();
                            }
                            auto add = work_begin, sub = add + block_word;
                            size_t add_len = block_word, sub_len = work_size - block_word;
                            multiplication::abs_mul64(in1, len1, in2, len2, add, work_begin + work_size, work_end);
                            bool borrow = addition_binary::abs_sub_binary(add, add_len, sub, sub_len, out);
                            out[block_word] = Word(0) - borrow;
                            if (borrow)
                            {
                                addMod(out);
                            }
                            return;
                        }
                        size_t next_block_word = nextBlockWord(), next_fft_len = nextFFTLength(), split_word = nextSplitWord();
                        assert(next_block_word >= split_word * 2 + 1); // 保证MOD足够大而无溢出
                        assert(block_word > next_block_word);          // 保证carry无溢出
                        FermatFFT<Word> next_fft(next_block_word);
                        size_t num_len = next_fft.getNumLen(), word_len = num_len * next_fft_len;
                        std::vector<Word> work_mem;
                        const size_t work_size = word_len * 2 + num_len + getNumLen();
                        if (work_begin + work_size > work_end)
                        {
                            work_mem.resize(work_size * 2);
                            work_begin = work_mem.data();
                            work_end = work_begin + work_mem.size();
                        }
                        auto fft1 = work_begin, fft2 = work_begin + word_len, carry = fft2 + word_len, sub = carry + num_len;
                        splitNumToArr(in1, len1, fft1, num_len, split_word);
                        splitNumToArr(in2, len2, fft2, num_len, split_word);

                        next_fft.negaCyclicConvolution(fft1, fft2, next_fft_len, work_begin + work_size, work_end);

                        bool carry_sign = false;
                        for (size_t i = 0, j = 0; i < block_word; i += split_word, j += num_len)
                        {
                            auto fft_ele = fft1 + j;
                            bool fft_sign = false;
                            if (next_fft.compareWithHalfMod(fft_ele) >= 0)
                            {
                                next_fft.modSub(fft_ele, fft_ele);
                                fft_sign = true;
                            }
                            carry_sign = next_fft.addWithSign(fft_ele, fft_sign, carry, carry_sign, carry);
                            next_fft.divRemBlockPow(carry, carry_sign, out + i, split_word);
                        }
                        for (size_t i = 0; i < block_word; i += split_word)
                        {
                            next_fft.divRemBlockPow(carry, carry_sign, sub + i, split_word);
                        }
                        if (carry_sign)
                        {
                            arithm::addition_binary::abs_add_num_binary_half(sub, getNumLen(), carry[0], sub);
                        }
                        else
                        {
                            arithm::addition_binary::abs_add_num_binary_half(out, getNumLen(), carry[0], out);
                        }
                        subInMod(out, sub, out);
                    }

                    void mul(const Word in1[], size_t len1, const Word in2[], size_t len2, Word out[], Word *work_begin = nullptr, Word *work_end = nullptr)
                    {
                        assert(len1 + len2 <= block_word);
                        if (len1 < len2)
                        {
                            std::swap(len1, len2);
                            std::swap(in1, in2);
                        }
                        if (0 == len2)
                        {
                            return;
                        }
                        if (len2 <= 128)
                        {
                            multiplication::abs_mul64(in1, len1, in2, len2, out);
                            return;
                        }
                        size_t next_block_word = nextBlockWord(), split_word = nextSplitWord();
                        size_t in1_split_block = (len1 + split_word - 1) / split_word;
                        size_t in2_split_block = (len2 + split_word - 1) / split_word;
                        size_t res_split_block = in1_split_block + in2_split_block, next_fft_len = int_ceil2(res_split_block - 1);
                        size_t conv_word = (res_split_block - 1) * split_word;

                        assert(next_block_word >= split_word * 2 + 1); // 保证MOD足够大而无溢出
                        assert(block_word > next_block_word);          // 保证carry无溢出

                        FermatFFT<Word> next_fft(next_block_word);
                        size_t num_len = next_fft.getNumLen(), word_len = num_len * next_fft_len;

                        std::vector<Word> work_mem;
                        const size_t work_size = word_len * 2 + num_len;
                        if (work_begin + work_size > work_end)
                        {
                            work_mem.resize(work_size * 2);
                            work_begin = work_mem.data();
                            work_end = work_begin + work_mem.size();
                        }
                        auto fft1 = work_begin, fft2 = fft1 + word_len, carry = fft2 + word_len;

                        bool carry_sign = false;
                        splitNumToArr(in1, len1, fft1, num_len, split_word);
                        splitNumToArr(in2, len2, fft2, num_len, split_word);

                        next_fft.cyclicConvolution(fft1, fft2, next_fft_len, work_begin + work_size, work_end);
                        size_t i = 0, j = 0;
                        for (; i < conv_word - split_word; i += split_word, j += num_len)
                        {
                            auto fft_ele = fft1 + j;
                            next_fft.addInMod(fft_ele, carry, carry);
                            next_fft.divRemBlockPow(carry, carry_sign, out + i, split_word);
                        }
                        next_fft.addInMod(fft1 + j, carry, carry);
                        std::copy(carry, carry + len1 + len2 - i, out + i);
                    }

                    void cyclicConvolution(Word in1_out[], Word in2[], size_t fft_len, Word *work_begin = nullptr, Word *work_end = nullptr)
                    {
                        assert(is_2pow(fft_len));
                        size_t len_inv_bit = blockBit() * 2 - hint_log2(fft_len);
                        len_inv_bit %= (blockBit() * 2);
                        difForward(in1_out, fft_len, work_begin, work_end);
                        difForward(in2, fft_len, work_begin, work_end);
                        const size_t num_len = getNumLen(), word_len = fft_len * num_len;
                        for (size_t i = 0; i < word_len; i += num_len)
                        {
                            mulMod(in1_out + i, in2 + i, in1_out + i, work_begin, work_end);
                            mul2Pow(in1_out + i, len_inv_bit, work_begin, work_end);
                        }
                        ditBackward(in1_out, fft_len, work_begin, work_end);
                    }

                    void weight(Word in_out[], size_t theta_bit, size_t fft_len, Word *work_begin = nullptr, Word *work_end = nullptr)
                    {
                        size_t theta = 0;
                        theta_bit %= (blockBit() * 2);
                        const size_t num_len = getNumLen(), word_len = fft_len * num_len;
                        for (size_t i = 0; i < word_len; i += num_len)
                        {
                            mul2Pow(in_out + i, theta, work_begin, work_end);
                            theta = omegaAddBit(theta, theta_bit);
                        }
                    }

                    void negaCyclicConvolution(Word in1_out[], Word in2[], size_t fft_len, Word *work_begin = nullptr, Word *work_end = nullptr)
                    {
                        assert(is_2pow(fft_len));
                        const size_t theta_bit = blockBit() / fft_len;
                        const size_t inv_theta_bit = blockBit() * 2 - theta_bit;
                        weight(in1_out, theta_bit, fft_len, work_begin, work_end);
                        difForward(in1_out, fft_len, work_begin, work_end);
                        weight(in2, theta_bit, fft_len, work_begin, work_end);
                        difForward(in2, fft_len, work_begin, work_end);
                        const size_t num_len = getNumLen(), word_len = fft_len * num_len;
                        for (size_t i = 0; i < word_len; i += num_len)
                        {
                            mulMod(in1_out + i, in2 + i, in1_out + i, work_begin, work_end);
                        }
                        ditBackward(in1_out, fft_len, work_begin, work_end);
                        size_t len_inv_bit = blockBit() * 2 - hint_log2(fft_len);
                        weight(in1_out, inv_theta_bit + len_inv_bit, fft_len, work_begin, work_end);
                    }

                private:
                    size_t block_word;
                };
            }
            inline void abs_mul64_ssa(const uint64_t in1[], size_t len1, const uint64_t in2[], size_t len2, uint64_t out[], uint64_t *work_begin = nullptr, uint64_t *work_end = nullptr)
            {
                ssa::FermatFFT<uint64_t> fft(int_ceil2(len1 + len2)); // TODO OPTI
                fft.mul(in1, len1, in2, len2, out);
            }
        }
        namespace division
        {
            using namespace multiplication;
            using namespace addition_binary;

            /// @brief 2^64 base long integer divided by 64bit number.
            /// @param in Input long integer.
            /// @param len Length of input long integer.
            /// @param out Output long integer, equals to input / divisor.
            /// @param divisor the 64 bit number to divide by.
            /// @return The remainder of the division, equals to input % divisor.
            inline uint64_t abs_div_rem_num64(const uint64_t in[], size_t len, uint64_t out[], uint64_t divisor)
            {
                uint64_t hi64 = 0;
                const DivSupporter<uint64_t, Uint128> divisor_supporter(divisor);
                while (len > 0)
                {
                    len--;
                    Uint128 n(in[len], hi64);
                    out[len] = divisor_supporter.divMod(n);
                    hi64 = uint64_t(divisor);
                }
                return hi64;
            }

            template <typename T>
            inline T divisor_normalize(T divisor[], size_t len, T base)
            {
                if (0 == len || nullptr == divisor)
                {
                    return 0;
                }
                if (divisor[len - 1] >= base / 2)
                {
                    return 1;
                }
                T first = divisor[len - 1];
                T factor = 1;
                while (first < base / 2)
                {
                    factor *= 2;
                    first *= 2;
                }
                const DivSupporter<uint64_t, Uint128> base_supporter(base);
                uint64_t carry = 0;
                for (size_t i = 0; i < len; i++)
                {
                    Uint128 temp = Uint128(divisor[i]) * factor + carry;
                    carry = base_supporter.divMod(temp);
                    divisor[i] = T(temp);
                }
                assert(carry == 0);
                assert(divisor[len - 1] >= base / 2);
                return factor;
            }

            constexpr int divisor_normalize64(const uint64_t in[], size_t len, uint64_t out[])
            {
                constexpr uint64_t BASE_HALF = uint64_t(1) << 63;
                if (0 == len || nullptr == in)
                {
                    return 0;
                }
                if (in[len - 1] >= BASE_HALF)
                {
                    std::copy(in, in + len, out);
                    return 0;
                }
                uint64_t first = in[len - 1];
                assert(first > 0);
                const int leading_zeros = hint_clz(first);
                lshift_in_word_half(in, len, out, leading_zeros);
                assert(out[len - 1] >= BASE_HALF);
                return leading_zeros;
            }

            // 只处理商的长度为dividend_len-divisor_len的情况
            inline void abs_div64_classic_core(uint64_t dividend[], size_t dividend_len, const uint64_t divisor[], size_t divisor_len, uint64_t quotient[], uint64_t *work_begin = nullptr, uint64_t *work_end = nullptr)
            {
                if (nullptr == dividend || dividend_len <= divisor_len)
                {
                    return;
                }
                assert(divisor[divisor_len - 1] >= uint64_t(1) << 63); // 除数最高位为1
                assert(divisor_len > 0);
                size_t pre_dividend_len = dividend_len;
                remove_leading_zeros(dividend, dividend_len); // 去除前导0
                while (dividend_len < pre_dividend_len && dividend_len >= divisor_len)
                {
                    size_t quot_i = dividend_len - divisor_len;
                    if (abs_compare(dividend + quot_i, divisor_len, divisor, divisor_len) >= 0)
                    {
                        quotient[quot_i] = 1;
                        abs_sub_binary(dividend + quot_i, divisor_len, divisor, divisor_len, dividend + quot_i);
                    }
                    pre_dividend_len = dividend_len;
                    remove_leading_zeros(dividend, dividend_len); // 去除前导0
                }
                if (dividend_len < divisor_len)
                {
                    return;
                }
                if (1 == divisor_len)
                {
                    uint64_t rem = abs_div_rem_num64(dividend, dividend_len, quotient, divisor[0]);
                    dividend[0] = rem;
                    return;
                }
                // Get enough work memory
                std::vector<uint64_t> work_mem;
                const size_t work_size = divisor_len + 1;
                if (work_begin + work_size > work_end)
                {
                    work_mem.resize(work_size);
                    work_begin = work_mem.data();
                    work_end = work_begin + work_mem.size();
                }
                const uint64_t divisor_1 = divisor[divisor_len - 1];
                const uint64_t divisor_0 = divisor[divisor_len - 2];
                const DivSupporter<uint64_t, Uint128> div(divisor_1);
                size_t i = dividend_len - divisor_len;
                while (i > 0)
                {
                    i--;
                    uint64_t qhat = 0, rhat = 0;
                    const uint64_t dividend_2 = dividend[divisor_len + i];
                    const uint64_t dividend_1 = dividend[divisor_len + i - 1];
                    assert(dividend_2 <= divisor_1);
                    Uint128 dividend_num = Uint128(dividend_1, dividend_2);
                    if (dividend_2 == divisor_1)
                    {
                        qhat = UINT64_MAX;
                        rhat = uint64_t(dividend_num - Uint128(divisor_1) * qhat);
                    }
                    else
                    {
                        qhat = div.divMod(dividend_num);
                        rhat = uint64_t(dividend_num);
                    }
                    { // 3 words / 2 words refine
                        dividend_num = Uint128(dividend[divisor_len + i - 2], rhat);
                        Uint128 prod = Uint128(divisor_0) * qhat;
                        if (prod > dividend_num)
                        {
                            qhat--;
                            prod -= divisor_0;
                            dividend_num += Uint128(0, divisor_1);
                            if (dividend_num > Uint128(0, divisor_1) && prod > dividend_num) [[unlikely]]
                            {
                                qhat--;
                            }
                        }
                    }
                    if (qhat > 0)
                    {
                        auto prod = work_begin;
                        auto rem = dividend + i;
                        abs_mul_add_num64(divisor, divisor_len, prod, 0, qhat);
                        size_t len_prod = count_ture_length(prod, divisor_len + 1);
                        size_t len_rem = count_ture_length(rem, divisor_len + 1);
                        int count = 0;
                        while (abs_compare(prod, len_prod, rem, len_rem) > 0)
                        {
                            qhat--;
                            abs_sub_binary(prod, len_prod, divisor, divisor_len, prod);
                            remove_leading_zeros(prod, len_prod);
                            assert(count < 2);
                            count++;
                        }
                        if (qhat > 0)
                        {
                            abs_sub_binary(rem, len_rem, prod, len_prod, rem);
                        }
                    }
                    quotient[i] = qhat;
                }
            }

            inline void div_check(const uint64_t in1[], size_t len1, const uint64_t in2[], size_t len2, const uint64_t quot[], const uint64_t rem[])
            {
                std::vector<uint64_t> quot_v(len1 - len2), rem_v(in1, in1 + len1);
                abs_div64_classic_core(rem_v.data(), len1, in2, len2, quot_v.data());
                for (size_t i = 0; i < len1 - len2; i++)
                {
                    if (quot_v[i] != quot[i])
                    {
                        std::cout << "quot check fail" << std::endl;
                        exit(1);
                    }
                }
                std::cout << "quot check success" << std::endl;
                for (size_t i = 0; i < len2; i++)
                {
                    if (rem_v[i] != rem[i])
                    {
                        std::cout << "rem check fail" << std::endl;
                        exit(1);
                    }
                }
                std::cout << "rem check success" << std::endl;
            }

            // 只处理商的长度为dividend_len-divisor_len的情况
            inline void abs_div64_recursive_core(uint64_t dividend[], size_t dividend_len, const uint64_t divisor[], size_t divisor_len, uint64_t quotient[], uint64_t *work_begin = nullptr, uint64_t *work_end = nullptr)
            {
                if (nullptr == dividend || dividend_len <= divisor_len)
                {
                    return;
                }
                assert(divisor_len > 0);
                assert(divisor[divisor_len - 1] >= uint64_t(1) << 63); // 除数最高位为1
                if (divisor_len <= 32 || (dividend_len <= divisor_len + 16))
                {
                    abs_div64_classic_core(dividend, dividend_len, divisor, divisor_len, quotient, work_begin, work_end);
                    return;
                }
                // 被除数分段处理
                if (dividend_len > divisor_len * 2)
                {
                    size_t rem = dividend_len % divisor_len, shift = dividend_len - divisor_len - rem;
                    if (rem > 0)
                    {
                        abs_div64_recursive_core(dividend + shift, divisor_len + rem, divisor, divisor_len, quotient + shift, work_begin, work_end);
                    }
                    while (shift > 0)
                    {
                        shift -= divisor_len;
                        abs_div64_recursive_core(dividend + shift, divisor_len * 2, divisor, divisor_len, quotient + shift, work_begin, work_end);
                    }
                }
                else if (dividend_len < divisor_len * 2)
                {
                    assert(dividend_len > divisor_len);
                    // 进行估商处理
                    size_t shift_len = divisor_len * 2 - dividend_len; // 进行截断处理,使得截断后被除数的长度为除数的两倍,以进行估商
                    size_t next_divisor_len = divisor_len - shift_len;
                    size_t next_dividend_len = dividend_len - shift_len;
                    size_t quot_len = next_dividend_len - next_divisor_len;
                    // Get enough work memory
                    std::vector<uint64_t> work_mem;
                    const size_t work_size = quot_len + shift_len;
                    if (work_begin + work_size > work_end)
                    {
                        work_mem.resize(work_size * 3);
                        work_begin = work_mem.data();
                        work_end = work_begin + work_mem.size();
                    }
                    auto prod = work_begin;
                    abs_div64_recursive_core(dividend + shift_len, next_dividend_len, divisor + shift_len, next_divisor_len, quotient, work_begin, work_end);
                    abs_mul64(divisor, shift_len, quotient, quot_len, prod, work_begin + work_size, work_end);
                    size_t prod_len = count_ture_length(prod, quot_len + shift_len), rem_len = count_ture_length(dividend, divisor_len);
                    // 修正
                    int count = 0;
                    while (abs_compare(prod, prod_len, dividend, rem_len) > 0)
                    {
                        abs_sub_num_binary(quotient, quot_len, uint64_t(1), quotient);
                        abs_sub_binary(prod, prod_len, divisor, shift_len, prod);
                        bool carry = abs_add_binary_half(dividend + shift_len, rem_len - shift_len, divisor + shift_len, quot_len, dividend + shift_len);
                        if (carry)
                        {
                            if (rem_len < dividend_len) // 防止溢出
                            {
                                dividend[rem_len] = 1;
                                rem_len++;
                            }
                            else
                            {
                                // 说明此时rem_len = dividend_len + 1
                                // prod_len <= quot_len + shift_len = dividend_len - divisor_len + divisor_len * 2 - dividend_len = divisor_len <= dividend_len
                                // 故此时prod_len < rem_len
                                break;
                            }
                        }
                        remove_leading_zeros(dividend, rem_len);
                        remove_leading_zeros(prod, prod_len);
                        assert(count < 2);
                        count++;
                    }
                    abs_sub_binary(dividend, rem_len, prod, prod_len, dividend);
                }
                else
                {
                    // 进行两次递归处理
                    size_t quot_lo_len = divisor_len / 2, quot_hi_len = divisor_len - quot_lo_len;
                    size_t dividend_len1 = divisor_len + quot_hi_len, dividend_len2 = divisor_len + quot_lo_len;
                    abs_div64_recursive_core(dividend + quot_lo_len, dividend_len1, divisor, divisor_len, quotient + quot_lo_len, work_begin, work_end); // 求出高位的商
                    abs_div64_recursive_core(dividend, dividend_len2, divisor, divisor_len, quotient, work_begin, work_end);                             // 求出低位的商
                }
            }

            // BASE = 2^64,BASE ^ base_len / in
            inline void abs_reciprocal64(const uint64_t in[], size_t in_len, size_t base_len, uint64_t out[])
            {
            }

            inline void abs_div64(uint64_t dividend[], size_t dividend_len, const uint64_t divisor[], size_t divisor_len, uint64_t quotient[])
            {
                std::fill_n(quotient, get_div_len(dividend_len, divisor_len), uint64_t(0));
                remove_leading_zeros(dividend, dividend_len);
                remove_leading_zeros(divisor, divisor_len);
                if (dividend_len < divisor_len)
                {
                    return;
                }
                size_t norm_dividend_len = dividend_len + 1, norm_divisor_len = divisor_len;
                std::vector<uint64_t> dividend_v(norm_dividend_len), divisor_v(norm_divisor_len);

                auto norm_dividend = dividend_v.data(), norm_divisor = divisor_v.data();
                int leading_zero = divisor_normalize64(divisor, divisor_len, norm_divisor);
                lshift_in_word(dividend, dividend_len, norm_dividend, leading_zero);
                remove_leading_zeros(norm_dividend, norm_dividend_len);
                // 解决商最高位为1的情况
                if (norm_dividend_len == dividend_len)
                {
                    size_t quot_len = norm_dividend_len - norm_divisor_len;
                    if (abs_compare(norm_dividend + quot_len, norm_divisor_len, norm_divisor, norm_divisor_len) >= 0)
                    {
                        quotient[quot_len] = 1;
                        abs_sub_binary(norm_dividend + quot_len, norm_divisor_len, norm_divisor, norm_divisor_len, norm_dividend + quot_len);
                    }
                }
                // 剩余的商长度为quot_len
                abs_div64_recursive_core(norm_dividend, norm_dividend_len, norm_divisor, norm_divisor_len, quotient);

                remove_leading_zeros(norm_dividend, norm_divisor_len); // 余数在被除数上
                rshift_in_word(norm_dividend, norm_divisor_len, dividend, leading_zero);
                std::fill(dividend + norm_divisor_len, dividend + dividend_len, uint64_t(0));
            }
        }
    }
}
#endif