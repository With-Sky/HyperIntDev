#include "hint_decimal.hpp"
#include <chrono>
inline uint64_t abs_mul_add_num64_half(const uint64_t in[], size_t len, uint64_t out[], uint64_t num_add, uint64_t num_mul)
{
    using namespace hint::extend_int;
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
// in * num_mul + in_out -> in_out
inline void mul64_sub_proc(const uint64_t in[], size_t len, uint64_t in_out[], uint64_t num_mul)
{
    using namespace hint::extend_int;
    uint64_t carry = 0;
    size_t i = 0;
    for (const size_t rem_len = len - len % 4; i < rem_len; i += 4)
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
// 小学乘法
inline void abs_mul64_classic(const uint64_t in1[], size_t len1, const uint64_t in2[], size_t len2, uint64_t out[], uint64_t *work_begin = nullptr, uint64_t *work_end = nullptr)
{
    using namespace hint::arithm::support;
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

inline void testbase()
{
    using namespace hint::arithm::multiplication;
    constexpr uint64_t BASE = 1e18;
    constexpr BaseExecutor<uint64_t> exec(BASE);
    size_t len1 = 1 << 15, len2 = 1 << 15;
    std::vector<uint64_t> a(len1, BASE - 1), b(len2, BASE - 1), c(len1 + len2, 0), buffer(len1 + len2, 0);
    auto t1 = std::chrono::high_resolution_clock::now();
    abs_mul_classic(a.data(), len1, b.data(), len2, c.data(), exec, buffer.data(), buffer.data() + len1 + len2);
    auto t2 = std::chrono::high_resolution_clock::now();

    // hint::arithm::support::ary_print(c.data(), len1 + len2);
    std::cout << "Time: " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() << "us" << std::endl;
}

inline auto testbin()
{
    size_t len1 = 1 << 15, len2 = 1 << 15;
    std::vector<uint64_t> a(len1, UINT64_MAX), b(len2, UINT64_MAX), c(len1 + len2, 0), buffer(len1 + len2, 0);
    auto t1 = std::chrono::high_resolution_clock::now();
    abs_mul64_classic(a.data(), len1, b.data(), len2, c.data(), buffer.data(), buffer.data() + len1 + len2);
    auto t2 = std::chrono::high_resolution_clock::now();

    // hint::arithm::support::ary_print(c.data(), len1 + len2);
    std::cout << "Time: " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() << "us" << std::endl;
    return c;
}

inline auto testbin2()
{
    using namespace hint::arithm::multiplication;
    constexpr BaseExecutorBinary exec;
    size_t len1 = 1 << 15, len2 = 1 << 15;
    std::vector<uint64_t> a(len1, UINT64_MAX), b(len2, UINT64_MAX), c(len1 + len2, 0), buffer(len1 + len2, 0);
    auto t1 = std::chrono::high_resolution_clock::now();
    abs_mul_classic(a.data(), len1, b.data(), len2, c.data(), exec, buffer.data(), buffer.data() + len1 + len2);
    auto t2 = std::chrono::high_resolution_clock::now();

    // hint::arithm::support::ary_print(c.data(), len1 + len2);
    std::cout << "Time: " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() << "us" << std::endl;
    return c;
}

inline void perf_mul_basic()
{
    testbase();
    auto bin1 = testbin();
    auto bin2 = testbin2();
    std::cout << std::boolalpha << (bin1 == bin2) << std::endl;
}

inline auto test_add()
{
    using namespace hint::arithm::addition_binary;
    size_t len1 = 1 << 25, len2 = 1 << 25;
    std::vector<uint64_t> a(len1, 0), b(len2, UINT64_MAX), c(std::max(len1, len2), 0);
    auto t1 = std::chrono::high_resolution_clock::now();
    abs_sub_binary(a.data(), len1, b.data(), len2, c.data());
    auto t2 = std::chrono::high_resolution_clock::now();

    // hint::arithm::support::ary_print(c.data(), c.size(), true);
    std::cout << "Time: " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() << "us" << std::endl;
    return c;
}

inline auto test_add_binbase()
{
    using namespace hint::arithm::addition_binary;
    constexpr BaseExecutorBinary exec;
    size_t len1 = 1 << 25, len2 = 1 << 25;
    std::vector<uint64_t> a(len1, 0), b(len2, UINT64_MAX), c(std::max(len1, len2), 0);
    auto t1 = std::chrono::high_resolution_clock::now();
    abs_sub_long_long(a.data(), len1, b.data(), len2, c.data(), exec);
    auto t2 = std::chrono::high_resolution_clock::now();

    // hint::arithm::support::ary_print(c.data(), c.size(), true);
    std::cout << "Time: " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() << "us" << std::endl;
    return c;
}

inline auto test_add_base()
{
    using namespace hint::arithm::addition_binary;
    constexpr uint64_t BASE = 1e19;

    constexpr BaseExecutor<uint64_t> exec(BASE);
    size_t len1 = 1 << 25, len2 = len1;
    std::vector<uint64_t> a(len1, BASE - 1), b(len2, BASE - 1), c(std::max(len1, len2), 0);
    auto t1 = std::chrono::high_resolution_clock::now();
    abs_sub_long_long(a.data(), len1, b.data(), len2, c.data(), exec);
    auto t2 = std::chrono::high_resolution_clock::now();

    // hint::arithm::support::ary_print(c.data(), c.size(), true);
    std::cout << "Time: " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() << "us" << std::endl;
    return c;
}

void perf_add()
{
    test_add_base();
    auto bin2 = test_add_binbase();
    auto bin1 = test_add();
    std::cout << std::boolalpha << (bin1 == bin2) << std::endl;
}

int main()
{
    std::ios::sync_with_stdio(false);
    std::cin.tie(0);
    // perf_mul_basic();
    perf_add();
    return 0;
}