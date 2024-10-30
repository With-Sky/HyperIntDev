#include "hint_arithm.hpp"
#include <bitset>
#include <algorithm>
void test_mul64()
{
    int shift = 0;
    // std::cin >> shift;
    size_t len = 1 << shift;
    size_t len1 = 64, len2 = len1;
    std::vector<uint64_t> a(len1, UINT64_MAX), b(len2, UINT64_MAX);
    static std::mt19937 gen(1);
    std::uniform_int_distribution<uint64_t> dist(0, UINT64_MAX);
    // for (auto &&i : a)
    // {
    //     i = dist(gen);
    // }
    // for (auto &&i : b)
    // {
    //     i = dist(gen);
    // }
    // b = a;
    // a = hint::arithm::support::load_arr("in1.txt");
    // b = hint::arithm::support::load_arr("in2.txt");
    std::vector<uint64_t> c(a.size() + b.size()), d = c, e = c, f = c;
    auto t1 = std::chrono::steady_clock::now();
    hint::arithm::multiplication::abs_mul64_balanced(a.data(), a.size(), b.data(), b.size(), c.data());
    auto t2 = std::chrono::steady_clock::now();
    hint::arithm::multiplication::abs_mul64_ntt(a.data(), a.size(), b.data(), b.size(), d.data());
    auto t3 = std::chrono::steady_clock::now();
    hint::arithm::multiplication::abs_mul64_ssa(a.data(), a.size(), b.data(), b.size(), e.data());
    auto t4 = std::chrono::steady_clock::now();
    hint::arithm::multiplication::abs_mul64_karatsuba(a.data(), a.size(), b.data(), b.size(), f.data());
    auto t5 = std::chrono::steady_clock::now();

    // hint::arithm::support::ary_print(d.data(), d.size());
    // hint::arithm::support::ary_print(e.data(), e.size());

    std::cout << (len1 + len2) * 64 << "bits\n";
    std::cout << std::boolalpha << (e == c) << "\n";
    std::cout << std::boolalpha << (c == d) << "\n";
    std::cout << std::boolalpha << (d == f) << "\n";

    std::cout << std::dec << "balance time:" << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() << "us\n";
    std::cout << std::dec << "ntt time:" << std::chrono::duration_cast<std::chrono::microseconds>(t3 - t2).count() << "us\n";
    std::cout << std::dec << "ssa time:" << std::chrono::duration_cast<std::chrono::microseconds>(t4 - t3).count() << "us\n";
    std::cout << std::dec << "karatsuba time:" << std::chrono::duration_cast<std::chrono::microseconds>(t5 - t4).count() << "us\n";
}

inline void test_add_VVW()
{
    size_t len1 = 200000000, len2 = len1 + 1;
    std::vector<uint64_t> a(len1, UINT64_MAX), b(len1, UINT64_MAX);
    b.resize(len2);
    static std::mt19937 gen(1);
    std::uniform_int_distribution<uint64_t> dist(0, UINT64_MAX);
    for (auto &&i : a)
    {
        // i = dist(gen);
    }
    for (auto &&i : b)
    {
        // i = dist(gen);
    }
    // b = a;
    // a = hint::arithm::support::load_arr("in1.txt");
    // b = hint::arithm::support::load_arr("in2.txt");
    auto t1 = std::chrono::steady_clock::now();
    hint::arithm::multiplication::mul64_sub_proc(a.data(), a.size(), b.data(), UINT64_MAX);
    auto t2 = std::chrono::steady_clock::now();

    // hint::arithm::support::ary_print(b.data(), b.size());

    std::cout << std::dec << "balance time:" << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() << "us\n";
}

inline void test_shift()
{
    constexpr size_t len = 10;
    uint64_t arr[len]{0b11001111110000011, 0b11111111111111110000, 0b0000000000111};
    uint64_t arr1[len]{};
    int lead = 64;
    std::cout << lead << "\n";
    hint::arithm::division::lshift_in_word(arr, len, arr1, lead);
    hint::arithm::division::rshift_in_word(arr1, len, arr1, lead);
    for (int64_t i = 2; i >= 0; i--)
    {
        std::cout << std::bitset<64>(arr[i]) << "\t";
    }
    std::cout << "\n";
    for (int64_t i = 2; i >= 0; i--)
    {
        std::cout << std::bitset<64>(arr1[i]) << "\t";
    }
}
namespace simple_bint
{
    inline uint32_t char_to_num(char c)
    {
        if (c >= 'a' && c <= 'z')
        {
            return c - 'a' + 36;
        }
        else if (c >= 'A' && c <= 'Z')
        {
            return c - 'A' + 10;
        }
        return c - '0';
    }

    inline char num_to_char(uint32_t n)
    {
        if (n >= 36 && n <= 61)
        {
            return n - 36 + 'a';
        }
        else if (n >= 10 && n <= 35)
        {
            return n - 10 + 'A';
        }
        return n + '0';
    }

    using Bint = std::vector<uint64_t>;

    Bint from_string(const std::string &s, int base = 10)
    {
        const int digit_bits = std::ceil(std::log2(base));
        size_t total_bits = s.size() * digit_bits;
        size_t total_words = (total_bits + 63) / 64;
        Bint res(total_words + 1);
        for (auto &&i : s)
        {
            hint::arithm::multiplication::abs_mul_add_num64(res.data(), total_words, res.data(), char_to_num(i), base);
        }
        hint::arithm::support::remove_leading_zeros(res.data(), total_words);
        res.resize(total_words);
        return res;
    }

    std::string to_string(Bint b, int base = 10)
    {
        if (b.empty())
        {
            return "0";
        }
        std::string s;
        size_t word_len = hint::arithm::support::count_ture_length(b.data(), b.size());
        while (word_len > 0)
        {
            uint32_t rem = hint::arithm::division::abs_div_rem_num64(b.data(), word_len, b.data(), base);
            s += num_to_char(rem);
            hint::arithm::support::remove_leading_zeros(b.data(), word_len);
        }
        std::reverse(s.begin(), s.end());
        return s;
    }
} // namespace simple_bint
#include <fstream>

void test_div()
{
    using namespace hint::arithm::support;
    using namespace simple_bint;
    size_t len1 = 40000, len2 = 20000;
    Bint b1(len1, UINT64_MAX), b2(len2);
    b2[len2 - 1] = UINT64_MAX;
    static std::mt19937 gen(1);
    std::uniform_int_distribution<uint64_t> dist(0, UINT64_MAX);
    for (auto &&i : b1)
    {
        i = dist(gen);
    }
    for (auto &&i : b2)
    {
        i = dist(gen);
    }
    // hint::arithm::support::ary_print(b1.data(), b1.size());
    // hint::arithm::support::ary_print(b2.data(), b2.size());
    Bint q(hint::arithm::support::get_div_len(b1.size(), b2.size()));
    auto dividend_copy = b1;

    std::vector<uint64_t> div_copy(b2.size() + q.size());
    auto t1 = std::chrono::high_resolution_clock::now();
    hint::arithm::division::abs_div64(b1.data(), b1.size(), b2.data(), b2.size(), q.data());
    auto t2 = std::chrono::high_resolution_clock::now();
    hint::arithm::multiplication::abs_mul64(q.data(), q.size(), b2.data(), b2.size(), div_copy.data());
    auto t3 = std::chrono::high_resolution_clock::now();
    assert(div_copy.size() >= b2.size());
    hint::arithm::addition_binary::abs_add_binary_half(div_copy.data(), div_copy.size(), b1.data(), b2.size(), div_copy.data());

    b1.resize(b2.size());
    // ary_print(q.data(), q.size());
    // ary_print(b1.data(), b1.size());

    // std::fstream f;
    // f.open("div.txt", std::ios::out);
    // f << to_string(q);
    // f.close();
    // f.open("rem.txt", std::ios::out);
    // f << to_string(b1);
    // f.close();
    div_copy.resize(count_ture_length(div_copy.data(), div_copy.size()));
    dividend_copy.resize(count_ture_length(dividend_copy.data(), dividend_copy.size()));
    std::cout << abs_compare(b2.data(), b2.size(), b1.data(), b1.size()) << std::endl;
    std::cout << abs_compare(div_copy.data(), div_copy.size(), dividend_copy.data(), dividend_copy.size()) << std::endl;
    std::cout << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() << "us" << std::endl;
    std::cout << std::chrono::duration_cast<std::chrono::microseconds>(t3 - t2).count() << "us" << std::endl;
}

template <uint64_t ROOT, typename ModInt>
void ntt_dit(ModInt in_out[], size_t ntt_len)
{
    for (size_t rank = 2; rank <= ntt_len; rank *= 2)
    {
        ModInt unit_omega = hint::qpow(ModInt(ROOT), (ModInt::mod() - 1) / rank);
        size_t dis = rank / 2;
        for (auto begin = in_out; begin < in_out + ntt_len; begin += rank)
        {
            ModInt omega = 1;
            for (auto p = begin; p < begin + dis; p++)
            {
                auto temp0 = p[0], temp1 = p[dis] * omega;
                p[0] = temp0 + temp1;
                p[dis] = temp0 - temp1;
                omega = omega * unit_omega;
            }
        }
    }
}

template <uint64_t ROOT, typename ModInt>
void ntt_dif(ModInt in_out[], size_t ntt_len)
{
    for (size_t rank = ntt_len; rank >= 2; rank /= 2)
    {
        ModInt unit_omega = hint::qpow(ModInt(ROOT), (ModInt::mod() - 1) / rank);
        size_t dis = rank / 2;
        for (auto begin = in_out; begin < in_out + ntt_len; begin += rank)
        {
            ModInt omega = 1;
            for (auto p = begin; p < begin + dis; p++)
            {
                auto temp0 = p[0], temp1 = p[dis];
                p[0] = temp0 + temp1;
                p[dis] = (temp0 - temp1) * omega;
                omega = omega * unit_omega;
            }
        }
    }
}

void ntt_check(int s = 24)
{
    using namespace hint;
    using namespace transform;
    using namespace ntt;
    constexpr uint64_t mod1 = 1945555039024054273, root1 = 5;
    constexpr uint64_t mod2 = 4179340454199820289, root2 = 3;
    constexpr uint64_t mod = 998244353, root = 3;

    // using ModInt = MontInt32Lazy<mod>;
    using ntt = split_radix::NTT<mod, root>;
    using ntt2 = split_radix::NTT<mod, root>;
    using ModInt = ntt::ModIntType;

    // std::cin >> s;
    size_t len = 1 << s;
    size_t times = 1000; // std::max<size_t>(1, (1 << 25) / len);
    std::vector<ModInt> a(len);
    std::vector<ModInt> b(len);
    std::vector<uint32_t> c(len);
    for (size_t i = 0; i < len; i++)
    {
        a[i] = uint64_t(i);
        b[i] = uint64_t(i);
        c[i] = uint64_t(i);
    }

    auto t1 = std::chrono::steady_clock::now();
    for (size_t i = 0; i < times; i++)
    {
        // ntt_dif<root>(a.data(), len);
        // ntt_dif<root>(a.data(), len);
        ntt_dit<root>(a.data(), len);
        // ntt2::dit2488(a.data(), len);
        // poly::fast_number_theoretic_transform_core::NTT(c.data(), len);
        // poly::fast_number_theoretic_transform_core::NTT(c.data(), len);
        // poly::fast_number_theoretic_transform_core::INTT<true>(c.data(), len);
    }
    auto t2 = std::chrono::steady_clock::now();
    for (size_t i = 0; i < times; i++)
    {
        // ntt::dif244(b.data(), len);
        // ntt::dif244(b.data(), len);
        ntt::dit244(b.data(), len);
        // ntt_dif<root>(a.data(), len);
        // ntt_dif<root>(a.data(), len);
        // ntt_dit<root>(b.data(), len);
    }
    auto t3 = std::chrono::steady_clock::now();
    auto time1 = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1).count();
    auto time2 = std::chrono::duration_cast<std::chrono::duration<double>>(t3 - t2).count();
    for (size_t i = 0; i < std::min<size_t>(len, 1024); i++)
    {
        if (uint64_t(a[i]) != uint64_t(b[i]))
        {
            std::cout << i << ":\t" << uint64_t(a[i]) << "\t" << uint64_t(b[i]) << "\n";
            return;
        }
    }
    std::cout << s << ":\n";
    std::cout << time1 << "\t" << time2 << "\t" << time1 / time2 << "\n";
}

template <typename T>
std::vector<T> poly_multiply(const std::vector<T> &in1, const std::vector<T> &in2)
{
    using namespace hint::transform::ntt;
    using NTT = NTT1;
    size_t len1 = in1.size(), len2 = in2.size(), out_len = len1 + len2;
    std::vector<T> result(out_len);
    size_t ntt_len = hint::int_floor2(out_len);
    std::vector<NTT::ModIntType> buffer1(ntt_len), buffer2(ntt_len);
    std::copy(in1.begin(), in1.end(), buffer1.begin());
    std::copy(in2.begin(), in2.end(), buffer2.begin());
    NTT::convolutionRecursion(buffer1.data(), buffer2.data(), buffer1.data(), ntt_len);
    for (size_t i = 0; i < out_len; i++)
    {
        result[i] = static_cast<T>(buffer1[i]);
    }
    return result;
}

template <typename T>
void result_test(const std::vector<T> &res, uint64_t ele)
{
    size_t len = res.size();
    for (size_t i = 0; i < len / 2; i++)
    {
        uint64_t x = (i + 1) * ele * ele;
        uint64_t y = res[i];
        if (x != y)
        {
            std::cout << "fail:" << i << "\t" << (i + 1) * ele * ele << "\t" << y << "\n";
            return;
        }
    }
    for (size_t i = len / 2; i < len; i++)
    {
        uint64_t x = (len - i - 1) * ele * ele;
        uint64_t y = res[i];
        if (x != y)
        {
            std::cout << "fail:" << i << "\t" << x << "\t" << y << "\n";
            return;
        }
    }
    std::cout << "success\n";
}
void test_ntt()
{
    int n = 18;
    std::cin >> n;
    size_t len = size_t(1) << n; // 变换长度
    uint64_t ele = 9;
    std::vector<uint64_t> in1(len / 2, ele);
    std::vector<uint64_t> in2(len / 2, ele); // 计算两个长度为len/2，每个元素为ele的卷积
    auto t1 = std::chrono::steady_clock::now();
    std::vector<uint64_t> res = poly_multiply(in1, in2);
    auto t2 = std::chrono::steady_clock::now();
    result_test<uint64_t>(res, ele); // 结果校验
    std::cout << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() << "us\n";
}

void ntt_perf(int s = 24)
{
    using namespace hint;
    using namespace transform;
    using namespace ntt;
    constexpr uint64_t mod1 = 1945555039024054273, root1 = 5;
    constexpr uint64_t mod2 = 4179340454199820289, root2 = 3;
    constexpr uint64_t mod = mod2, root = root2;

    using ntt = split_radix::NTT<mod, root, Uint128>;
    using ModInt = ntt::ModIntType;

    size_t len = 1 << s;
    size_t times = std::max<size_t>(1, (1 << 22) / len);
    std::vector<ModInt> a(len);

    auto t1 = std::chrono::steady_clock::now();
    for (size_t i = 0; i < len; i++)
    {
        a[i] = uint64_t(i);
    }
    for (size_t i = 0; i < times; i++)
    {
        ntt::dif244(a.data(), len);
        ntt::dif244(a.data(), len);
        ntt::dit244(a.data(), len);
    }
    auto t2 = std::chrono::steady_clock::now();
    for (size_t i = 0; i < len; i++)
    {
        a[i] = uint64_t(i);
    }
    for (size_t i = 0; i < times; i++)
    {
        ntt_dif<root>(a.data(), len);
        ntt_dif<root>(a.data(), len);
        ntt_dit<root>(a.data(), len);
    }
    auto t3 = std::chrono::steady_clock::now();
    auto time1 = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1).count();
    auto time2 = std::chrono::duration_cast<std::chrono::duration<double>>(t3 - t2).count();
    std::cout << s << ":\n";
    std::cout << time1 << "\t" << time2 << "\t" << time2 / time1 << "\n";
}

void ntt_perf_loop()
{
    for (int i = 10; i <= 30; i++)
    {
        ntt_perf(i);
    }
}

void test_add()
{
    std::cout << "test_add\n";
    using namespace hint;
    size_t len = 1e8;
    constexpr uint64_t BASE = 1e18;
    std::vector<uint64_t> a(len, BASE - 1);
    std::vector<uint64_t> b(len, BASE - 1);
    std::vector<uint64_t> c(len + 1);
    constexpr hint::arithm::support::BaseExecutor<uint64_t> executor(BASE);
    auto t1 = std::chrono::steady_clock::now();
    hint::arithm::addition_binary::abs_add_long_long(a.data(), a.size(), b.data(), b.size(), c.data(), executor);
    auto t2 = std::chrono::steady_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() << "us\n";
}

int main()
{
    // test_add_VVW();
    // test_ntt();
    // test_mul64();
    test_add();
    // test_div();
}