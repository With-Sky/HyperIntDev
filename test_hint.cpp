#include <chrono>
#include <random>
#include "hint.hpp"

#include <algorithm>

// out = in * num_mul + num_add, return carry
template <typename NumTy, typename Executor>
inline NumTy abs_mul_add_num_half(const NumTy in[], size_t len, NumTy out[], NumTy num_add, NumTy num_mul, const Executor &exec)
{
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
inline NumTy abs_mul_add_num_add_self_half(const NumTy in[], size_t len, NumTy out[], NumTy num_add, NumTy num_mul, const Executor &exec)
{
    size_t i = 0;
    for (const size_t rem_len = len - len % 4; i < rem_len; i += 4)
    {
        // num_add = exec.mulAddPlusSelfX4(in + i, num_mul, num_add, out + i);
    }
    for (; i < len; i++)
    {
        exec.mulAddPlusSelf(in[i], num_mul, out[i], num_add);
    }
    return num_add;
}

void test_abs_mul_add_num_half()
{
    uint64_t base = 0, num = base - 1;
    hint::utility::BaseExecutorBinary<uint64_t> exec;
    constexpr size_t len1 = 1e8, len2 = len1, loop = 1e6;
    std::vector<uint64_t> a(len1, base - 1);
    std::vector<uint64_t> b(len2, base - 1);
    auto t1 = std::chrono::high_resolution_clock::now();
    auto c = abs_mul_add_num_add_self_half(a.data(), len1, b.data(), num, num, exec);
    auto t2 = std::chrono::high_resolution_clock::now();
    for (auto &i : b)
    {
        // std::cout << i << " ";
    }
    std::cout << std::endl;
    std::cout << c << std::endl;
    std::cout << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() << "us" << std::endl;
}

inline uint64_t rand_num(uint64_t mod)
{
    constexpr uint64_t prime = 1000000000000002493ull;
    static uint64_t seed = 0;
    seed = seed * prime + 1;
    return seed % mod;
}

void test_div()
{
    using ui128 = __uint128_t;
    uint64_t base = 1000'000'000'000'000'0000ull;
    // uint64_t base = INT64_MAX + 1;
    std::cin >> base;
    hint::utility::DivExecutor<uint64_t> exec(base);
    size_t len = 1e8;
    std::vector<uint64_t> a(len), b(len);
    for (auto &i : a)
    {
        i = rand_num(base);
        // i = base - 1;
    }
    uint64_t sum_quot1 = 0, sum_rem1 = 0;
    uint64_t sum_quot2 = 0, sum_rem2 = 0;
    auto t1 = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < len - 1; i++)
    {
        uint64_t quot1, rem1;
        ui128 prod = (ui128(a[i]) << 64) | a[i + 1];
        // quot1 = exec.divRemNorm(prod, rem1);
        quot1 = exec.divRem(prod >> 64, prod, rem1);
        sum_quot1 += quot1;
        sum_rem1 += rem1;
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < len - 1; i++)
    {
        uint64_t quot2, rem2;
        ui128 prod = (ui128(a[i]) << 64) | a[i + 1];
        // quot2 = prod / base;
        // rem2 = prod % base;
        sum_quot2 += quot2;
        sum_rem2 += rem2;
    }
    auto t3 = std::chrono::high_resolution_clock::now();
    std::cout << sum_quot1 << " " << sum_rem1 << std::endl;
    std::cout << sum_quot2 << " " << sum_rem2 << std::endl;
    std::cout << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() << "us" << std::endl;
    std::cout << std::chrono::duration_cast<std::chrono::microseconds>(t3 - t2).count() << "us" << std::endl;
}

template <uint32_t BASE>
void mul_basic(uint32_t in1[], size_t len1, uint32_t in2[], size_t len2, uint32_t out[])
{
    for (size_t i = 0; i < len1; i++)
    {
        uint64_t carry = 0;
        for (size_t j = 0; j < len2; j++)
        {
            uint64_t prod = (uint64_t)in1[i] * in2[j];
            prod += carry + out[i + j];
            out[i + j] = prod % BASE;
            carry = prod / BASE;
        }
        out[i + len2] = carry;
    }
}

void test_mul_basic10()
{
    constexpr uint64_t base = 1e9, num = base - 1;
    // utility::BaseExecutorBinary exec;
    size_t len1 = 2e4, len2 = len1, loop = 1e5;
    // std::cin >> len1 >> len2 >> loop;
    std::vector<uint32_t> a(len1, num);
    std::vector<uint32_t> b(len2, num);
    std::vector<uint32_t> c(len1 + len2);
    static std::mt19937 gen(1);
    std::uniform_int_distribution<uint64_t> dist(0, num);
    // for (auto &&i : a)
    // {
    //     i = dist(gen);
    // }
    // for (auto &&i : b)
    // {
    //     i = dist(gen);
    // }
    auto t1 = std::chrono::steady_clock::now();
    // for (size_t i = 0; i < loop; i++)
    {
        mul_basic<base>(a.data(), len1, b.data(), len2, c.data());
    }
    auto t2 = std::chrono::steady_clock::now();
    for (auto it = c.rbegin(); it != c.rend(); it++)
    {
        std::cout << *it << " ";
    }
    std::cout << std::endl;
    std::cout << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() << "us" << std::endl;
}

void test_mul_basic()
{
    using namespace hint;
    using namespace arithm::multiplication;
    using namespace arithm::addition;
    uint64_t base = 1e19, num = base - 1;
    // BaseExecutor<uint64_t> exec(base);
    utility::BaseExecutorBinary<uint64_t> exec;
    size_t len1 = (1e5), len2 = len1, loop = 1 << 20;
    // std::cin >> len1 >> len2 >> loop;
    std::vector<uint64_t> a(len1, num);
    std::vector<uint64_t> b(len2, num);
    // auto a_p = new uint64_t[len1];
    // auto b_p = new uint64_t[len2];
    // std::cout << a_p << " " << b_p << std::endl;
    std::vector<uint64_t> c(len1 + len2);
    std::vector<uint64_t> d(len1 + len2);
    static std::mt19937 gen(1);
    // b = a = {12271144805854719029, 14901432754109524217, 3976023265224744767, 3691743531038793554, 10518674766591855369, 16040479651847935293, 6547057120208965452, 13555195664566830802, 6588984070404368609, 95284951371483940, 14616649413911854801, 908727287399969004, 1122504066081417718, 12597235699414803597, 6487115003628942535, 9016559605287246797, 6501391612429803521, 9335840991529667284, 7536087634734216999, 1096566204796952383, 2521600334672742513, 5327178105255534790, 6563796726199367596, 4264676419002213608, 16152462487545995875, 18022713201819445811, 14989106788728424597, 11889449673523465139, 816205900141760993, 7603667619536851764, 9544754529937541261, 3405619203695849313, 9253668332735781264, 3772883358476565758, 5767975708667686435, 17411129618902945802, 2747094638218966479, 6555570056851783241, 4429752750818721395, 18098432011149379461, 7092532507664853421, 10226239231214421834, 9673825486394464490, 7756997698450048382, 4619110562266398829, 14709072737518605339, 1006481701186435726, 13959947564169318589, 10318185079568189768, 2342163127147520799, 8381760586784675853, 12450498320086370081, 215252575380578374, 5097522528455575746, 4466973776777661050, 12580805000889164358, 5166050641316974184, 4745928933147667131, 9507146309724095926, 10625921511001148772, 3473640986122561753, 7459972929883376396, 10696605722387200250, 5007028558846880153, 624310316918951752, 15344318165583501092, 17221732859175648399, 18277164935840980329, 13688036639861191966, 2422004203947623829, 3337298988363341972, 10237377498548984045, 9528633094567476014, 14751330182012594658, 14144291570796327767, 8635407794992039412, 13929108628708672478, 14793142760908167813, 1258472876168276427, 9347237275004853734};
    std::uniform_int_distribution<uint64_t> dist(0, num);
    for (auto &&i : a)
    {
        i = dist(gen);
    }
    for (auto &&i : b)
    {
        i = dist(gen);
    }
    // a = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 18446744073709551615};
    // b = {18446744073709551615, 18446744073709551615, 18446744073709551615, 18446744073709551615, 18446744073709551615, 18446744073709551615, 18446744073709551615, 18446744073709551615, 18446744073709551615, 18446744073709551615, 18446744073709551615, 18446744073709551615, 18446744073709551615, 18446744073709551615, 18446744073709551615, 18446744073709551615, 18446744073709551615, 18446744073709551615, 18446744073709551615, 18446744073709551615, 18446744073709551615, 18446744073709551615, 18446744073709551615, 18446744073709551615, 18446744073709551615, 18446744073709551615, 18446744073709551615, 18446744073709551615, 18446744073709551615, 18446744073709551615, 18446744073709551615, 18446744073709551615, 18446744073709551615};
    uint64_t carry = 0;
    auto t1 = std::chrono::steady_clock::now();
    // for (size_t i = 0; i < loop; i++)
    {
        abs_mul_karatusba(a.data(), a.size(), b.data(), b.size(), d.data(), exec);
        // abs_sub(a.data(), a.size(), b.data(), b.size(), a.data(), exec);
        // carry = abs_add_equal_mul_add_num(a.data(), len1, b.data(), uint64_t(0), num, exec);
        // abs_mul_basic(a.data(), a.size(), b.data(), b.size(), c.data(), exec);
    }
    auto t2 = std::chrono::steady_clock::now();
    // for (size_t i = 0; i < loop; i++)
    {
        // abs_mul_basic(a.data(), a.size(), b.data(), b.size(), c.data(), exec);
        abs_mul_ntt(a.data(), a.size(), b.data(), b.size(), c.data(), exec);
        // abs_mul_karatusba1(a.data(), a.size(), b.data(), b.size(), c.data(), exec);
        // abs_mul_basic_bin(a.data(), 40, a.data() + 40, 40, a.data());
        // abs_mul_karatusba(a.data(), a.size(), b.data(), b.size(), d.data(), exec);
        // abs_mul_basic(a.data(), a.size(), b.data(), b.size(), d.data(), exec, false);
        // abs_mul_basic(a.data(), a.size(), b.data(), b.size(), d.data(), exec);
        // abs_mul_add_num(a.data(), a.size(), d.data(), 0ull, base - 2, exec);
    }
    auto t3 = std::chrono::steady_clock::now();

    for (auto it = c.rbegin(); it != c.rend(); it++)
    {
        // std::cout << *it << " ";
    }
    // std::cout << carry << std::endl;
    std::cout << "\n";
    for (auto it = d.rbegin(); it != d.rend(); it++)
    {
        // std::cout << *it << " ";
    }
    // std::cout << carry << std::endl;
    // std::cout << "Size" << b.size() << std::endl;
    for (auto it = b.rbegin(); it != b.rend(); it++)
    {
        // std::cout << *it << " ";
    }
    // for(size_t i = 0; i < len2; i++)
    // {
    //     std::cout << b_p[i] << " ";
    // }
    for (auto i : a)
    {
        // std::cout << i << " ";
    }
    std::cout << "\n";

    std::cout << std::boolalpha << (c == d) << "\n";
    // std::cout << std::boolalpha << (a == b) << "\n";

    std::cout << "\n";
    std::cout << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() << "us" << std::endl;
    std::cout << std::chrono::duration_cast<std::chrono::microseconds>(t3 - t2).count() << "us" << std::endl;

    // std::cout<< a_p << " " << b_p << std::endl;
    // delete[] a_p;
    // delete[] b_p;
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

template <typename T>
std::vector<T> poly_multiply(const std::vector<T> &in1, const std::vector<T> &in2)
{
    using namespace hint::transform::ntt;
    using NTT = hint::transform::ntt::NTT64_1;
    // using NTT = hint::transform::ntt::NTT<998244353,3>;
    size_t len1 = in1.size(), len2 = in2.size(), out_len = len1 + len2;
    std::vector<T> result(out_len);
    size_t ntt_len = NTT::findFitLen(out_len);
    std::vector<NTT::IntType> buffer1(ntt_len), buffer2(ntt_len);
    std::copy(in1.begin(), in1.end(), buffer1.begin());
    std::copy(in2.begin(), in2.end(), buffer2.begin());
    NTT::convolution(buffer1.data(), buffer2.data(), ntt_len);
    std::copy(buffer1.begin(), buffer1.begin() + out_len, result.begin());
    return result;
}
void test_ntt()
{
    int n = 18;
    std::cin >> n;
    size_t len = size_t(1) << n; // 变换长度
    // len = 117440512;
    uint64_t ele = 2;
    std::vector<uint32_t> in1(len / 2, ele);
    std::vector<uint32_t> in2(len / 2, ele); // 计算两个长度为len/2，每个元素为ele的卷积
    auto t1 = std::chrono::steady_clock::now();
    std::vector<uint32_t> res = poly_multiply(in1, in2);
    auto t2 = std::chrono::steady_clock::now();
    result_test(res, ele); // 结果校验
    std::cout << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() << "us\n";
}

#include "../../test/bind_cpu.hpp"

void test_hint_str()
{
    hint::HyperIntHex a, b;
    a.fromString(std::string(1000000, '9'));
    b.fromString(std::string(1000000, '7'));
    std::cout << a.limbSize() << std::endl;
    auto t1 = std::chrono::steady_clock::now();
    a *= b;
    auto t2 = std::chrono::steady_clock::now();
    // auto s = a.toString();
    // std::cout << s << std::endl;
    std::cout << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() << "us\n";
}

void test_branch(uint64_t p[], size_t len)
{
    using NTT = hint::transform::ntt::NTT64_1;
    using Mint = NTT::ModInt;
    auto t1 = std::chrono::steady_clock::now();
    size_t stride = len / 4;
    auto it0 = reinterpret_cast<Mint *>(p);
    auto it1 = it0 + stride, it2 = it1 + stride, it3 = it2 + stride, end = it1;
    Mint omega = stride;
    for (; it3 < end; it0++, it1++, it2++, it3++)
    {
        auto t0 = it0[0].norm2(), t1 = it1[0].norm2();
        auto t2 = it2[0].norm2(), t3 = it3[0].norm2();
        // hint::transform::transform2(t0, t2);
        // auto diff = t1.sub(t3);
        // t1 = t1 + t3;
        // t3 = diff * omega;
        it0[0] = t0.add(t1), it1[0] = t0.sub(t1), it2[0] = t2.add(t3), it3[0] = t2.sub(t3);
    }
    auto t2 = std::chrono::steady_clock::now();
}

int main()
{
    bind_cpu(0);
    // test_hint_str();
    test_mul_basic();
    // test_ntt();

    // test_div();
    // test_mul_basic();
    // test_mul_basic10();
    // test_abs_mul_add_num_half();
    std::cin.get();
}