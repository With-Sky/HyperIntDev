#include "hint_decimal.hpp"
#include <chrono>

void testbase()
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

void testbin()
{
    using namespace hint::arithm::multiplication;

    size_t len1 = 1 << 15, len2 = 1 << 15;
    std::vector<uint64_t> a(len1, UINT64_MAX), b(len2, UINT64_MAX), c(len1 + len2, 0), buffer(len1 + len2, 0);
    auto t1 = std::chrono::high_resolution_clock::now();
    abs_mul64_classic(a.data(), len1, b.data(), len2, c.data(), buffer.data(), buffer.data() + len1 + len2);
    auto t2 = std::chrono::high_resolution_clock::now();

    // hint::arithm::support::ary_print(c.data(), len1 + len2);
    std::cout << "Time: " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() << "us" << std::endl;
}

void testbin2()
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
}

int main()
{
    std::ios::sync_with_stdio(false);
    std::cin.tie(0);
    testbase();
    testbin2();
    testbin();
    return 0;
}