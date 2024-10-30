#include "hint_decimal.hpp"
#include <chrono>
int main()
{
    std::string s(1800000000, '9');
    hint::HyperDecimal h = s;
    auto t1 = std::chrono::high_resolution_clock::now();
    h = h + h;
    auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << h.wordLen() << std::endl;
    std::cout << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() << "us" << std::endl;
    return 0;
}