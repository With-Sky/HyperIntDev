#include "hint_decimal.hpp"

int main()
{
    std::string s = "12345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890";
    hint::HyperDecimal h = s;
    std::cout << std::boolalpha << (h.toString() == s) << std::endl;
    return 0;
}