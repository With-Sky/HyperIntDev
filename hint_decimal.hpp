#ifndef HINT_DECIMAL_HPP
#define HINT_DECIMAL_HPP
#include "hint_arithm.hpp"

namespace hint
{
    template <typename WordTy>
    constexpr int getBaseDigits()
    {
        switch (sizeof(WordTy))
        {
        case 1:
            return 2;
            break;
        case 2:
            return 4;
            break;
        case 4:
            return 9;
            break;
        case 8:
            return 18;
            break;
        default:
            return 18;
            break;
        }
        return 0;
    }
    class HyperDecimal
    {
    public:
        HyperDecimal() = default;
        HyperDecimal(const HyperDecimal &other) = default;
        HyperDecimal(HyperDecimal &&other) = default;
        HyperDecimal &operator=(const HyperDecimal &other) = default;
        HyperDecimal &operator=(HyperDecimal &&other) = default;
        ~HyperDecimal() = default;

        HyperDecimal(const std::string &str)
        {
            size_t len = checkStringIn(str.data(), str.size());
            data.resize(strLenToWordLen(len));
            strToVec(str.data(), len, data.data());
        }
        HyperDecimal(const char *str)
        {
            size_t len = checkStringIn(str);
            data.resize(strLenToWordLen(len));
            strToVec(str, len, data.data());
        }
        HyperDecimal(const std::vector<uint64_t> &data_in)
        {
            data = data_in;
        }

        std::string toString() const
        {
            size_t len = wordLen();
            if (0 == len)
            {
                return "0";
            }
            std::string str = std::to_string(data[len - 1]);
            len--;
            while (len > 0)
            {
                len--;
                str += hint::ui64to_string_base10(data[len], BASE_DIGITS); // TODO: optimize
            }
            return str;
        }

    private:
        using Word = uint64_t;
        using DataVector = std::vector<uint64_t>;
        static constexpr size_t strLenToWordLen(size_t str_len)
        {
            return (str_len + BASE_DIGITS - 1) / BASE_DIGITS;
        }
        static constexpr size_t checkStringIn(const char str[], size_t len)
        {
            for (size_t i = 0; i < len; i++)
            {
                if (str[i] < '0' || str[i] > '9')
                {
                    return i;
                }
            }
            return len;
        }
        static constexpr size_t checkStringIn(const char str[])
        {
            if (nullptr == str)
            {
                return 0;
            }
            const char *p = str;
            while (*p)
            {
                if (*p < '0' || *p > '9')
                {
                    return p - str;
                }
                p++;
            }
            return p - str;
        }

        static constexpr Word strToWord(const char *str, const char *end)
        {
            Word res = 0;
            while (str < end)
            {
                res = res * 10 + (*str - '0');
                str++;
            }
            return res;
        }

        static void strToVec(const char str[], size_t len, Word out[])
        {
            if (0 == len)
            {
                return;
            }
            auto end = str + len;
            size_t i = 0;
            while (end - BASE_DIGITS > str)
            {
                out[i] = strToWord(end - BASE_DIGITS, end);
                i++;
                end -= BASE_DIGITS;
            }
            out[i] = strToWord(str, end);
        }

        static constexpr int BASE_DIGITS = getBaseDigits<Word>();
        static constexpr Word BASE = hint::qpow<Word>(10, BASE_DIGITS);

        size_t wordLen() const
        {
            return data.size();
        }

        DataVector data;
    };

    constexpr int HyperDecimal::BASE_DIGITS;
    constexpr HyperDecimal::Word HyperDecimal::BASE;
}
#endif