#ifndef INCLUDE_UTIL_FORMATTING_HPP
#define INCLUDE_UTIL_FORMATTING_HPP

#include <ostream>

namespace gmml
{
    namespace util
    {
        enum class textAlignment
        {
            left,
            right
        };

        struct decimalFormat
        {
            uint precision;
            uint digits;
        };

        struct floatFormat
        {
            textAlignment alignment;
            uint width;
            decimalFormat decimals;
        };

        void writeFloat(std::ostream& stream, const floatFormat& format, double value);

    } // namespace util
} // namespace gmml
#endif
