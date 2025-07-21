#ifndef INCLUDES_CODEUTILS_FORMATTING_HPP
#define INCLUDES_CODEUTILS_FORMATTING_HPP

#include <iomanip>
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

        inline void writeFloat(std::ostream& stream, const floatFormat& format, double value)
        {
            const decimalFormat& decimals = format.decimals;
            stream << (format.alignment == textAlignment::left ? std::left : std::right);
            stream << std::setw(format.width) << std::fixed;
            stream << std::setprecision(decimals.precision) << value;
            for (size_t n = 0; n < decimals.digits - decimals.precision; n++)
            {
                stream << "0";
            }
        }

    } // namespace util
} // namespace gmml
#endif
