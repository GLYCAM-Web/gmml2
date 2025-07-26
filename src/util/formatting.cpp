#include "include/util/formatting.hpp"

#include <iomanip>
#include <ostream>

namespace gmml
{
    namespace util
    {
        void writeFloat(std::ostream& stream, const floatFormat& format, double value)
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
