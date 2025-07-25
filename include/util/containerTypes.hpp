#ifndef INCLUDE_UTIL_CONTAINERTYPES_HPP
#define INCLUDE_UTIL_CONTAINERTYPES_HPP

#include <vector>

namespace gmml
{
    namespace util
    {
        template<typename T> struct SparseVector
        {
            std::vector<T> values;
            std::vector<bool> hasValue;
        };
    } // namespace util
} // namespace gmml
#endif
