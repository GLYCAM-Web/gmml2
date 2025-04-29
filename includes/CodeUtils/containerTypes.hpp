#ifndef INCLUDES_CODEUTILS_CONTAINERTYPES_HPP
#define INCLUDES_CODEUTILS_CONTAINERTYPES_HPP

#include <vector>

namespace codeUtils
{
    template<typename T> struct SparseVector
    {
        std::vector<T> values;
        std::vector<bool> hasValue;
    };
} // namespace codeUtils
#endif
