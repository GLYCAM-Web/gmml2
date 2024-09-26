#ifndef INCLUDES_CODEUTILS_CASTING_HPP
#define INCLUDES_CODEUTILS_CASTING_HPP

#include <stdexcept>
#include <string>

namespace codeUtils
{
    template<class T, class U> T erratic_cast(U obj, const std::string& msg = "")
    {
        T result = dynamic_cast<T>(obj);
        if (result == nullptr)
        {
            throw std::runtime_error("casting returned nullptr: " + msg);
        }
        return result;
    }
} // namespace codeUtils
#endif
