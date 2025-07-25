#ifndef INCLUDE_UTIL_CASTING_HPP
#define INCLUDE_UTIL_CASTING_HPP

#include <stdexcept>
#include <string>

namespace gmml
{
    namespace util
    {
        template<class T, class U> T erratic_cast(U obj, const std::string& msg = "")
        {
            T result = dynamic_cast<T>(obj);
            if (result == nullptr)
            {
                throw std::runtime_error("Dynamic casting returned nullptr: " + msg + "\n");
            }
            return result;
        }
    } // namespace util
} // namespace gmml
#endif
