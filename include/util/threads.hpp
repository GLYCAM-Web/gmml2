#ifndef INCLUDE_UTIL_THREADS_HPP
#define INCLUDE_UTIL_THREADS_HPP

namespace gmml
{
    namespace util
    {
        bool isOpenMpDefined();
        void setOpenMpNumberOfThreads(int num);
    } // namespace util
} // namespace gmml
#endif
