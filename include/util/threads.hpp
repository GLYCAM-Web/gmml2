#ifndef INCLUDES_CODEUTILS_THREADS_HPP
#define INCLUDES_CODEUTILS_THREADS_HPP

namespace gmml
{
    namespace util
    {
        bool isOpenMpDefined();
        void setOpenMpNumberOfThreads(int num);
    } // namespace util
} // namespace gmml
#endif
