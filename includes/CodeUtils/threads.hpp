#ifndef INCLUDES_CODEUTILS_THREADS_HPP
#define INCLUDES_CODEUTILS_THREADS_HPP

namespace codeUtils
{
    bool isOpenMpDefined();
    void setOpenMpNumberOfThreads(int num);
} // namespace codeUtils
#endif
