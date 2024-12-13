#include "includes/CodeUtils/threads.hpp"

#ifdef _OPENMP
#include <omp.h>

bool codeUtils::isOpenMpDefined()
{
    return true;
}

void codeUtils::setOpenMpNumberOfThreads(int num)
{
    omp_set_num_threads(num);
}

#else

bool codeUtils::isOpenMpDefined()
{
    return false;
}

void codeUtils::setOpenMpNumberOfThreads(int)
{}

#endif /* _OPENMP */
