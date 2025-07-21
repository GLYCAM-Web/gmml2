#include "include/util/threads.hpp"

namespace gmml
{
    namespace util
    {
#ifdef _OPENMP
#include <omp.h>

        bool isOpenMpDefined() { return true; }

        void setOpenMpNumberOfThreads(int num) { omp_set_num_threads(num); }

#else

        bool isOpenMpDefined() { return false; }

        void setOpenMpNumberOfThreads(int) {}

#endif /* _OPENMP */
    }  // namespace util
} // namespace gmml
