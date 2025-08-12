#ifndef INCLUDE_CARBOHYDRATE_PDBWRITER_HPP
#define INCLUDE_CARBOHYDRATE_PDBWRITER_HPP

#include "include/assembly/assemblyTypes.hpp"
#include "include/carbohydrate/carbohydrateTypes.hpp"

#include <ostream>
#include <string>
#include <vector>

namespace gmml
{
    void writePdb(
        std::ostream& stream,
        const carbohydrate::CarbohydrateData& data,
        const assembly::Graph& graph,
        const std::vector<std::string>& headerLines);
} // namespace gmml
#endif
