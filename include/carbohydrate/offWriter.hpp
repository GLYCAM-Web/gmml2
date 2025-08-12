#ifndef INCLUDE_CARBOHYDRATE_OFFWRITER_HPP
#define INCLUDE_CARBOHYDRATE_OFFWRITER_HPP

#include "include/carbohydrate/carbohydrateTypes.hpp"
#include "include/fileType/off/offFileData.hpp"

#include <ostream>
#include <string>
#include <vector>

namespace gmml
{
    std::vector<std::vector<size_t>> atomsConnectedToOtherResidues(const assembly::Graph& graph);
    off::OffFileData toOffFileData(const carbohydrate::CarbohydrateData& data, const assembly::Graph& graph);
} // namespace gmml
#endif
