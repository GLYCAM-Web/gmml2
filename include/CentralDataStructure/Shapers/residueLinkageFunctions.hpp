#ifndef INCLUDE_CENTRALDATASTRUCTURE_SHAPERS_RESIDUELINKAGEFUNCTIONS_HPP
#define INCLUDE_CENTRALDATASTRUCTURE_SHAPERS_RESIDUELINKAGEFUNCTIONS_HPP

#include "include/CentralDataStructure/Shapers/residueLinkageTypes.hpp"
#include "include/CentralDataStructure/cdsTypes.hpp"
#include "include/geometry/orientation.hpp"
#include "include/metadata/dihedralangledata.hpp"

#include <utility>
#include <vector>

namespace gmml
{
    std::vector<size_t> rotatableDihedralsWithMultipleRotamers(const std::vector<std::vector<size_t>>& metadata);
    std::vector<ResidueLinkage> nonDerivativeResidueLinkages(const std::vector<ResidueLinkage>& linkages);

    size_t numberOfShapes(
        const DihedralAngleDataTable& metadataTable,
        RotamerType rotamerType,
        const std::vector<std::vector<size_t>>& metadata);

    size_t numberOfLikelyShapes(
        const DihedralAngleDataTable& metadataTable,
        RotamerType rotamerType,
        const std::vector<std::vector<size_t>>& metadata);

    DihedralCoordinates dihedralCoordinates(const RotatableDihedral& dihedral);
    std::string print(const ResidueLink& link);
    std::string print(const RotatableDihedral& dihedral);
    std::string print(const DihedralAngleDataTable& table, const ResidueLinkage& linkage);
} // namespace gmml
#endif
