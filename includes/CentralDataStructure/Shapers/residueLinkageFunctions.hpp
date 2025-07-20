#ifndef INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_RESIDUELINKAGEFUNCTIONS_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_RESIDUELINKAGEFUNCTIONS_HPP

#include "includes/CentralDataStructure/Geometry/orientation.hpp"
#include "includes/CentralDataStructure/Shapers/residueLinkageTypes.hpp"
#include "includes/CentralDataStructure/cdsTypes.hpp"
#include "includes/MolecularMetadata/GLYCAM/dihedralangledata.hpp"

#include <utility>
#include <vector>

namespace cds
{

    std::vector<size_t> rotatableDihedralsWithMultipleRotamers(const std::vector<std::vector<size_t>>& metadata);
    std::vector<ResidueLinkage> nonDerivativeResidueLinkages(const std::vector<ResidueLinkage>& linkages);

    size_t numberOfShapes(
        const GlycamMetadata::DihedralAngleDataTable& metadataTable,
        GlycamMetadata::RotamerType rotamerType,
        const std::vector<std::vector<size_t>>& metadata);

    size_t numberOfLikelyShapes(
        const GlycamMetadata::DihedralAngleDataTable& metadataTable,
        GlycamMetadata::RotamerType rotamerType,
        const std::vector<std::vector<size_t>>& metadata);

    DihedralCoordinates dihedralCoordinates(const RotatableDihedral& dihedral);
    std::string print(const ResidueLink& link);
    std::string print(const RotatableDihedral& dihedral);
    std::string print(const GlycamMetadata::DihedralAngleDataTable& table, const ResidueLinkage& linkage);
} // namespace cds
#endif
