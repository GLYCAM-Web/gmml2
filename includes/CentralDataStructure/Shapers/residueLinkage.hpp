#ifndef INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_RESIDUELINKAGE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_RESIDUELINKAGE_HPP
/*
 * This class figures out the rotatable bonds between two residues
 * Starts/ends at the CA atoms in proteins. Looks for cycles (as they aren't rotatable).
 * Stores each rotatable bond as a RotatableDihedral object.
 */
#include "includes/CentralDataStructure/cdsTypes.hpp"
#include "includes/CentralDataStructure/Geometry/orientation.hpp"
#include "includes/MolecularMetadata/GLYCAM/dihedralangledata.hpp"
#include "includes/CodeUtils/logging.hpp"

#include <string>
#include <array>
#include <vector>
#include <utility>

namespace cds
{
    struct ResidueLink
    {
        std::pair<cds::Residue*, cds::Residue*> residues;
        std::pair<cds::Atom*, cds::Atom*> atoms;
    };

    struct ResidueLinkAttributes
    {
        std::pair<ResidueAttributes, ResidueAttributes> residues;
        std::pair<std::string, std::string> atoms;
    };

    typedef std::array<cds::Atom*, 4> DihedralAtoms;

    struct RotatableDihedral
    {
        // The four atoms that define the dihedral angle. The bond between atom2_ and atom3_ is what is rotated.
        std::array<Atom*, 4> atoms;
        // A vector of pointers to the atoms that are connected to atom2_ and atom3_, and will be rotated when that bond
        // is rotated.
        std::vector<Atom*> movingAtoms;
        size_t currentMetadataIndex;
    };

    struct ResidueLinkage
    {
        ResidueLink link;
        std::vector<RotatableDihedral> rotatableDihedrals;
        std::vector<std::vector<size_t>> dihedralMetadata;
        GlycamMetadata::RotamerType rotamerType;
        bool isDerivative;
        unsigned long long index = 0;
        std::string name         = ""; // e.g. "DGalpb1-6DGlcpNAc". It being empty works with GetName();
    };

    std::vector<size_t> rotatableDihedralsWithMultipleRotamers(const std::vector<std::vector<size_t>>& metadata);
    std::vector<ResidueLinkage> nonDerivativeResidueLinkages(const std::vector<cds::ResidueLinkage>& linkages);
    size_t numberOfShapes(const GlycamMetadata::DihedralAngleDataTable& metadataTable,
                          GlycamMetadata::RotamerType rotamerType, const std::vector<std::vector<size_t>>& metadata);
    size_t numberOfLikelyShapes(const GlycamMetadata::DihedralAngleDataTable& metadataTable,
                                GlycamMetadata::RotamerType rotamerType,
                                const std::vector<std::vector<size_t>>& metadata);
    DihedralCoordinates dihedralCoordinates(const cds::RotatableDihedral& dihedral);
    std::string print(const ResidueLink& link);
    std::string print(const RotatableDihedral& dihedral);
    std::string print(const GlycamMetadata::DihedralAngleDataTable& table, const ResidueLinkage& linkage);

} // namespace cds
#endif
