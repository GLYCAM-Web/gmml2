#ifndef INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_RESIDUELINKAGE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_RESIDUELINKAGE_HPP
/*
 * This class figures out the rotatable bonds between two residues
 * Starts/ends at the CA atoms in proteins. Looks for cycles (as they aren't rotatable).
 * Stores each rotatable bond as a RotatableDihedral object.
 */
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/Geometry/orientation.hpp"
#include "includes/MolecularMetadata/GLYCAM/dihedralangledata.hpp"
#include "includes/CodeUtils/logging.hpp"

#include <string>
#include <array>
#include <vector>
#include <utility>

namespace cds
{
    using GlycamMetadata::DihedralAngleData;
    using GlycamMetadata::DihedralAngleDataVector;

    struct ResidueLink
    {
        std::pair<cds::Residue*, cds::Residue*> residues;
        std::pair<cds::Atom*, cds::Atom*> atoms;
    };

    struct ResidueLinkNames
    {
        std::pair<std::string, std::string> residues;
        std::pair<std::string, std::string> atoms;
    };

    struct DihedralAtoms
    {
        bool isBranching;
        std::array<cds::Atom*, 4> atoms;
    };

    struct RotatableDihedral
    {
        bool isBranchingLinkage;
        // The four atoms that define the dihedral angle. The bond between atom2_ and atom3_ is what is rotated.
        std::array<Atom*, 4> atoms;
        // A vector of pointers to the atoms that are connected to atom2_ and atom3_, and will be rotated when that bond
        // is rotated.
        std::vector<Atom*> movingAtoms;
        size_t currentMetadataIndex;
    };

    typedef std::vector<GlycamMetadata::DihedralAngleDataVector> DihedralAngleMetadata;

    struct ResidueLinkage
    {
        ResidueLink link;
        std::vector<RotatableDihedral> rotatableDihedrals;
        std::vector<DihedralAngleDataVector> dihedralMetadata;
        GlycamMetadata::RotamerType rotamerType;
        unsigned long long index = 0;
        std::string name         = ""; // e.g. "DGalpb1-6DGlcpNAc". It being empty works with GetName();
        std::vector<cds::Residue*> reducingOverlapResidues;
        std::vector<cds::Residue*> nonReducingOverlapResidues;
    };

    std::vector<size_t> rotatableDihedralsWithMultipleRotamers(const std::vector<DihedralAngleDataVector>& metadata);
    size_t numberOfShapes(GlycamMetadata::RotamerType rotamerType,
                          const std::vector<DihedralAngleDataVector>& metadata);
    size_t numberOfLikelyShapes(GlycamMetadata::RotamerType rotamerType,
                                const std::vector<DihedralAngleDataVector>& metadata);
    DihedralCoordinates dihedralCoordinates(const cds::RotatableDihedral& dihedral);
    std::string print(const ResidueLink& link);
    std::string print(const RotatableDihedral& dihedral);
    std::string print(const ResidueLinkage& linkage);

} // namespace cds
#endif
