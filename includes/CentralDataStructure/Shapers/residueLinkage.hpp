#ifndef INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_RESIDUELINKAGE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_RESIDUELINKAGE_HPP
/*
 * This class figures out the rotatable bonds between two residues
 * Starts/ends at the CA atoms in proteins. Looks for cycles (as they aren't rotatable).
 * Stores each rotatable bond as a RotatableDihedral object.
 */
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralAngles.hpp"
#include "includes/CentralDataStructure/Geometry/orientation.hpp"
#include "includes/MolecularMetadata/GLYCAM/dihedralangledata.hpp"
#include "includes/CodeUtils/logging.hpp"

#include <string>
#include <array>
#include <vector>
#include <utility>

namespace cds
{
    using gmml::MolecularMetadata::GLYCAM::DihedralAngleData;
    using gmml::MolecularMetadata::GLYCAM::DihedralAngleDataVector;

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
        RotatableDihedral(bool isBranchingLinkage, const std::array<Atom*, 4>& atoms_,
                          const std::vector<DihedralAngleData>& metadata)
            : isBranchingLinkage(isBranchingLinkage), atoms(atoms_), metadataVector(metadata),
              currentMetadataIndex(0) {};

        bool isBranchingLinkage;
        // The four atoms that define the dihedral angle. The bond between atom2_ and atom3_ is what is rotated.
        std::array<cds::Atom*, 4> atoms;
        // A vector of pointers to the atoms that are connected to atom2_ and atom3_, and will be rotated when that bond
        // is rotated.
        std::vector<cds::Coordinate*> movingCoordinates;
        DihedralAngleDataVector metadataVector;
        size_t currentMetadataIndex;
    };

    typedef std::vector<gmml::MolecularMetadata::GLYCAM::DihedralAngleDataVector> DihedralAngleMetadata;

    struct ResidueLinkage
    {
        ResidueLinkage(ResidueLink link_, std::vector<RotatableDihedral> dihedrals,
                       gmml::MolecularMetadata::GLYCAM::RotamerType rotamerType_, unsigned long long index_,
                       std::string name_, std::vector<Residue*> reducingOverlapResidues_,
                       std::vector<Residue*> nonReducingOverlapResidues_)
            : link(link_), rotatableDihedrals(dihedrals), rotamerType(rotamerType_), index(index_), name(name_),
              reducingOverlapResidues(reducingOverlapResidues_),
              nonReducingOverlapResidues(nonReducingOverlapResidues_) {};

        ResidueLink link;
        std::vector<RotatableDihedral> rotatableDihedrals;
        gmml::MolecularMetadata::GLYCAM::RotamerType rotamerType;
        unsigned long long index = 0;
        std::string name         = ""; // e.g. "DGalpb1-6DGlcpNAc". It being empty works with GetName();
        std::vector<cds::Residue*> reducingOverlapResidues;
        std::vector<cds::Residue*> nonReducingOverlapResidues;
    };

    std::vector<RotatableDihedral>
    rotatableDihedralsWithMultipleRotamers(const std::vector<RotatableDihedral>& dihedrals);
    size_t numberOfShapes(gmml::MolecularMetadata::GLYCAM::RotamerType rotamerType,
                          const std::vector<RotatableDihedral>& dihedrals);
    size_t numberOfLikelyShapes(gmml::MolecularMetadata::GLYCAM::RotamerType rotamerType,
                                const std::vector<RotatableDihedral>& dihedrals);
    DihedralCoordinates dihedralCoordinates(const cds::RotatableDihedral& dihedral);
    std::string print(const ResidueLinkage& linkage);
    std::string print(const RotatableDihedral& dihedral);

} // namespace cds
#endif
