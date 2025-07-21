#ifndef INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_RESIDUELINKAGETYPES_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_SHAPERS_RESIDUELINKAGETYPES_HPP

#include "include/CentralDataStructure/cdsTypes.hpp"
#include "include/metadata/dihedralangledata.hpp"

#include <array>
#include <string>
#include <utility>
#include <vector>

namespace gmml
{
    struct ResidueLink
    {
        std::pair<Residue*, Residue*> residues;
        std::pair<Atom*, Atom*> atoms;
    };

    struct ResidueLinkAttributes
    {
        std::pair<ResidueAttributes, ResidueAttributes> residues;
        std::pair<std::string, std::string> atoms;
    };

    typedef std::array<Atom*, 4> DihedralAtoms;

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
        RotamerType rotamerType;
        bool isDerivative;
        unsigned long long index = 0;
        std::string name = ""; // e.g. "DGalpb1-6DGlcpNAc". It being empty works with GetName();
    };
} // namespace gmml
#endif
