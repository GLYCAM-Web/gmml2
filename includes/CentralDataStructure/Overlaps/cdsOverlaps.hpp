#ifndef GMML_INCLUDES_CENTRALDATASTRUCTURE_OVERLAPS_OVERLAPS_HPP
#define GMML_INCLUDES_CENTRALDATASTRUCTURE_OVERLAPS_OVERLAPS_HPP

#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/residue.hpp"

#include <vector>
#include <utility>

namespace cds
{
    struct ResidueAtomOverlapInput
    {
        bool isFixed;
        bool isPartOfDihedral;
        Coordinate geometricCenter;
        std::vector<Coordinate*> coordinates;
    };

    typedef std::pair<std::vector<ResidueAtomOverlapInput>, std::vector<ResidueAtomOverlapInput>>
        ResidueAtomOverlapInputPair;

    void setGeometricCenters(cds::ResidueAtomOverlapInputPair& pair);

    ResidueAtomOverlapInputPair toResidueAtomOverlapInput(const std::vector<Residue*>& residuesA,
                                                          const std::vector<Residue*>& residuesB,
                                                          bool assumeFirstSetStaysFixed);

    unsigned int CountOverlappingAtoms(const std::vector<Atom*>& atomsA, const std::vector<Atom*>& atomsB);
    unsigned int CountOverlappingAtoms(const std::vector<ResidueAtomOverlapInput>& residuesA,
                                       const std::vector<ResidueAtomOverlapInput>& residuesB);
    unsigned int CountOverlappingCoordinates(const std::vector<Coordinate*>& coordsA,
                                             const std::vector<Coordinate*>& coordsB);
} // namespace cds
#endif
