#ifndef GMML_INCLUDES_CENTRALDATASTRUCTURE_OVERLAPS_OVERLAPS_HPP
#define GMML_INCLUDES_CENTRALDATASTRUCTURE_OVERLAPS_OVERLAPS_HPP

#include "includes/CentralDataStructure/boundingSphere.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/residue.hpp"

#include <vector>
#include <utility>

namespace cds
{
    struct ResidueAtomOverlapInput
    {
        std::vector<Coordinate> atomCoordinates;
        std::vector<Sphere> boundingSpheres;
        const std::vector<std::pair<size_t, size_t>> residueAtoms;
    };

    struct ResidueAtomOverlapInputReference
    {
        std::vector<Coordinate>& atomCoordinates;
        std::vector<Sphere> boundingSpheres;
        const std::vector<std::pair<size_t, size_t>>& residueAtoms;
    };

    ResidueAtomOverlapInput toOverlapInput(const std::vector<Residue*>& residues);
    unsigned int CountOverlappingAtoms(const std::vector<Atom*>& atomsA, const std::vector<Atom*>& atomsB);
    unsigned int CountOverlappingAtoms(const ResidueAtomOverlapInputReference& mostlyFixed,
                                       const ResidueAtomOverlapInputReference& moving);
    unsigned int CountOverlappingAtoms(const std::vector<Residue*>& residuesA, const std::vector<Residue*>& residuesB);
    unsigned int CountOverlappingCoordinates(const std::vector<Coordinate*>& coordsA,
                                             const std::vector<Coordinate*>& coordsB);
} // namespace cds
#endif
