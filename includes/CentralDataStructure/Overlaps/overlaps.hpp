#ifndef INCLUDES_CENTRALDATASTRUCTURE_OVERLAPS_OVERLAPS_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_OVERLAPS_OVERLAPS_HPP

#include "includes/CentralDataStructure/Geometry/boundingSphere.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/residue.hpp"

#include <vector>
#include <utility>

namespace cds
{
    struct ResidueAtomOverlapInput
    {
        std::vector<Sphere> atomCoordinates;
        std::vector<Sphere> boundingSpheres;
        const std::vector<std::pair<size_t, size_t>> residueAtoms;
    };

    struct ResidueAtomOverlapInputReference
    {
        std::vector<Sphere>& atomCoordinates;
        std::vector<Sphere>& boundingSpheres;
        const std::vector<std::pair<size_t, size_t>>& residueAtoms;
    };

    ResidueAtomOverlapInput toOverlapInput(const std::vector<Residue*>& residues);
    unsigned int CountOverlappingAtoms(const std::vector<Atom*>& atomsA, const std::vector<Atom*>& atomsB);
    unsigned int CountOverlappingAtoms(const ResidueAtomOverlapInputReference& mostlyFixed,
                                       const ResidueAtomOverlapInputReference& moving);
    unsigned int CountOverlappingAtoms(const std::vector<Residue*>& residuesA, const std::vector<Residue*>& residuesB);
} // namespace cds
#endif
