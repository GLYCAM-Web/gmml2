#ifndef INCLUDES_CENTRALDATASTRUCTURE_OVERLAPS_ATOMOVERLAPS_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_OVERLAPS_ATOMOVERLAPS_HPP

#include "includes/CentralDataStructure/Geometry/types.hpp"
#include "includes/CentralDataStructure/Geometry/overlap.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/residue.hpp"

#include <vector>
#include <utility>

namespace cds
{
    struct ResiduesWithOverlapWeight
    {
        std::vector<Residue*> residues;
        std::vector<double> weights;
    };

    struct ResidueAtomOverlapInput
    {
        std::vector<Sphere> atomCoordinates;
        std::vector<Sphere> boundingSpheres;
        const std::vector<std::pair<size_t, size_t>> residueAtoms;
        const std::vector<double> residueWeights;
    };

    struct ResidueAtomOverlapInputReference
    {
        std::vector<Sphere>& atomCoordinates;
        std::vector<Sphere>& boundingSpheres;
        const std::vector<std::pair<size_t, size_t>>& residueAtoms;
        const std::vector<double>& residueWeights;
    };

    ResidueAtomOverlapInput toOverlapInput(const ResiduesWithOverlapWeight& input);
    Overlap CountOverlappingAtoms(const std::vector<Atom*>& atomsA, const std::vector<Atom*>& atomsB);
    Overlap CountOverlappingAtoms(bool ignoreNeighboringResidues, const ResidueAtomOverlapInputReference& mostlyFixed,
                                  const ResidueAtomOverlapInputReference& moving);
    Overlap CountOverlappingAtoms(bool ignoreNeighboringResidues, const ResiduesWithOverlapWeight& residuesA,
                                  const ResiduesWithOverlapWeight& residuesB);
} // namespace cds
#endif
