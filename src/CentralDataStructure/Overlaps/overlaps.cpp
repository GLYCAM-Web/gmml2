#include "includes/CentralDataStructure/Overlaps/overlaps.hpp"
#include "includes/CodeUtils/constants.hpp" // maxcutoff
#include "includes/CentralDataStructure/Measurements/measurements.hpp"

#include <numeric>

namespace
{
    unsigned int coordinateOverlaps(const std::vector<cds::Coordinate>& coordsA,
                                    const std::vector<cds::Coordinate>& coordsB, const std::pair<size_t, size_t> rangeA,
                                    const std::pair<size_t, size_t> rangeB)
    {
        unsigned int count = 0;
        for (size_t n = rangeA.first; n < rangeA.second; n++)
        {
            for (size_t k = rangeB.first; k < rangeB.second; k++)
            {
                count += cds::withinDistance(constants::maxCutOff, coordsA[n], coordsB[k]);
            }
        }
        return count;
    }
} // namespace

cds::ResidueAtomOverlapInput cds::toOverlapInput(const std::vector<Residue*>& residues)
{
    size_t atomCount = 0;
    for (auto res : residues)
    {
        atomCount += res->atomCount();
    }
    std::vector<Coordinate> coordinates;
    coordinates.reserve(atomCount);
    std::vector<std::pair<size_t, size_t>> residueAtoms;
    residueAtoms.reserve(residues.size());
    size_t currentAtom = 0;
    for (auto& res : residues)
    {
        size_t startAtom = currentAtom;
        res->insertCoordinatesInto(coordinates);
        currentAtom += res->atomCount();
        residueAtoms.push_back({startAtom, currentAtom});
    }
    std::vector<Coordinate> geometricCenters;
    geometricCenters.reserve(residues.size());
    for (size_t n = 0; n < residueAtoms.size(); n++)
    {
        auto range        = residueAtoms[n];
        Coordinate center = std::accumulate(coordinates.begin() + range.first, coordinates.begin() + range.second,
                                            Coordinate(0.0, 0.0, 0.0));
        geometricCenters.push_back(scaleBy(1.0 / (range.second - range.first), center));
    }
    return {coordinates, geometricCenters, residueAtoms};
}

unsigned int cds::CountOverlappingAtoms(const std::vector<cds::Atom*>& atomsA, const std::vector<cds::Atom*>& atomsB)
{
    unsigned int overlapCount = 0;
    for (auto& atomA : atomsA)
    {
        for (auto& atomB : atomsB)
        {
            overlapCount += cds::withinDistance(constants::maxCutOff, *atomA->getCoordinate(), *atomB->getCoordinate());
        }
    }
    return overlapCount;
}

unsigned int cds::CountOverlappingAtoms(const ResidueAtomOverlapInputReference& mostlyFixed,
                                        const ResidueAtomOverlapInputReference& moving)
{
    unsigned int overlapCount = 0;
    for (size_t n = 0; n < mostlyFixed.geometricCenters.size(); n++)
    {
        for (size_t k = 0; k < moving.geometricCenters.size(); k++)
        {
            if ((n > 0 || k > 0) && cds::withinDistance(constants::residueDistanceOverlapCutoff + constants::maxCutOff,
                                                        mostlyFixed.geometricCenters[n], moving.geometricCenters[k]))
            {
                overlapCount += coordinateOverlaps(mostlyFixed.atomCoordinates, moving.atomCoordinates,
                                                   mostlyFixed.residueAtoms[n], moving.residueAtoms[k]);
            }
        }
    }
    return overlapCount;
}

unsigned int cds::CountOverlappingAtoms(const std::vector<Residue*>& residuesA, const std::vector<Residue*>& residuesB)
{
    auto inputA = toOverlapInput(residuesA);
    auto inputB = toOverlapInput(residuesB);

    return CountOverlappingAtoms({inputA.atomCoordinates, inputA.geometricCenters, inputA.residueAtoms},
                                 {inputB.atomCoordinates, inputB.geometricCenters, inputB.residueAtoms});
}

unsigned int cds::CountOverlappingCoordinates(const std::vector<cds::Coordinate*>& coordsA,
                                              const std::vector<cds::Coordinate*>& coordsB)
{
    unsigned int overlapCount = 0;
    for (auto& coordA : coordsA)
    {
        for (auto& coordB : coordsB)
        {
            overlapCount += cds::withinDistance(constants::maxCutOff, *coordA, *coordB);
        }
    }
    return overlapCount;
}