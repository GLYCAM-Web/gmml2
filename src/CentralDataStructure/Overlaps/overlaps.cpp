#include "includes/CentralDataStructure/Overlaps/overlaps.hpp"
#include "includes/CodeUtils/constants.hpp" // maxcutoff
#include "includes/CentralDataStructure/Measurements/measurements.hpp"

#include <numeric>

namespace
{
    unsigned int coordinateOverlaps(const std::vector<cds::Coordinate>& coordsA,
                                    const std::vector<cds::Coordinate>& coordsB)
    {
        unsigned int count = 0;
        for (auto& a : coordsA)
        {
            for (auto& b : coordsB)
            {
                count += cds::withinDistance(constants::maxCutOff, a, b);
            }
        }
        return count;
    }

    void insertCoordinatesWithinDistance(std::vector<cds::Coordinate>& result, double distance,
                                         const cds::Coordinate& pt, const std::vector<cds::Coordinate>& coords,
                                         const std::pair<size_t, size_t>& range)
    {
        result.clear();
        for (size_t n = range.first; n < range.second; n++)
        {
            auto& a = coords[n];
            if (cds::withinDistance(distance, pt, a))
            {
                result.push_back(a);
            }
        }
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
    std::vector<Coordinate> coordsA;
    std::vector<Coordinate> coordsB;
    double residueCutoff      = constants::residueDistanceOverlapCutoff + constants::maxCutOff;
    double radialCutoff       = 0.5 * constants::residueDistanceOverlapCutoff + constants::maxCutOff;
    unsigned int overlapCount = 0;
    for (size_t n = 0; n < mostlyFixed.geometricCenters.size(); n++)
    {
        auto& centerFixed = mostlyFixed.geometricCenters[n];
        for (size_t k = 0; k < moving.geometricCenters.size(); k++)
        {
            auto& centerMoving = moving.geometricCenters[k];
            if ((n > 0 || k > 0) && cds::withinDistance(residueCutoff, centerFixed, centerMoving))
            {
                insertCoordinatesWithinDistance(coordsA, radialCutoff, centerMoving, mostlyFixed.atomCoordinates,
                                                mostlyFixed.residueAtoms[n]);
                insertCoordinatesWithinDistance(coordsB, radialCutoff, centerFixed, moving.atomCoordinates,
                                                moving.residueAtoms[k]);
                overlapCount += coordinateOverlaps(coordsA, coordsB);
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
