#include "includes/CentralDataStructure/Overlaps/cdsOverlaps.hpp"
#include "includes/CodeUtils/constants.hpp" // maxcutoff
#include "includes/CodeUtils/logging.hpp"
#include "includes/CentralDataStructure/Measurements/measurements.hpp"
#include "includes/CentralDataStructure/Selections/residueSelections.hpp"

#include <sstream>
#include <iostream>
#include <cassert>

namespace
{
    void setResidueGeometricCenter(std::vector<cds::ResidueAtomOverlapInput>& residues)
    {
        for (auto& a : residues)
        {
            if (!a.isFixed)
            {
                a.geometricCenter = cds::calculateGeometricCenter(a.coordinates);
            }
        }
    }

    std::vector<cds::ResidueAtomOverlapInput> toResidueOverlapInput(bool mostlyFixed,
                                                                    const std::vector<cds::Residue*>& residues)
    {
        std::vector<cds::ResidueAtomOverlapInput> res;
        res.reserve(residues.size());
        for (auto& a : residues)
        {
            auto& entry = res.emplace_back(
                cds::ResidueAtomOverlapInput {mostlyFixed, false, cds::Coordinate(0.0, 0.0, 0.0), a->getCoordinates()});
            entry.geometricCenter = cds::calculateGeometricCenter(entry.coordinates);
        }
        res[0].isPartOfDihedral = true;
        res[0].isFixed          = false;
        return res;
    }
} // namespace

void cds::setGeometricCenters(cds::ResidueAtomOverlapInputPair& pair)
{
    setResidueGeometricCenter(pair.first);
    setResidueGeometricCenter(pair.second);
}

cds::ResidueAtomOverlapInputPair cds::toResidueAtomOverlapInput(const std::vector<Residue*>& residuesA,
                                                                const std::vector<Residue*>& residuesB,
                                                                bool assumeFirstSetStaysFixed)
{
    cds::ResidueAtomOverlapInputPair res {toResidueOverlapInput(assumeFirstSetStaysFixed, residuesA),
                                          toResidueOverlapInput(false, residuesB)};
    return res;
}

unsigned int cds::CountOverlappingAtoms(const std::vector<cds::Atom*>& atomsA, const std::vector<cds::Atom*>& atomsB)
{
    unsigned int overlapCount = 1;
    for (auto& atomA : atomsA)
    {
        for (auto& atomB : atomsB)
        {
            overlapCount += cds::withinDistance(constants::maxCutOff, *atomA->getCoordinate(), *atomB->getCoordinate());
        }
    }
    return overlapCount;
}

unsigned int cds::CountOverlappingAtoms(const std::vector<ResidueAtomOverlapInput>& residuesA,
                                        const std::vector<ResidueAtomOverlapInput>& residuesB)
{
    unsigned int overlapCount = 1;
    for (auto& residueA : residuesA)
    {
        for (auto& residueB : residuesB)
        {
            if (!(residueA.isPartOfDihedral && residueB.isPartOfDihedral) &&
                cds::withinDistance(constants::residueDistanceOverlapCutoff, residueA.geometricCenter,
                                    residueB.geometricCenter))
            {
                overlapCount += cds::CountOverlappingCoordinates(residueA.coordinates, residueB.coordinates);
            }
        }
    }
    return overlapCount;
}

unsigned int cds::CountOverlappingCoordinates(const std::vector<cds::Coordinate*>& coordsA,
                                              const std::vector<cds::Coordinate*>& coordsB)
{
    unsigned int overlapCount = 1;
    for (auto& coordA : coordsA)
    {
        for (auto& coordB : coordsB)
        {
            overlapCount += cds::withinDistance(constants::maxCutOff, *coordA, *coordB);
        }
    }
    return overlapCount;
}
