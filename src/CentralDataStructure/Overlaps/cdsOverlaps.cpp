#include "includes/CentralDataStructure/Overlaps/cdsOverlaps.hpp"
#include "includes/CodeUtils/constants.hpp" // maxcutoff
#include "includes/CentralDataStructure/Measurements/measurements.hpp"
#include "includes/CentralDataStructure/Selections/residueSelections.hpp"

#include <cassert>

namespace
{
    void setResidueGeometricCenter(std::vector<cds::ResidueAtomOverlapInput>& residues)
    {
        for (auto& a : residues)
        {
            a.geometricCenter = cds::calculateGeometricCenter(a.coordinates);
        }
    }

    std::vector<cds::ResidueAtomOverlapInput> toResidueOverlapInput(const std::vector<cds::Residue*>& residues)
    {
        std::vector<cds::ResidueAtomOverlapInput> res;
        res.reserve(residues.size());
        for (auto& a : residues)
        {
            res.emplace_back(cds::ResidueAtomOverlapInput {
                false, Coordinate {0.0, 0.0, 0.0},
                  a->getCoordinates()
            });
        }
        return res;
    }

    void setNeighbors(const std::vector<cds::Residue*>& residuesA, const std::vector<cds::Residue*>& residuesB,
                      std::vector<cds::ResidueAtomOverlapInput>& inputA,
                      std::vector<cds::ResidueAtomOverlapInput>& inputB)
    {
        int neighborCount = 0;
        for (size_t n = 0; n < residuesA.size(); n++)
        {
            for (size_t k = 0; k < residuesB.size(); k++)
            {
                if (cdsSelections::areNeighbors(residuesA[n], residuesB[k]))
                {
                    inputA[n].isPartOfDihedral = true;
                    inputB[k].isPartOfDihedral = true;
                    neighborCount++;
                }
            }
        }
        assert(neighborCount <= 1);
    }
} // namespace

void cds::setGeometricCenters(cds::ResidueAtomOverlapInputPair& pair)
{
    setResidueGeometricCenter(pair.first);
    setResidueGeometricCenter(pair.second);
}

cds::ResidueAtomOverlapInputPair cds::toResidueAtomOverlapInput(const std::vector<Residue*>& residuesA,
                                                                const std::vector<Residue*>& residuesB)
{
    cds::ResidueAtomOverlapInputPair res {toResidueOverlapInput(residuesA), toResidueOverlapInput(residuesB)};
    setNeighbors(residuesA, residuesB, res.first, res.second);
    setGeometricCenters(res);
    return res;
}

double cds::CalculateAtomicOverlaps(cds::Atom* atomA, cds::Atom* atomB, double radiusA, double radiusB)
{
    double distance = atomA->getCoordinate()->Distance(atomB->getCoordinate());
    if (radiusA == 0.0) // default value is 0.0, but user can provide.
    { // element info not usually set, so I look at first letter of atom name. This may be why you're reading this.
        if (atomA->getName().at(0) == 'C')
        {
            radiusA = 1.70; // Rowland and Taylor modification to vdW.
        }
        if (atomA->getName().at(0) == 'O')
        {
            radiusA = 1.52;
        }
        if (atomA->getName().at(0) == 'N')
        {
            radiusA = 1.55;
        }
        if (atomA->getName().at(0) == 'S')
        {
            radiusA = 1.80;
        }
        if (atomA->getName().at(0) == 'P')
        {
            radiusA = 1.80;
        }
        if (atomA->getName().at(0) == 'H')
        {
            radiusA = 1.09;
        }
    }
    if (radiusB == 0.0)
    {
        if (atomB->getName().at(0) == 'C')
        {
            radiusB = 1.70;
        }
        if (atomB->getName().at(0) == 'O')
        {
            radiusB = 1.52;
        }
        if (atomB->getName().at(0) == 'N')
        {
            radiusB = 1.55;
        }
        if (atomB->getName().at(0) == 'S')
        {
            radiusB = 1.80;
        }
        if (atomB->getName().at(0) == 'P')
        {
            radiusB = 1.80;
        }
        if (atomB->getName().at(0) == 'H')
        {
            radiusB = 1.09;
        }
    }
    //    std::cout << "Distance: " << distance << " radiusA: " << radiusA << " radiusB: " << radiusB << std::endl;
    double overlap = 0.0;
    if (radiusA + radiusB > distance) // Close enough to overlap
    {
        if (std::abs(radiusA - radiusB) > distance) // If one sphere is completely inside the other
        {                                           // then calculate the surface area of the buried (smaller) sphere.
            if (radiusA < radiusB)
            {
                overlap = 4 * constants::PI_RADIAN * (radiusA * radiusA);
            }
            else
            {
                overlap = 4 * constants::PI_RADIAN * (radiusB * radiusB);
            }
        }
        else // Normal situation, partial overlap we need to calculate
        {    // Eqn 1, Rychkov and Petukhov, J. Comput. Chem., 2006, Joint Neighbours. Each atom against each atom, so
             // overlap can be "double" counted. See paper.
            overlap = (2 * (constants::PI_RADIAN)*radiusA *
                       (radiusA - distance / 2 - (((radiusA * radiusA) - (radiusB * radiusB)) / (2 * distance))));
        }
    }
    if ((overlap < 0.0) || (radiusA == -0.1) || (radiusB == -0.1))
    { // Either the user didn't specify the radius or the element isn't one of the above
        //        std::cout << "Neggie: " << overlap << " d: " << distance << ", A: " << atomA->getName() << ", rA: " <<
        //        radiusA << ", B: " << atomB->getName() << ", rB: " << radiusB << std::endl;
        return 0.0; // negative overlap isn't a thing.
    }
    //    std::cout << "Non-normalized Overlap=" << totalOverlap << std::endl;
    return overlap;
}

double cds::CalculateAtomicOverlaps(std::vector<cds::Atom*> atomsA, std::vector<cds::Atom*> atomsB, bool print)
{
    double totalOverlap   = 0.0;
    double currentOverlap = 0.0;
    for (auto& atomA : atomsA)
    {
        for (auto& atomB : atomsB)
        { // if not the same atom (index is unique)
            if ((atomA->getIndex() != atomB->getIndex()) &&
                (cds::CheckIfOtherCoordinateIsWithinDistance(atomA->getCoordinate(), atomB->getCoordinate(),
                                                             constants::maxCutOff * 2)))
            {
                currentOverlap = cds::CalculateAtomicOverlaps(atomA, atomB);
                totalOverlap   += currentOverlap;
                if (print)
                {
                    //                    std::cout << atomA->getId() << "::" << atomB->getId() << ": " <<
                    //                    (currentOverlap / constants::CARBON_SURFACE_AREA) << "\n";
                }
            }
        }
    }
    if (totalOverlap < 0.0) // Negative number fail
    {
        std::stringstream ss;
        ss << "Negative overlap should not happen, this is a bug: " << totalOverlap;
        gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
        throw std::runtime_error(ss.str());
    }
    return (totalOverlap / constants::CARBON_SURFACE_AREA); // Normalise to area of a buried carbon
}

double cds::CalculateAtomicOverlapsBetweenNonBondedAtoms(std::vector<cds::Atom*>& atomsA,
                                                         std::vector<cds::Atom*>& atomsB)
{
    double totalOverlap = 0.0;
    for (auto& atomA : atomsA)
    {
        for (auto& atomB : atomsB)
        {
            bool isNeighbor                   = false;
            std::vector<cds::Atom*> neighbors = atomA->getNeighbors();
            for (auto& neighbor : neighbors)
            {
                if (atomB->getIndex() == neighbor->getIndex())
                {
                    isNeighbor = true;
                }
            }
            if ((isNeighbor == false) && (atomA->getIndex() != atomB->getIndex()) &&
                (cds::CheckIfOtherCoordinateIsWithinDistance(atomA->getCoordinate(), atomB->getCoordinate(),
                                                             constants::maxCutOff)))
            {
                totalOverlap += cds::CalculateAtomicOverlaps(atomA, atomB);
            }
        }
    }
    return (totalOverlap / constants::CARBON_SURFACE_AREA); // Normalise to area of a buried carbon
}

unsigned int cds::CountOverlappingResidues(const std::vector<cds::Residue*>& residuesA,
                                           const std::vector<cds::Residue*>& residuesB)
{
    unsigned int overlapCount = 1;
    //    std::cout << "Number of A : B residues is " << residuesA.size() << " : " << residuesB.size() << std::endl <<
    //    std::flush;
    for (auto& residueA : residuesA)
    {
        const Coordinate* residueA_Center = residueA->calculateGeometricCenter();
        for (auto& residueB : residuesB)
        {
            overlapCount += cds::CheckIfOtherCoordinateIsWithinDistance(residueA_Center, residueB->getGeometricCenter(),
                                                                        constants::residueDistanceOverlapCutoff);
        }
    }
    return overlapCount;
}

unsigned int cds::CountOverlappingAtoms(const std::vector<cds::Atom*>& atomsA, const std::vector<cds::Atom*>& atomsB)
{
    unsigned int overlapCount = 1;
    for (auto& atomA : atomsA)
    {
        for (auto& atomB : atomsB)
        {
            overlapCount += cds::CheckIfOtherCoordinateIsWithinDistance(atomA->getCoordinate(), atomB->getCoordinate(),
                                                                        constants::maxCutOff);
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
                cds::CheckIfOtherCoordinateIsWithinDistance(&residueA.geometricCenter, &residueB.geometricCenter,
                                                            constants::residueDistanceOverlapCutoff))
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
            overlapCount += cds::CheckIfOtherCoordinateIsWithinDistance(coordA, coordB, constants::maxCutOff);
        }
    }
    return overlapCount;
}
