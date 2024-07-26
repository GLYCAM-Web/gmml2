#include "includes/CentralDataStructure/InternalPrograms/WiggleToSite/wiggleToSite.hpp"
#include "includes/CentralDataStructure/InternalPrograms/WiggleToSite/wiggleToSite.hpp"

#include "includes/CentralDataStructure/Selections/templatedSelections.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbFile.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbSelections.hpp" //select
#include "includes/CentralDataStructure/Editors/superimposition.hpp"
#include "includes/CodeUtils/metropolisCriterion.hpp"
#include "includes/CodeUtils/random.hpp"
#include "includes/CentralDataStructure/Overlaps/overlaps.hpp"
#include "includes/CentralDataStructure/Shapers/residueLinkage.hpp"
#include "includes/CentralDataStructure/Shapers/residueLinkageCreation.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralShape.hpp"
#include "includes/External_Libraries/PCG/pcg_random.h"

// Prototype: Working and producing useful data in 1.5 days. Included fixing some things in the CDS.
using cds::Atom;
using gmmlPrograms::WiggleToSite;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
WiggleToSite::WiggleToSite(WiggleToSiteInputs inputStruct)
    : substrate_(inputStruct.substrateFile_), carbohydrate_(inputStruct.carbohydrateSequence_)
{
    this->getCarbohydrate().Generate3DStructureFiles("./", "initial");
    std::vector<Coordinate*> carbohydrateCoordinates = cds::getCoordinatesFromAtoms(this->getCarbohydrate().getAtoms());
    // const Residue* superimpositionTarget = codeUtils::findElementWithNumber(this->getSubstrate().getResidues(),
    // inputStruct.superimpositionTargetResidue_);
    const Residue* superimpositionTarget             = pdb::residueSelector(
        this->getSubstrate(), inputStruct.superimpositionTargetResidue_, inputStruct.substrateModelNumber_);
    const Residue* wigglingTarget = pdb::residueSelector(this->getSubstrate(), inputStruct.wigglingTargetResidue_,
                                                         inputStruct.substrateModelNumber_);
    if (superimpositionTarget == nullptr || wigglingTarget == nullptr)
    {
        std::stringstream ss;
        ss << "Selection for superimposition target: " << inputStruct.superimpositionTargetResidue_
           << "\nOr selection for wiggling target: " << inputStruct.wigglingTargetResidue_ << " was not found\n";
        throw std::runtime_error(ss.str());
    }
    //    std::cout << "Super target is " << superimpositionTarget->getStringId() << std::endl;
    Residue* superimposeMe = codeUtils::findElementWithNumber(this->getCarbohydrate().getResidues(),
                                                              inputStruct.carbohydrateSuperimpositionResidue_);
    Residue* wiggleMe      = codeUtils::findElementWithNumber(this->getCarbohydrate().getResidues(),
                                                              inputStruct.carbohydrateWigglingResidue_);
    //    std::cout << "Residue to wiggle is " << wiggleMe->getStringId();
    this->superimpose(carbohydrateCoordinates, superimpositionTarget, superimposeMe);
    this->getCarbohydrate().Generate3DStructureFiles("./", "superimposed");
    this->determineWiggleLinkages(superimposeMe, wiggleMe);
    std::vector<cds::Atom*> substrateWithoutSuperimpositionAtoms =
        codeUtils::findElementsNotInVector(this->getSubstrate().getAtoms(), superimpositionTarget->getAtoms());
    std::vector<cds::Atom*> substrateAtomsToAvoidOverlappingWith =
        codeUtils::findElementsNotInVector(substrateWithoutSuperimpositionAtoms, wigglingTarget->getAtoms());
    this->atomsToAvoid_ = substrateAtomsToAvoidOverlappingWith;
    this->setCurrentOverlapCount(cds::CountOverlappingAtoms(atomsToAvoid_, this->getCarbohydrate().getAtoms()));
    // call below function.
    //    std::cout << "Finished reading and ready to rock captain" << std::endl;
    this->wiggleMeCoordinates_ = {wiggleMe->FindAtom("C1")->getCoordinate(), wiggleMe->FindAtom("C3")->getCoordinate(),
                                  wiggleMe->FindAtom("C5")->getCoordinate()};
    this->wiggleTargetCoordinates_ = {wigglingTarget->FindAtom("C1")->getCoordinate(),
                                      wigglingTarget->FindAtom("C3")->getCoordinate(),
                                      wigglingTarget->FindAtom("C5")->getCoordinate()};
    if (wiggleMeCoordinates_.size() < 3 || wiggleTargetCoordinates_.size() < 3)
    {
        throw std::runtime_error("Did not find the cooordinates of the atoms required for wiggling\n");
    }
    this->setCurrentDistance(this->calculateDistance());
    int structureCount = this->minimizeDistance(inputStruct.persistCycles_, !inputStruct.isDeterministic_);
    this->minimizeDistance(inputStruct.persistCycles_, false, structureCount);
    this->getCarbohydrate().Generate3DStructureFiles("./", "finished");
}

int WiggleToSite::minimizeDistance(int persistCycles, bool useMonteCarlo, int structureCount)
{
    uint64_t seed = codeUtils::generateRandomSeed();
    pcg32 rng(seed);
    auto randomMetadata = [&rng](gmml::MolecularMetadata::GLYCAM::DihedralAngleDataVector metadataVector)
    {
        return codeUtils::uniformRandomVectorIndex(rng, metadataVector);
    };
    auto randomAngle = [&rng](gmml::MolecularMetadata::GLYCAM::DihedralAngleData metadata)
    {
        auto range = cds::angleBounds(metadata);
        return codeUtils::uniformRandomDoubleWithinRange(rng, range.lower, range.upper);
    };

    int cycle = 0;
    while (cycle < persistCycles)
    {
        ++cycle;
        std::stringstream ss;
        ss << "Cycle " << cycle << "/" << persistCycles << "\n";
        gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
        for (auto& linkage : codeUtils::shuffleVector(rng, this->getWiggleLinkages()))
        {
            auto recordedShape = cds::currentShape(linkage.rotatableDihedrals);
            cds::setRandomShapeUsingMetadata(randomMetadata, randomAngle, linkage.rotatableDihedrals);
            double acceptance = codeUtils::uniformRandomDoubleWithinRange(rng, 0, 1);
            if (this->acceptDistance(useMonteCarlo, acceptance) && this->acceptOverlaps())
            {
                cycle = 0; // reset when it improves
            }
            else
            {
                cds::setShape(linkage.rotatableDihedrals, recordedShape);
            }
        }
    }
    return structureCount;
}

//////////////////////////////////////////////////////////
//                  PRIVATE FUNCTIONS                   //
//////////////////////////////////////////////////////////
void WiggleToSite::superimpose(std::vector<Coordinate*>& carbohydrateCoordinates, const Residue* superimpositionTarget,
                               Residue* superimposeMe)
{
    // Limiting the selection to just these atoms as sometimes hydrogens or an oxygen is missing from xtal. That's ok.
    std::vector<Coordinate*> superimposeMeCoordinates = {superimposeMe->FindAtom("C1")->getCoordinate(),
                                                         superimposeMe->FindAtom("C3")->getCoordinate(),
                                                         superimposeMe->FindAtom("C5")->getCoordinate()};
    std::vector<Coordinate*> superTargetCoordinates   = {superimpositionTarget->FindAtom("C1")->getCoordinate(),
                                                         superimpositionTarget->FindAtom("C3")->getCoordinate(),
                                                         superimpositionTarget->FindAtom("C5")->getCoordinate()};
    cds::Superimpose(superimposeMeCoordinates, superTargetCoordinates,
                     carbohydrateCoordinates); // "alsoMoving" are the carbohydrate Coordinates
    return;
}

std::vector<cds::ResidueLinkage>& WiggleToSite::determineWiggleLinkages(Residue* startResidue, Residue* endResidue)
{
    std::vector<Residue*> residuesInPath;
    bool targetFound = false;
    std::vector<Residue*> visitedResidues;
    std::cout << "Gonna find path between " << startResidue->getName() << " and " << endResidue->getName() << "\n";
    codeUtils::findPathBetweenElementsInGraph(startResidue, endResidue, visitedResidues, residuesInPath, targetFound);
    Residue* previousResidue = nullptr; // wanna skip the first iteration
    for (auto& residue : residuesInPath)
    {
        std::cout << residue->getName() << "_" << residue->getNumber() << ", ";
        if (previousResidue != nullptr)
        {
            cds::ResidueLink link = cds::findResidueLink({previousResidue, residue});
            wiggleLinkages_.emplace_back(cds::createResidueLinkage(link));
        }
        previousResidue = residue;
    }
    std::cout << "\nLinkages I behold:\n" << std::endl;
    for (auto& linkage : this->getWiggleLinkages())
    {
        std::cout << linkage.name << ": " << cds::numberOfShapes(linkage.rotatableDihedrals) << std::endl;
    }
    return this->getWiggleLinkages();
}

double WiggleToSite::calculateDistance()
{
    return distance(*wiggleTargetCoordinates_.at(0), *wiggleMeCoordinates_.at(0));
}

bool WiggleToSite::acceptOverlaps()
{
    cds::Overlap overlapCount = cds::CountOverlappingAtoms(atomsToAvoid_, getCarbohydrate().getAtoms());
    if (cds::compareOverlaps(overlapCount, this->getCurrentOverlapCount()) > 0)
    {
        return false;
    }
    this->setCurrentOverlapCount(overlapCount);
    return true;
}

bool WiggleToSite::acceptDistance(bool useMonteCarlo, double acceptance)
{
    double distance = this->calculateDistance();
    if (distance < this->getCurrentDistance() ||
        (useMonteCarlo &&
         monte_carlo::accept_via_metropolis_criterion(acceptance, distance - this->getCurrentDistance())))
    {
        this->setCurrentDistance(distance);
        return true;
    }
    return false;
}
