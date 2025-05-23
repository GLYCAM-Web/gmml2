#include "includes/CentralDataStructure/InternalPrograms/WiggleToSite/wiggleToSite.hpp"

#include "includes/CentralDataStructure/Selections/templatedSelections.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbFile.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbModel.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbSelections.hpp" //select
#include "includes/CentralDataStructure/Editors/superimposition.hpp"
#include "includes/CentralDataStructure/Geometry/overlap.hpp"
#include "includes/CentralDataStructure/Overlaps/atomOverlaps.hpp"
#include "includes/CentralDataStructure/Parameters/parameterManager.hpp"
#include "includes/CentralDataStructure/cdsFunctions/cdsFunctions.hpp"
#include "includes/CentralDataStructure/CondensedSequence/sequenceParser.hpp"
#include "includes/CentralDataStructure/Shapers/residueLinkageFunctions.hpp"
#include "includes/CentralDataStructure/Shapers/residueLinkageCreation.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralShape.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralAngleSearch.hpp"
#include "includes/MolecularMetadata/GLYCAM/dihedralangledata.hpp"
#include "includes/MolecularMetadata/elements.hpp"
#include "includes/External_Libraries/PCG/pcg_random.h"
#include "includes/CodeUtils/constants.hpp"
#include "includes/CodeUtils/containerTypes.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/metropolisCriterion.hpp"
#include "includes/CodeUtils/random.hpp"

// Prototype: Working and producing useful data in 1.5 days. Included fixing some things in the CDS.
using cds::Atom;
using gmmlPrograms::WiggleToSite;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
WiggleToSite::WiggleToSite(const cdsParameters::ParameterManager& parameterManager, WiggleToSiteInputs inputStruct)
    : substrate_(inputStruct.substrateFile_, {pdb::InputType::modelsAsMolecules, false}),
      carbohydrate_(parameterManager, MolecularMetadata::vanDerWaalsRadii(),
                    cdsCondensedSequence::parseAndReorder(inputStruct.carbohydrateSequence_))
{
    this->getCarbohydrate().Generate3DStructureFiles("./", "initial", {});
    pdb::PdbData& pdbData = this->getSubstrate().data;
    size_t superimpositionTargetId =
        pdb::residueSelector(pdbData, inputStruct.superimpositionTargetResidue_, inputStruct.substrateModelNumber_);
    cds::Residue* superimpositionTarget = pdbData.objects.residues[superimpositionTargetId];
    size_t wigglingTargetId =
        pdb::residueSelector(pdbData, inputStruct.wigglingTargetResidue_, inputStruct.substrateModelNumber_);
    cds::Residue* wigglingTarget = pdbData.objects.residues[wigglingTargetId];
    if (superimpositionTarget == nullptr || wigglingTarget == nullptr)
    {
        std::stringstream ss;
        ss << "Selection for superimposition target: " << inputStruct.superimpositionTargetResidue_
           << "\nOr selection for wiggling target: " << inputStruct.wigglingTargetResidue_ << " was not found\n";
        throw std::runtime_error(ss.str());
    }
    auto atoms                                           = getCarbohydrate().mutableAtoms();
    std::vector<cds::Coordinate> carbohydrateCoordinates = cds::atomCoordinates(atoms);
    std::vector<cds::Residue*> residues                  = this->getCarbohydrate().getResidues();
    std::vector<uint> residueNumbers                     = cds::residueNumbers(residues);
    size_t superimposedIndex = codeUtils::indexOf(residueNumbers, inputStruct.carbohydrateSuperimpositionResidue_);
    size_t wiggleIndex       = codeUtils::indexOf(residueNumbers, inputStruct.carbohydrateWigglingResidue_);
    if (superimposedIndex == residues.size())
    {
        throw std::runtime_error("Requested residue number not found in structure: " +
                                 std::to_string(inputStruct.carbohydrateSuperimpositionResidue_));
    }
    if (wiggleIndex == residues.size())
    {
        throw std::runtime_error("Requested residue number not found in structure: " +
                                 std::to_string(inputStruct.carbohydrateWigglingResidue_));
    }
    Residue* superimposeMe                                      = residues[superimposedIndex];
    Residue* wiggleMe                                           = residues[wiggleIndex];
    const GlycamMetadata::DihedralAngleDataTable& metadataTable = GlycamMetadata::dihedralAngleDataTable();
    this->superimpose(carbohydrateCoordinates, superimpositionTarget, superimposeMe);
    cds::setAtomCoordinates(atoms, carbohydrateCoordinates);
    this->getCarbohydrate().Generate3DStructureFiles("./", "superimposed", {});
    this->determineWiggleLinkages(metadataTable, superimposeMe, wiggleMe);
    std::vector<cds::Atom*> substrateWithoutSuperimpositionAtoms = codeUtils::findElementsNotInVector(
        pdb::getAtoms(this->getSubstrate().getAssemblies()), superimpositionTarget->getAtoms());
    std::vector<cds::Atom*> substrateAtomsToAvoidOverlappingWith =
        codeUtils::findElementsNotInVector(substrateWithoutSuperimpositionAtoms, wigglingTarget->getAtoms());
    this->atomsToAvoid_                                 = substrateAtomsToAvoidOverlappingWith;
    const codeUtils::SparseVector<double>& elementRadii = MolecularMetadata::vanDerWaalsRadii();
    this->setCurrentOverlap(cds::CountOverlappingAtoms(elementRadii, constants::overlapTolerance, atomsToAvoid_,
                                                       this->getCarbohydrate().getAtoms()));
    this->wiggleMeAtoms_     = {wiggleMe->FindAtom("C1"), wiggleMe->FindAtom("C3"), wiggleMe->FindAtom("C5")};
    this->wiggleTargetAtoms_ = {wigglingTarget->FindAtom("C1"), wigglingTarget->FindAtom("C3"),
                                wigglingTarget->FindAtom("C5")};
    cds::Atom* nullAtom      = nullptr;
    if (codeUtils::contains(wiggleMeAtoms_, nullAtom) || codeUtils::contains(wiggleTargetAtoms_, nullAtom))
    {
        throw std::runtime_error("Did not find the cooordinates of the atoms required for wiggling\n");
    }
    this->setCurrentDistance(this->calculateDistance());
    int structureCount =
        this->minimizeDistance(elementRadii, metadataTable, inputStruct.persistCycles_, !inputStruct.isDeterministic_);
    this->minimizeDistance(elementRadii, metadataTable, inputStruct.persistCycles_, false, structureCount);
    this->getCarbohydrate().Generate3DStructureFiles("./", "finished", {});
}

int WiggleToSite::minimizeDistance(const codeUtils::SparseVector<double>& elementRadii,
                                   const GlycamMetadata::DihedralAngleDataTable& metadataTable, int persistCycles,
                                   bool useMonteCarlo, int structureCount)
{
    uint64_t seed = codeUtils::generateRandomSeed();
    pcg32 rng(seed);
    auto randomMetadata =
        [&rng](const GlycamMetadata::DihedralAngleDataTable& table, const std::vector<size_t>& indices)
    {
        return codeUtils::weightedRandomOrder(rng, codeUtils::indicesToValues(table.weights, indices));
    };
    double angleStandardDeviation = 2.0;
    auto randomAngle              = [&rng, &angleStandardDeviation](GlycamMetadata::DihedralAngleData metadata)
    {
        auto random = [&rng, &metadata](double stdCutoff, double lower, double upper)
        {
            double num = codeUtils::normalDistributionRandomDoubleWithCutoff(rng, -stdCutoff, stdCutoff);
            return metadata.default_angle + num * (num < 0 ? lower : upper);
        };
        std::function<double(const GlycamMetadata::AngleLimit&)> onLimit = [&](const GlycamMetadata::AngleLimit& dev)
        {
            return random(1.0, dev.lowerDeviationLimit, dev.upperDeviationLimit);
        };
        std::function<double(const GlycamMetadata::AngleStd&)> onStd = [&](const GlycamMetadata::AngleStd& dev)
        {
            return random(angleStandardDeviation, dev.lowerDeviationStd, dev.upperDeviationStd);
        };
        return cds::onAngleDeviation(onLimit, onStd, metadata.angle_deviation);
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
            auto recordedShape = cds::currentShape(metadataTable, linkage.rotatableDihedrals, linkage.dihedralMetadata);
            cds::setShapeToPreference(linkage,
                                      cds::linkageShapePreference(randomMetadata, randomAngle, metadataTable,
                                                                  linkage.rotamerType, linkage.dihedralMetadata));
            double acceptance = codeUtils::uniformRandomDoubleWithinRange(rng, 0, 1);
            if (this->acceptDistance(useMonteCarlo, acceptance) && this->acceptOverlaps(elementRadii))
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
void WiggleToSite::superimpose(std::vector<cds::Coordinate>& carbohydrateCoordinates,
                               const Residue* superimpositionTarget, Residue* superimposeMe)
{
    // Limiting the selection to just these atoms as sometimes hydrogens or an oxygen is missing from xtal. That's ok.
    std::vector<cds::Atom*> superimposeMeAtoms = {superimposeMe->FindAtom("C1"), superimposeMe->FindAtom("C3"),
                                                  superimposeMe->FindAtom("C5")};
    std::vector<cds::Atom*> superTargetAtoms   = {superimpositionTarget->FindAtom("C1"),
                                                  superimpositionTarget->FindAtom("C3"),
                                                  superimpositionTarget->FindAtom("C5")};
    std::vector<cds::Coordinate> superimposeMeCoordinates = cds::atomCoordinates(superimposeMeAtoms);
    std::vector<cds::Coordinate> superTargetCoordinates   = cds::atomCoordinates(superTargetAtoms);
    cds::Superimpose(superimposeMeCoordinates, superTargetCoordinates,
                     carbohydrateCoordinates); // "alsoMoving" are the carbohydrate Coordinates
    cds::setAtomCoordinates(superimposeMeAtoms, superimposeMeCoordinates);
    cds::setAtomCoordinates(superTargetAtoms, superTargetCoordinates);
}

std::vector<cds::ResidueLinkage>&
WiggleToSite::determineWiggleLinkages(const GlycamMetadata::DihedralAngleDataTable& metadataTable,
                                      Residue* startResidue, Residue* endResidue)
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
            wiggleLinkages_.emplace_back(cds::createResidueLinkage(metadataTable, link));
        }
        previousResidue = residue;
    }
    std::cout << "\nLinkages I behold:\n" << std::endl;
    for (auto& linkage : this->getWiggleLinkages())
    {
        std::cout << linkage.name << ": "
                  << cds::numberOfShapes(metadataTable, linkage.rotamerType, linkage.dihedralMetadata) << std::endl;
    }
    return this->getWiggleLinkages();
}

double WiggleToSite::calculateDistance()
{
    return distance(wiggleTargetAtoms_[0]->coordinate(), wiggleMeAtoms_[0]->coordinate());
}

bool WiggleToSite::acceptOverlaps(const codeUtils::SparseVector<double>& elementRadii)
{
    double overlapCount = cds::CountOverlappingAtoms(elementRadii, constants::overlapTolerance, atomsToAvoid_,
                                                     getCarbohydrate().getAtoms());
    if (cds::compareOverlaps(overlapCount, this->getCurrentOverlap()) > 0)
    {
        return false;
    }
    this->setCurrentOverlap(overlapCount);
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
