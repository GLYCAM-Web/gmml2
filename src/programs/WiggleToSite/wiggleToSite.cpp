#include "include/programs/WiggleToSite/wiggleToSite.hpp"

#include "include/CentralDataStructure/Selections/templatedSelections.hpp"
#include "include/CentralDataStructure/Shapers/dihedralAngleSearch.hpp"
#include "include/CentralDataStructure/Shapers/dihedralShape.hpp"
#include "include/CentralDataStructure/Shapers/residueLinkageCreation.hpp"
#include "include/CentralDataStructure/Shapers/residueLinkageFunctions.hpp"
#include "include/CentralDataStructure/cdsFunctions.hpp"
#include "include/external/pcg/pcg_random.h"
#include "include/geometry/overlap.hpp"
#include "include/geometry/superimposition.hpp"
#include "include/metadata/dihedralangledata.hpp"
#include "include/metadata/elements.hpp"
#include "include/pdb/pdbFile.hpp"
#include "include/pdb/pdbModel.hpp"
#include "include/readers/parameterManager.hpp"
#include "include/sequence/sequenceParser.hpp"
#include "include/structure/atomOverlaps.hpp"
#include "include/util/constants.hpp"
#include "include/util/containerTypes.hpp"
#include "include/util/logging.hpp"
#include "include/util/metropolisCriterion.hpp"
#include "include/util/random.hpp"

namespace gmml
{
    // Prototype: Working and producing useful data in 1.5 days. Included fixing some things in the CDS.

    namespace
    {
        double CountOverlappingAtoms(
            const util::SparseVector<double>& elementRadii,
            double overlapTolerance,
            const std::vector<Atom*>& atomsA,
            const std::vector<Atom*>& atomsB)
        {
            std::vector<Sphere> coordsA = atomCoordinatesWithRadii(elementRadii, atomsA);
            std::vector<Sphere> coordsB = atomCoordinatesWithRadii(elementRadii, atomsB);
            std::vector<Element> elementsA = atomElements(atomsA);
            std::vector<Element> elementsB = atomElements(atomsB);
            const PotentialTable potentialTable =
                gmml::potentialTable(elementRadii, util::vectorOr(foundElements(elementsA), foundElements(elementsB)));

            double overlap = 0.0;
            for (size_t n = 0; n < atomsA.size(); n++)
            {
                for (size_t k = 0; k < atomsB.size(); k++)
                {
                    PotentialFactor factor = potentialFactor(potentialTable, elementsA[n], elementsB[k]);
                    overlap += overlapAmount(factor, overlapTolerance, coordsA[n], coordsB[k]);
                }
            }
            return overlap;
        }
    } // namespace

    //////////////////////////////////////////////////////////
    //                       CONSTRUCTOR                    //
    //////////////////////////////////////////////////////////
    WiggleToSite::WiggleToSite(const ParameterManager& parameterManager, WiggleToSiteInputs inputStruct)
        : substrate_(inputStruct.substrateFile_, {pdb::InputType::modelsAsMolecules, false}),
          carbohydrate_(
              parameterManager, vanDerWaalsRadii(), sequence::parseAndReorder(inputStruct.carbohydrateSequence_))
    {
        this->getCarbohydrate().Generate3DStructureFiles("./", "initial", {});
        pdb::PdbData& pdbData = this->getSubstrate().data;
        size_t superimpositionTargetId =
            residueSelector(pdbData, inputStruct.superimpositionTargetResidue_, inputStruct.substrateModelNumber_);
        Residue* superimpositionTarget = pdbData.objects.residues[superimpositionTargetId];
        size_t wigglingTargetId =
            residueSelector(pdbData, inputStruct.wigglingTargetResidue_, inputStruct.substrateModelNumber_);
        Residue* wigglingTarget = pdbData.objects.residues[wigglingTargetId];
        if (superimpositionTarget == nullptr || wigglingTarget == nullptr)
        {
            std::stringstream ss;
            ss << "Selection for superimposition target: " << inputStruct.superimpositionTargetResidue_
               << "\nOr selection for wiggling target: " << inputStruct.wigglingTargetResidue_ << " was not found\n";
            throw std::runtime_error(ss.str());
        }
        auto atoms = getCarbohydrate().mutableAtoms();
        std::vector<Coordinate> carbohydrateCoordinates = atomCoordinates(atoms);
        std::vector<Residue*> residues = this->getCarbohydrate().getResidues();
        std::vector<uint> residueNumbers = gmml::residueNumbers(residues);
        size_t superimposedIndex = util::indexOf(residueNumbers, inputStruct.carbohydrateSuperimpositionResidue_);
        size_t wiggleIndex = util::indexOf(residueNumbers, inputStruct.carbohydrateWigglingResidue_);
        if (superimposedIndex == residues.size())
        {
            throw std::runtime_error(
                "Requested residue number not found in structure: " +
                std::to_string(inputStruct.carbohydrateSuperimpositionResidue_));
        }
        if (wiggleIndex == residues.size())
        {
            throw std::runtime_error(
                "Requested residue number not found in structure: " +
                std::to_string(inputStruct.carbohydrateWigglingResidue_));
        }
        Residue* superimposeMe = residues[superimposedIndex];
        Residue* wiggleMe = residues[wiggleIndex];
        const DihedralAngleDataTable& metadataTable = dihedralAngleDataTable();
        this->superimpose(carbohydrateCoordinates, superimpositionTarget, superimposeMe);
        setAtomCoordinates(atoms, carbohydrateCoordinates);
        this->getCarbohydrate().Generate3DStructureFiles("./", "superimposed", {});
        this->determineWiggleLinkages(metadataTable, superimposeMe, wiggleMe);
        std::vector<Atom*> substrateWithoutSuperimpositionAtoms =
            findElementsNotInVector(getAtoms(this->getSubstrate().getAssemblies()), superimpositionTarget->getAtoms());
        std::vector<Atom*> substrateAtomsToAvoidOverlappingWith =
            findElementsNotInVector(substrateWithoutSuperimpositionAtoms, wigglingTarget->getAtoms());
        this->atomsToAvoid_ = substrateAtomsToAvoidOverlappingWith;
        const util::SparseVector<double>& elementRadii = vanDerWaalsRadii();
        this->setCurrentOverlap(CountOverlappingAtoms(
            elementRadii, constants::overlapTolerance, atomsToAvoid_, this->getCarbohydrate().getAtoms()));
        this->wiggleMeAtoms_ = {wiggleMe->FindAtom("C1"), wiggleMe->FindAtom("C3"), wiggleMe->FindAtom("C5")};
        this->wiggleTargetAtoms_ = {
            wigglingTarget->FindAtom("C1"), wigglingTarget->FindAtom("C3"), wigglingTarget->FindAtom("C5")};
        Atom* nullAtom = nullptr;
        if (util::contains(wiggleMeAtoms_, nullAtom) || util::contains(wiggleTargetAtoms_, nullAtom))
        {
            throw std::runtime_error("Did not find the cooordinates of the atoms required for wiggling\n");
        }
        this->setCurrentDistance(this->calculateDistance());
        int structureCount = this->minimizeDistance(
            elementRadii, metadataTable, inputStruct.persistCycles_, !inputStruct.isDeterministic_);
        this->minimizeDistance(elementRadii, metadataTable, inputStruct.persistCycles_, false, structureCount);
        this->getCarbohydrate().Generate3DStructureFiles("./", "finished", {});
    }

    int WiggleToSite::minimizeDistance(
        const util::SparseVector<double>& elementRadii,
        const DihedralAngleDataTable& metadataTable,
        int persistCycles,
        bool useMonteCarlo,
        int structureCount)
    {
        uint64_t seed = util::generateRandomSeed();
        pcg32 rng(seed);
        auto randomMetadata = [&rng](const DihedralAngleDataTable& table, const std::vector<size_t>& indices)
        { return util::weightedRandomOrder(rng, util::indicesToValues(table.weights, indices)); };
        double angleStandardDeviation = 2.0;
        auto randomAngle = [&rng, &angleStandardDeviation](DihedralAngleData metadata)
        {
            auto random = [&rng, &metadata](double stdCutoff, double lower, double upper)
            {
                double num = util::normalDistributionRandomDoubleWithCutoff(rng, -stdCutoff, stdCutoff);
                return metadata.default_angle + num * (num < 0 ? lower : upper);
            };
            std::function<double(const AngleLimit&)> onLimit = [&](const AngleLimit& dev)
            { return random(1.0, dev.lowerDeviationLimit, dev.upperDeviationLimit); };
            std::function<double(const AngleStd&)> onStd = [&](const AngleStd& dev)
            { return random(angleStandardDeviation, dev.lowerDeviationStd, dev.upperDeviationStd); };
            return onAngleDeviation(onLimit, onStd, metadata.angle_deviation);
        };

        int cycle = 0;
        while (cycle < persistCycles)
        {
            ++cycle;
            std::stringstream ss;
            ss << "Cycle " << cycle << "/" << persistCycles << "\n";
            util::log(__LINE__, __FILE__, util::INF, ss.str());
            for (auto& linkage : util::shuffleVector(rng, this->getWiggleLinkages()))
            {
                auto recordedShape = currentShape(metadataTable, linkage.rotatableDihedrals, linkage.dihedralMetadata);
                setShapeToPreference(
                    linkage,
                    linkageShapePreference(
                        randomMetadata, randomAngle, metadataTable, linkage.rotamerType, linkage.dihedralMetadata));
                double acceptance = util::uniformRandomDoubleWithinRange(rng, 0, 1);
                if (this->acceptDistance(useMonteCarlo, acceptance) && this->acceptOverlaps(elementRadii))
                {
                    cycle = 0; // reset when it improves
                }
                else
                {
                    setShape(linkage.rotatableDihedrals, recordedShape);
                }
            }
        }
        return structureCount;
    }

    //////////////////////////////////////////////////////////
    //                  PRIVATE FUNCTIONS                   //
    //////////////////////////////////////////////////////////
    void WiggleToSite::superimpose(
        std::vector<Coordinate>& carbohydrateCoordinates, const Residue* superimpositionTarget, Residue* superimposeMe)
    {
        // Limiting the selection to just these atoms as sometimes hydrogens or an oxygen is missing from xtal. That's
        // ok.
        std::vector<Atom*> superimposeMeAtoms = {
            superimposeMe->FindAtom("C1"), superimposeMe->FindAtom("C3"), superimposeMe->FindAtom("C5")};
        std::vector<Atom*> superTargetAtoms = {
            superimpositionTarget->FindAtom("C1"),
            superimpositionTarget->FindAtom("C3"),
            superimpositionTarget->FindAtom("C5")};
        std::vector<Coordinate> superimposeMeCoordinates = atomCoordinates(superimposeMeAtoms);
        std::vector<Coordinate> superTargetCoordinates = atomCoordinates(superTargetAtoms);
        Superimpose(
            superimposeMeCoordinates,
            superTargetCoordinates,
            carbohydrateCoordinates); // "alsoMoving" are the carbohydrate Coordinates
        setAtomCoordinates(superimposeMeAtoms, superimposeMeCoordinates);
        setAtomCoordinates(superTargetAtoms, superTargetCoordinates);
    }

    std::vector<ResidueLinkage>& WiggleToSite::determineWiggleLinkages(
        const DihedralAngleDataTable& metadataTable, Residue* startResidue, Residue* endResidue)
    {
        std::vector<Residue*> residuesInPath;
        bool targetFound = false;
        std::vector<Residue*> visitedResidues;
        std::cout << "Gonna find path between " << startResidue->getName() << " and " << endResidue->getName() << "\n";
        findPathBetweenElementsInGraph(startResidue, endResidue, visitedResidues, residuesInPath, targetFound);
        Residue* previousResidue = nullptr; // wanna skip the first iteration
        for (auto& residue : residuesInPath)
        {
            std::cout << residue->getName() << "_" << residue->getNumber() << ", ";
            if (previousResidue != nullptr)
            {
                ResidueLink link = findResidueLink({previousResidue, residue});
                wiggleLinkages_.emplace_back(createResidueLinkage(metadataTable, link));
            }
            previousResidue = residue;
        }
        std::cout << "\nLinkages I behold:\n" << std::endl;
        for (auto& linkage : this->getWiggleLinkages())
        {
            std::cout << linkage.name << ": "
                      << numberOfShapes(metadataTable, linkage.rotamerType, linkage.dihedralMetadata) << std::endl;
        }
        return this->getWiggleLinkages();
    }

    double WiggleToSite::calculateDistance()
    {
        return distance(wiggleTargetAtoms_[0]->coordinate(), wiggleMeAtoms_[0]->coordinate());
    }

    bool WiggleToSite::acceptOverlaps(const util::SparseVector<double>& elementRadii)
    {
        double overlapCount = CountOverlappingAtoms(
            elementRadii, constants::overlapTolerance, atomsToAvoid_, getCarbohydrate().getAtoms());
        if (compareOverlaps(overlapCount, this->getCurrentOverlap()) > 0)
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
            (useMonteCarlo && util::accept_via_metropolis_criterion(acceptance, distance - this->getCurrentDistance())))
        {
            this->setCurrentDistance(distance);
            return true;
        }
        return false;
    }
} // namespace gmml
