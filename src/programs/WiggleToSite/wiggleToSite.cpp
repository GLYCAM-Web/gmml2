#include "include/programs/WiggleToSite/wiggleToSite.hpp"

#include "include/CentralDataStructure/Selections/templatedSelections.hpp"
#include "include/CentralDataStructure/cdsFunctions.hpp"
#include "include/CentralDataStructure/residueLinkage/residueLinkageCreation.hpp"
#include "include/CentralDataStructure/residueLinkage/residueLinkageFunctions.hpp"
#include "include/assembly/assemblyBounds.hpp"
#include "include/assembly/assemblyOverlap.hpp"
#include "include/carbohydrate/dihedralAngleSearch.hpp"
#include "include/carbohydrate/dihedralShape.hpp"
#include "include/external/pcg/pcg_random.h"
#include "include/fileType/pdb/pdbFile.hpp"
#include "include/fileType/pdb/pdbModel.hpp"
#include "include/geometry/matrix.hpp"
#include "include/geometry/overlap.hpp"
#include "include/metadata/dihedralangledata.hpp"
#include "include/metadata/elements.hpp"
#include "include/preprocess/parameterManager.hpp"
#include "include/sequence/sequenceParser.hpp"
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
            const std::vector<Atom*>& atomsA,
            const std::vector<Atom*>& atomsB)
        {
            std::vector<Sphere> coordsA =
                assembly::toAtomBounds(elementRadii, atomElements(atomsA), atomCoordinates(atomsA));
            std::vector<Sphere> coordsB =
                assembly::toAtomBounds(elementRadii, atomElements(atomsB), atomCoordinates(atomsB));
            std::vector<Element> elementsA = atomElements(atomsA);
            std::vector<Element> elementsB = atomElements(atomsB);
            const PotentialTable potentialTable =
                gmml::potentialTable(elementRadii, util::vectorOr(foundElements(elementsA), foundElements(elementsB)));

            double overlapTolerance = 0.0;
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
    WiggleToSite::WiggleToSite(
        const preprocess::ParameterManager& parameterManager,
        const util::SparseVector<double>& elementRadii,
        const DihedralAngleDataTable& dihedralAngleData,
        const WiggleToSiteInputs& inputStruct)
        : substrate_(pdb::toPdbFile(inputStruct.substrateFile_, {pdb::InputType::modelsAsMolecules, false}))
    {
        initializeCarbohydrate(
            carbohydrate_,
            glycosidicLinkages_,
            parameterManager,
            elementRadii,
            dihedralAngleData,
            sequence::parseAndReorder(inputStruct.carbohydrateSequence_));
        generate3DStructureFiles(structured({&carbohydrate_}), carbohydrate_.getName(), "./", "initial", {});
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
            ss << "Selection for superimposition target: " << toString(inputStruct.superimpositionTargetResidue_)
               << "\nOr selection for wiggling target: " << toString(inputStruct.wigglingTargetResidue_)
               << " was not found\n";
            throw std::runtime_error(ss.str());
        }
        auto atoms = carbohydrate_.mutableAtoms();
        std::vector<Coordinate> carbohydrateCoordinates = atomCoordinates(atoms);
        std::vector<Residue*> residues = carbohydrate_.getResidues();
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
        this->superimpose(carbohydrateCoordinates, superimpositionTarget, superimposeMe);
        setAtomCoordinates(atoms, carbohydrateCoordinates);
        generate3DStructureFiles(structured({&carbohydrate_}), carbohydrate_.getName(), "./", "superimposed", {});
        this->determineWiggleLinkages(dihedralAngleData, superimposeMe, wiggleMe);
        std::vector<Atom*> substrateWithoutSuperimpositionAtoms =
            findElementsNotInVector(getAtoms(getAssemblies(getSubstrate())), superimpositionTarget->getAtoms());
        std::vector<Atom*> substrateAtomsToAvoidOverlappingWith =
            findElementsNotInVector(substrateWithoutSuperimpositionAtoms, wigglingTarget->getAtoms());
        this->atomsToAvoid_ = substrateAtomsToAvoidOverlappingWith;
        this->setCurrentOverlap(CountOverlappingAtoms(elementRadii, atomsToAvoid_, carbohydrate_.getAtoms()));
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
            elementRadii, dihedralAngleData, inputStruct.persistCycles_, !inputStruct.isDeterministic_);
        this->minimizeDistance(elementRadii, dihedralAngleData, inputStruct.persistCycles_, false, structureCount);
        generate3DStructureFiles(structured({&carbohydrate_}), carbohydrate_.getName(), "./", "finished", {});
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
                auto recordedShape = currentShape(metadataTable, linkage.rotatableBonds, linkage.dihedralMetadata);
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
                    setShape(linkage.rotatableBonds, recordedShape);
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
            superimposeMe->FindAtom("C5"), superimposeMe->FindAtom("C3"), superimposeMe->FindAtom("C1")};
        std::vector<Atom*> superTargetAtoms = {
            superimpositionTarget->FindAtom("C5"),
            superimpositionTarget->FindAtom("C3"),
            superimpositionTarget->FindAtom("C1")};
        auto toArray = [](const std::vector<Atom*>& atoms)
        {
            return std::array<Coordinate, 3> {
                {atoms[0]->coordinate(), atoms[1]->coordinate(), atoms[2]->coordinate()}
            };
        };
        std::array<Coordinate, 3> superimpose = toArray(superimposeMeAtoms);
        std::array<Coordinate, 3> target = toArray(superTargetAtoms);
        Matrix4x4 mat = superimposition(target, superimpose);
        setAtomCoordinates(superimposeMeAtoms, transform(mat, atomCoordinates(superimposeMeAtoms)));
        carbohydrateCoordinates = transform(mat, carbohydrateCoordinates);
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
        double overlapCount = CountOverlappingAtoms(elementRadii, atomsToAvoid_, carbohydrate_.getAtoms());
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
