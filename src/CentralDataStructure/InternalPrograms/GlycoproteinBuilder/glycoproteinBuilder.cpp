#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycoproteinBuilder.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/gpInputStructs.hpp"
#include "includes/CentralDataStructure/cdsFunctions/cdsFunctions.hpp"
#include "includes/CentralDataStructure/cdsFunctions/bondByDistance.hpp"
#include "includes/CentralDataStructure/cdsFunctions/atomicConnectivity.hpp"
#include "includes/CentralDataStructure/cdsFunctions/atomCoordinates.hpp"
#include "includes/CentralDataStructure/Writers/offWriter.hpp"
#include "includes/CentralDataStructure/Writers/pdbWriter.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbFile.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbResidue.hpp"
#include "includes/CentralDataStructure/Selections/residueSelections.hpp"
#include "includes/CentralDataStructure/Shapers/residueLinkage.hpp"
#include "includes/CentralDataStructure/Shapers/residueLinkageCreation.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralShape.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralAngleSearch.hpp"
#include "includes/CentralDataStructure/Overlaps/atomOverlaps.hpp"
#include "includes/CodeUtils/casting.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/files.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/random.hpp"
#include "includes/CodeUtils/strings.hpp"
#include "includes/External_Libraries/PCG/pcg_random.h"
#include "includes/External_Libraries/PCG/pcg_extras.h"

#include <variant>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <stdexcept>

struct GlycositeData
{
    std::vector<std::vector<cds::ResidueLinkage>> glycosidicLinkages;
    std::vector<std::vector<Residue*>> glycositeResidues;
    std::vector<cds::ResiduesWithOverlapWeight> glycositeResiduesWithWeights;
};

namespace
{
    GlycositeData toGlycositeData(std::vector<GlycosylationSite>& glycosites)
    {
        std::vector<std::vector<cds::ResidueLinkage>> glycosidicLinkages;
        std::vector<std::vector<Residue*>> glycositeResidues;
        for (auto& a : glycosites)
        {
            auto glycan = a.GetGlycan();
            glycosidicLinkages.push_back(glycan->GetGlycosidicLinkages());
            glycositeResidues.push_back(glycan->getResidues());
        }

        std::vector<cds::ResiduesWithOverlapWeight> glycositeResiduesWithWeights;
        for (auto& a : glycositeResidues)
        {
            glycositeResiduesWithWeights.push_back({a, std::vector<double>(a.size(), 1.0)});
        }
        return GlycositeData {glycosidicLinkages, glycositeResidues, glycositeResiduesWithWeights};
    }

    std::vector<OverlapResidues> toOverlapResidues(const std::vector<cds::Residue*>& proteinResidues,
                                                   std::vector<GlycosylationSite>& glycosites,
                                                   const std::vector<std::vector<Residue*>>& glycositeResidues)
    {
        std::vector<OverlapResidues> result;
        for (size_t n = 0; n < glycosites.size(); n++)
        {
            std::vector<Residue*> protein = proteinResidues;
            auto glycosite                = glycosites[n].GetResidue();
            protein.erase(std::remove(protein.begin(), protein.end(), glycosite), protein.end());
            protein.insert(protein.begin(), glycosite);

            std::vector<Residue*> glycan;
            for (auto& otherSite : codeUtils::withoutNth(n, glycositeResidues))
            {
                codeUtils::insertInto(glycan, otherSite);
            }

            result.push_back({protein, glycan});
        }
        return result;
    }

    Residue* selectResidueFromInput(cds::Assembly* glycoprotein, const std::string userSelection)
    { // Chain_residueNumber_insertionCode* *optional.
        std::vector<std::string> splitUserSelection = codeUtils::split(userSelection, '_');
        if (splitUserSelection.size() < 2)
        {
            throw std::runtime_error(
                "userSelection (" + userSelection +
                ") for residue to glycosylate is incorrect format.\nMust be "
                "chain_residueNumber_insertionCode.\nInsertionCode is optional. Chain can be ? if no "
                "chain numbers are in input.\nExamples: ?_24_? or ?_24 will use the first residue it "
                "encounters numbered 24. A_24_B is A chain, residue 24, insertion code B");
        }
        std::string userSelectedChain = "";
        if (splitUserSelection.at(0) != "?")
        {
            userSelectedChain = splitUserSelection.at(0);
        }
        std::string userSelectedResidue = splitUserSelection.at(1);
        if (splitUserSelection.size() == 3)
        {
            userSelectedResidue += splitUserSelection.at(2); // So will be 24A or just 24.
        }
        gmml::log(__LINE__, __FILE__, gmml::INF, "We working with " + userSelectedChain + "_" + userSelectedResidue);
        for (auto& residue : glycoprotein->getResidues())
        {
            pdb::PdbResidue* pdbResidue = codeUtils::erratic_cast<pdb::PdbResidue*>(residue);
            // std::cout << pdbResidue->getChainId() << "_";
            //        std::cout << pdbResidue->getNumberAndInsertionCode() << "\n";
            if ((pdbResidue->getChainId() == userSelectedChain) &&
                (pdbResidue->getNumberAndInsertionCode() == userSelectedResidue))
            {
                gmml::log(__LINE__, __FILE__, gmml::INF, "Id of selected glycosite: " + pdbResidue->printId());
                return residue;
            }
        }
        return nullptr;
    }

    std::vector<GlycosylationSite> createGlycosites(cds::Assembly* glycoprotein,
                                                    std::vector<glycoprotein::GlycositeInput> glycositesInputVector)
    {
        std::vector<GlycosylationSite> glycosites;
        for (auto& glycositeInput : glycositesInputVector)
        {
            gmml::log(__LINE__, __FILE__, gmml::INF,
                      "Creating glycosite on residue " + glycositeInput.proteinResidueId + " with glycan " +
                          glycositeInput.glycanInput);
            Carbohydrate* carb = codeUtils::erratic_cast<Carbohydrate*>(
                glycoprotein->addMolecule(std::make_unique<Carbohydrate>(glycositeInput.glycanInput)));
            Residue* glycositeResidue = selectResidueFromInput(glycoprotein, glycositeInput.proteinResidueId);
            if (glycositeResidue == nullptr)
            {
                throw std::runtime_error("Did not find a residue with id matching " + glycositeInput.proteinResidueId);
            }
            gmml::log(__LINE__, __FILE__, gmml::INF,
                      "About to emplace_back to glycosites with: " + glycositeInput.proteinResidueId + " and glycan " +
                          glycositeInput.glycanInput);
            unsigned int highestResidueNumber = cdsSelections::findHighestResidueNumber(glycoprotein->getResidues());
            glycosites.emplace_back(glycositeResidue, carb, highestResidueNumber);
            for (auto& linkage : carb->GetGlycosidicLinkages())
            {
                cds::determineResiduesForOverlapCheck(linkage); // Now that the protein residue is attached.
            }
            //	    std::cout << "Done with glycan" << std::endl;
            gmml::log(__LINE__, __FILE__, gmml::INF,
                      "Completed creating glycosite on residue " + glycositeInput.proteinResidueId + " with glycan " +
                          glycositeInput.glycanInput);
        }
        return glycosites;
    }
} // namespace

GlycoproteinBuilder::GlycoproteinBuilder(glycoprotein::GlycoproteinBuilderInputs inputStruct,
                                         pdb::PreprocessorOptions preprocessingOptions)
    : settings(inputStruct)
{
    try
    {
        pdbFile = pdb::PdbFile(inputStruct.substrateFileName);
        if (!inputStruct.skipMDPrep)
        {
            gmml::log(__LINE__, __FILE__, gmml::INF, "Performing MDPrep aka preprocessing.");
            pdbFile.PreProcess(preprocessingOptions);
        }
        glycoprotein_                                = &pdbFile.mutableAssemblies().front();
        std::vector<cds::Residue*> gpInitialResidues = glycoprotein_->getResidues();
        cds::setIntraConnectivity(gpInitialResidues);
        gmml::log(__LINE__, __FILE__, gmml::INF, "Attaching Glycans To Glycosites.");
        proteinResidues_ = glycoprotein_->getResidues();
        glycosites_      = createGlycosites(glycoprotein_, inputStruct.glycositesInputVector);
        cds::setInterConnectivity(gpInitialResidues); // do the inter here, so that the whole protein isn't included as
                                                      // overlap residues in the glycan linkages.
    }
    catch (const std::string& errorMessage)
    {
        gmml::log(__LINE__, __FILE__, gmml::ERR, errorMessage);
        throw std::runtime_error(errorMessage);
    }
    gmml::log(__LINE__, __FILE__, gmml::INF, "Initialization of Glycoprotein builder complete!");
    return;
}

void GlycoproteinBuilder::ResolveOverlaps(std::string outputDir)
{
    double selfWeight    = 1000000.0;
    double proteinWeight = 1000.0;
    double glycanWeight  = 1.0;
    auto overlapWeight   = OverlapWeight {proteinWeight, glycanWeight, selfWeight};

    uint64_t seed = settings.isDeterministic ? settings.seed : codeUtils::generateRandomSeed();
    pcg32 rng(seed);

    auto randomMetadata = [&rng](GlycamMetadata::DihedralAngleDataVector metadataVector)
    {
        auto weights = GlycamMetadata::dihedralAngleDataWeights(metadataVector);
        return codeUtils::weightedRandomOrder(rng, weights);
    };
    bool freezeGlycositeResidueConformation = settings.freezeGlycositeResidueConformation;
    double preferenceDeviation              = 2.0;
    double searchDeviation                  = 0.5;
    double searchIncrement                  = 1.0;
    auto standardDeviation = [&preferenceDeviation, &searchDeviation](const GlycamMetadata::DihedralAngleData& metadata)
    {
        auto angle = metadata.angle_deviation;
        if (std::holds_alternative<GlycamMetadata::AngleLimit>(angle))
        {
            auto dev         = std::get<GlycamMetadata::AngleLimit>(angle);
            double max_std   = preferenceDeviation + searchDeviation;
            double lower_std = dev.lowerDeviationLimit / max_std;
            double upper_std = dev.upperDeviationLimit / max_std;
            return std::pair<double, double> {lower_std, upper_std};
        }
        else if (std::holds_alternative<GlycamMetadata::AngleStd>(angle))
        {
            auto dev = std::get<GlycamMetadata::AngleStd>(angle);
            return std::pair<double, double> {dev.lowerDeviationStd, dev.upperDeviationStd};
        }
        else
        {
            throw std::runtime_error("unknown angle deviation type");
        }
    };
    auto searchAngles = [&standardDeviation, &searchIncrement](const GlycamMetadata::DihedralAngleData& metadata,
                                                               double preference, double deviation)
    {
        auto std = standardDeviation(metadata);
        return cds::evenlySpacedAngles(preference, deviation * std.first, deviation * std.second, searchIncrement);
    };
    auto searchSettings = cds::AngleSearchSettings {searchDeviation, searchAngles};
    auto randomAngle    = [&rng, &standardDeviation, &preferenceDeviation](GlycamMetadata::DihedralAngleData metadata)
    {
        double stdCutoff = preferenceDeviation;
        double num       = codeUtils::normalDistributionRandomDoubleWithCutoff(rng, -stdCutoff, stdCutoff);
        auto std         = standardDeviation(metadata);
        return metadata.default_angle + num * (num < 0 ? std.first : std.second);
    };
    auto glycositeResidueConformationFrozen =
        [](std::vector<cds::ResidueLinkageShapePreference> preference, const std::vector<cds::ResidueLinkage>& linkages)
    {
        if (std::holds_alternative<cds::ConformerShapePreference>(preference[0]))
        {
            auto& firstLinkage = linkages[0];
            auto& pref         = std::get<cds::ConformerShapePreference>(preference[0]);
            auto shape         = cds::currentShape(firstLinkage.rotatableDihedrals, firstLinkage.dihedralMetadata);
            for (size_t n = 0; n < firstLinkage.rotatableDihedrals.size(); n++)
            {
                auto& name = firstLinkage.dihedralMetadata[n][0].dihedral_angle_name_;
                if ((name == "Chi1") || (name == "Chi2"))
                {
                    pref.isFrozen[n]         = true;
                    size_t metadata          = shape[n].metadataIndex;
                    pref.angles[n][metadata] = shape[n].value;
                    pref.metadataOrder       = {metadata};
                }
            }
        }
        return preference;
    };
    auto randomizeShape = [&randomMetadata, &randomAngle, &freezeGlycositeResidueConformation,
                           &glycositeResidueConformationFrozen](const std::vector<cds::ResidueLinkage>& linkages)
    {
        auto preference = cds::linkageShapePreference(randomMetadata, randomAngle, linkages);
        return freezeGlycositeResidueConformation ? glycositeResidueConformationFrozen(preference, linkages)
                                                  : preference;
    };
    auto wiggleGlycan = [&searchSettings](OverlapWeight weight, std::vector<cds::ResidueLinkage>& linkages,
                                          const std::vector<cds::ResidueLinkageShapePreference>& preferences,
                                          const OverlapResidues& overlapResidues)
    {
        return wiggleGlycosite(searchSettings, weight, linkages, preferences, overlapResidues);
    };

    auto resolveOverlapsWithWiggler =
        [&](const std::vector<cds::Residue*>& proteinResidues, std::vector<GlycosylationSite>& glycosites)
    {
        auto data                = toGlycositeData(glycosites);
        auto& glycosidicLinkages = data.glycosidicLinkages;
        auto& weightedResidues   = data.glycositeResiduesWithWeights;
        auto overlapResidues     = toOverlapResidues(proteinResidues, glycosites, data.glycositeResidues);

        std::vector<std::vector<cds::ResidueLinkageShapePreference>> glycositePreferences;
        std::vector<std::vector<std::vector<cds::AngleWithMetadata>>> glycositeShape;
        for (auto& linkages : glycosidicLinkages)
        {
            auto preference = randomizeShape(linkages);
            cds::setShapeToPreference(linkages, preference);
            glycositePreferences.push_back(preference);
            glycositeShape.push_back(cds::currentShape(linkages));
        }
        for (size_t site : codeUtils::shuffleVector(rng, codeUtils::indexVector(glycosidicLinkages)))
        {
            glycositeShape[site] = wiggleGlycan(overlapWeight, glycosidicLinkages[site], glycositePreferences[site],
                                                overlapResidues[site]);
        }
        cds::Overlap initialOverlap =
            totalOverlaps(overlapWeight, overlapResidues, weightedResidues, glycosidicLinkages);
        auto overlapSites = determineSitesWithOverlap(codeUtils::indexVector(glycosidicLinkages), glycosidicLinkages,
                                                      overlapResidues, weightedResidues);
        auto initialState = GlycoproteinState {initialOverlap, overlapSites, glycositePreferences, glycositeShape};
        GlycoproteinState currentState =
            randomDescent(rng, randomizeShape, wiggleGlycan, settings.persistCycles, overlapWeight, glycosidicLinkages,
                          initialState, overlapResidues, weightedResidues);
        gmml::log(__LINE__, __FILE__, gmml::INF, "Overlap: " + std::to_string(currentState.overlap.count));
    };

    auto writeOffFile = [](const std::vector<cds::Residue*>& residues, const std::string& prefix)
    {
        std::string fileName = prefix + ".off";
        std::ofstream outFileStream;
        outFileStream.open(fileName.c_str());
        cds::WriteToOffFile(residues, outFileStream, "GLYCOPROTEINBUILDER");
        outFileStream.close();
    };

    auto writePdbFile = [](const std::vector<std::vector<size_t>>& residueIndices,
                           const std::vector<std::vector<bool>>& residueTER, const cds::PdbWriterData& data,
                           const std::vector<std::pair<int, int>>& connectionNumbers, const std::string& prefix)
    {
        std::string fileName = prefix + ".pdb";
        std::ofstream outFileStream;
        outFileStream.open(fileName.c_str());
        for (size_t n = 0; n < residueIndices.size(); n++)
        {
            cds::writeMoleculeToPdb(outFileStream, residueIndices[n], residueTER[n], data);
        }
        cds::writeConectCards(outFileStream, connectionNumbers);
        outFileStream.close();
    };

    auto printDihedralAnglesAndOverlapOfGlycosites =
        [](const std::vector<cds::Residue*>& proteinResidues, std::vector<GlycosylationSite>& glycosites)
    {
        auto data               = toGlycositeData(glycosites);
        auto overlapResiduesVec = toOverlapResidues(proteinResidues, glycosites, data.glycositeResidues);
        std::stringstream logss;
        for (size_t n = 0; n < glycosites.size(); n++)
        {
            auto& glycosite       = glycosites[n];
            auto glycan           = glycosite.GetGlycan();
            auto glycanResidues   = glycan->getResidues();
            auto& overlapResidues = overlapResiduesVec[n];
            auto glycanResiduesWithWeight =
                cds::ResiduesWithOverlapWeight {glycanResidues, std::vector<double>(glycanResidues.size(), 1.0)};
            auto proteinOverlaps = countOverlaps(overlapResidues.protein, glycanResiduesWithWeight);
            auto glycanOverlaps  = countOverlaps(overlapResidues.glycan, glycanResiduesWithWeight);
            auto selfOverlaps    = intraGlycanOverlaps(glycan->GetGlycosidicLinkages());
            logss << "Residue ID: " << glycosite.GetResidue()->getStringId()
                  << ", protein overlap: " << proteinOverlaps.count << ", glycan overlap: " << glycanOverlaps.count
                  << ", self overlap: " << selfOverlaps.count;
            gmml::log(__LINE__, __FILE__, gmml::INF, logss.str());
        }
    };

    Assembly* assembly = getGlycoprotein();
    std::vector<cds::Residue*> residues;
    std::vector<cds::Atom*> atoms;
    std::vector<std::vector<bool>> residueTER;
    std::vector<std::vector<size_t>> residueIndices;
    std::vector<std::vector<size_t>> atomIndices;

    for (auto& molecule : assembly->getMolecules())
    {
        std::vector<cds::Residue*> moleculeResidues = molecule->getResidues();
        residueIndices.push_back(codeUtils::indexVectorWithOffset(residues.size(), moleculeResidues));
        codeUtils::insertInto(residues, moleculeResidues);
        residueTER.push_back(cds::residueTER(cds::residueTypes(moleculeResidues)));
        for (auto& residue : moleculeResidues)
        {
            std::vector<Atom*> residueAtoms = residue->getAtoms();
            atomIndices.push_back(codeUtils::indexVectorWithOffset(atoms.size(), residueAtoms));
            codeUtils::insertInto(atoms, residueAtoms);
        }
    }

    std::vector<std::string> recordNames(atoms.size(), "ATOM");
    std::vector<std::string> chainIds(residues.size(), "");
    std::vector<std::string> insertionCodes(residues.size(), "");
    cds::ResiduePdbData residuePdbData(atomIndices, residueNumbers(residues), residueNames(residues), chainIds,
                                       insertionCodes);
    cds::AtomPdbData atomPdbData(atoms, recordNames);
    cds::PdbWriterData writerData {residuePdbData, atomPdbData};

    auto pdbResidues =
        cdsSelections::selectResiduesByType(residues, {cds::ResidueType::Sugar, cds::ResidueType::Derivative,
                                                       cds::ResidueType::Aglycone, cds::ResidueType::Undefined});

    std::vector<std::pair<int, int>> noConnections = {};
    writePdbFile(residueIndices, residueTER, writerData, noConnections, outputDir + "glycoprotein_initial");
    resolveOverlapsWithWiggler(proteinResidues_, glycosites_);
    printDihedralAnglesAndOverlapOfGlycosites(proteinResidues_, glycosites_);
    writerData.atoms.coordinates = cds::atomCoordinates(atoms);
    writePdbFile(residueIndices, residueTER, writerData, noConnections, outputDir + "glycoprotein");
    cds::serializeNumbers(atoms);
    cds::serializeNumbers(residues);
    writerData.atoms.numbers    = cds::atomNumbers(atoms);
    writerData.residues.numbers = cds::residueNumbers(residues);
    std::vector<std::pair<int, int>> connectionNumbers =
        atomPairNumbers(cds::atomPairsConnectedToOtherResidues(pdbResidues));
    writeOffFile(residues, outputDir + "glycoprotein");
    writePdbFile(residueIndices, residueTER, writerData, connectionNumbers, outputDir + "glycoprotein_serialized");

    for (size_t count = 0; count < settings.number3DStructures; count++)
    {
        printDihedralAnglesAndOverlapOfGlycosites(proteinResidues_, glycosites_);
        std::stringstream prefix;
        prefix << count << "_glycoprotein";
        resolveOverlapsWithWiggler(proteinResidues_, glycosites_);
        writerData.atoms.coordinates = cds::atomCoordinates(atoms);
        writePdbFile(residueIndices, residueTER, writerData, connectionNumbers, outputDir + prefix.str());
    }
}
