#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycoproteinBuilder.hpp"
#include "includes/CentralDataStructure/cdsFunctions/bondByDistance.hpp"
#include "includes/CentralDataStructure/cdsFunctions/atomicConnectivity.hpp"
#include "includes/CentralDataStructure/Writers/offWriter.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/gpInputStructs.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/files.hpp"
#include "includes/CodeUtils/strings.hpp" // split
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/random.hpp"
#include "includes/CentralDataStructure/Writers/pdbWriter.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbFile.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbResidue.hpp"
#include "includes/CentralDataStructure/Selections/residueSelections.hpp" // selectResiduesByType
#include "includes/CentralDataStructure/Shapers/residueLinkage.hpp"
#include "includes/CentralDataStructure/Shapers/residueLinkageCreation.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralShape.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralAngleSearch.hpp"
#include "includes/CentralDataStructure/Overlaps/atomOverlaps.hpp"
#include "includes/External_Libraries/PCG/pcg_random.h"
#include "includes/External_Libraries/PCG/pcg_extras.h"

#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <stdexcept>

using cds::Assembly;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
GlycoproteinBuilder::GlycoproteinBuilder(glycoprotein::GlycoproteinBuilderInputs inputStruct,
                                         pdb::PreprocessorOptions preprocessingOptions)
    : settings(inputStruct)
{
    try
    {
        pdb::PdbFile pdbFile(inputStruct.substrateFileName);
        if (!inputStruct.skipMDPrep)
        {
            gmml::log(__LINE__, __FILE__, gmml::INF, "Performing MDPrep aka preprocessing.");
            pdbFile.PreProcess(preprocessingOptions);
        }
        glycoprotein_                                = std::move(*(pdbFile.getAssemblies().front()));
        std::vector<cds::Residue*> gpInitialResidues = glycoprotein_.getResidues();
        cds::setIntraConnectivity(gpInitialResidues);
        gmml::log(__LINE__, __FILE__, gmml::INF, "Attaching Glycans To Glycosites.");
        proteinResidues_ = getGlycoprotein()->getResidues();
        this->CreateGlycosites(inputStruct.glycositesInputVector);
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

//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//                       FUNCTIONS                      //
//////////////////////////////////////////////////////////

void GlycoproteinBuilder::WriteOffFile(const std::string prefix)
{
    std::string fileName = prefix + ".off";
    std::ofstream outFileStream;
    outFileStream.open(fileName.c_str());
    cds::WriteAssemblyToOffFile(this->getGlycoprotein(), outFileStream, "GLYCOPROTEINBUILDER");
    outFileStream.close();
    return;
}

void GlycoproteinBuilder::WritePdbFile(const std::string prefix, const bool writeConectSection)
{
    std::string fileName = prefix + ".pdb";
    std::ofstream outFileStream;
    outFileStream.open(fileName.c_str());
    cds::writeAssemblyToPdb(outFileStream, this->getGlycoprotein()->getMolecules());
    if (writeConectSection)
    {
        cds::writeConectCards(outFileStream, cdsSelections::selectResiduesByType(
                                                 this->getGlycoprotein()->getResidues(),
                                                 {cds::ResidueType::Sugar, cds::ResidueType::Derivative,
                                                  cds::ResidueType::Aglycone, cds::ResidueType::Undefined}));
    }
    outFileStream.close();
    return;
}

void GlycoproteinBuilder::ResolveOverlaps()
{ // First time here so get default output files that might be deterministic for tests.
    std::vector<std::vector<cds::ResidueLinkage>> glycosidicLinkages;
    std::vector<std::vector<Residue*>> glycositeResidues;
    for (auto& a : glycosites_)
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

    double selfWeight    = 1000000.0;
    double proteinWeight = 1000.0;
    double glycanWeight  = 1.0;
    auto overlapWeight   = OverlapWeight {proteinWeight, glycanWeight, selfWeight};

    overlapResidues_.clear();
    for (size_t n = 0; n < glycosites_.size(); n++)
    {
        std::vector<Residue*> protein = proteinResidues_;
        auto glycosite                = glycosites_[n].GetResidue();
        protein.erase(std::remove(protein.begin(), protein.end(), glycosite), protein.end());
        protein.insert(protein.begin(), glycosite);

        std::vector<Residue*> glycan;
        for (auto& otherSite : codeUtils::withoutNth(n, glycositeResidues))
        {
            codeUtils::insertInto(glycan, otherSite);
        }

        overlapResidues_.push_back({protein, glycan});
    }

    uint64_t seed = settings.isDeterministic ? 0 : codeUtils::generateRandomSeed();
    pcg32 rng(seed);

    auto randomMetadata = [&rng](GlycamMetadata::DihedralAngleDataVector metadataVector)
    {
        auto weights = GlycamMetadata::dihedralAngleDataWeights(metadataVector);
        return codeUtils::weightedRandomOrder(rng, weights);
    };
    double preferenceDeviation = 2.0;
    double searchDeviation     = 0.5;
    double searchIncrement     = 1.0;
    auto searchAngles =
        [&searchIncrement](const GlycamMetadata::DihedralAngleData& metadata, double preference, double deviation)
    {
        return cds::evenlySpacedAngles(preference, deviation, searchIncrement, metadata);
    };
    auto searchSettings = cds::AngleSearchSettings {searchDeviation, searchAngles};
    auto randomAngle    = [&rng, &preferenceDeviation](GlycamMetadata::DihedralAngleData metadata)
    {
        double stdCutoff = preferenceDeviation;
        double num       = codeUtils::normalDistributionRandomDoubleWithCutoff(rng, -stdCutoff, stdCutoff);
        return metadata.default_angle_value_ + num * (num < 0 ? metadata.lower_deviation_ : metadata.upper_deviation_);
    };
    auto randomizeShape = [&randomMetadata, &randomAngle](std::vector<cds::ResidueLinkage> linkages)
    {
        return cds::linkageShapePreference(randomMetadata, randomAngle, linkages);
    };
    auto wiggleGlycan = [&searchSettings](OverlapWeight weight, std::vector<cds::ResidueLinkage>& linkages,
                                          const std::vector<cds::ResidueLinkageShapePreference>& preferences,
                                          const OverlapResidues& overlapResidues)
    {
        return wiggleGlycosite(searchSettings, weight, linkages, preferences, overlapResidues);
    };

    auto resolveOverlapsWithWiggler = [&](std::vector<std::vector<cds::ResidueLinkage>>& glycosidicLinkages,
                                          const std::vector<OverlapResidues>& overlapResidues,
                                          const std::vector<cds::ResiduesWithOverlapWeight>& glycositeResidues)
    {
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
            totalOverlaps(overlapWeight, overlapResidues, glycositeResidues, glycosidicLinkages);
        auto initialState = GlycoproteinState {initialOverlap, glycositePreferences, glycositeShape};
        GlycoproteinState currentState =
            randomDescent(rng, randomizeShape, wiggleGlycan, settings.persistCycles, overlapWeight, glycosidicLinkages,
                          initialState, overlapResidues, glycositeResidues);
        gmml::log(__LINE__, __FILE__, gmml::INF, "Overlap: " + std::to_string(currentState.overlap.count));
    };

    this->WritePdbFile("glycoprotein_initial");
    resolveOverlapsWithWiggler(glycosidicLinkages, overlapResidues_, glycositeResiduesWithWeights);
    this->PrintDihedralAnglesAndOverlapOfGlycosites();
    this->WritePdbFile("glycoprotein");
    this->WriteOffFile("glycoprotein");
    this->WritePdbFile("glycoprotein_serialized");

    for (size_t count = 0; count < settings.number3DStructures; count++)
    {
        resolveOverlapsWithWiggler(glycosidicLinkages, overlapResidues_, glycositeResiduesWithWeights);
        this->PrintDihedralAnglesAndOverlapOfGlycosites();
        std::stringstream prefix;
        prefix << count << "_glycoprotein";
        this->WritePdbFile(prefix.str(), true);
    }
}

void GlycoproteinBuilder::CreateGlycosites(std::vector<glycoprotein::GlycositeInput> glycositesInputVector)
{
    for (auto& glycositeInput : glycositesInputVector)
    {
        gmml::log(__LINE__, __FILE__, gmml::INF,
                  "Creating glycosite on residue " + glycositeInput.proteinResidueId + " with glycan " +
                      glycositeInput.glycanInput);
        Carbohydrate* carb = static_cast<Carbohydrate*>(
            glycoprotein_.addMolecule(std::make_unique<Carbohydrate>(glycositeInput.glycanInput)));
        Residue* glycositeResidue = this->SelectResidueFromInput(glycositeInput.proteinResidueId);
        if (glycositeResidue == nullptr)
        {
            throw std::runtime_error("Did not find a residue with id matching " + glycositeInput.proteinResidueId);
        }
        gmml::log(__LINE__, __FILE__, gmml::INF,
                  "About to emplace_back to glycosites with: " + glycositeInput.proteinResidueId + " and glycan " +
                      glycositeInput.glycanInput);
        unsigned int highestResidueNumber =
            cdsSelections::findHighestResidueNumber(this->getGlycoprotein()->getResidues());
        glycosites_.emplace_back(glycositeResidue, carb, highestResidueNumber);
        for (auto& linkage : carb->GetGlycosidicLinkages())
        {
            cds::determineResiduesForOverlapCheck(linkage); // Now that the protein residue is attached.
        }
        //	    std::cout << "Done with glycan" << std::endl;
        gmml::log(__LINE__, __FILE__, gmml::INF,
                  "Completed creating glycosite on residue " + glycositeInput.proteinResidueId + " with glycan " +
                      glycositeInput.glycanInput);
    }
    return;
}

Residue* GlycoproteinBuilder::SelectResidueFromInput(const std::string userSelection)
{ // Chain_residueNumber_insertionCode* *optional.
    std::vector<std::string> splitUserSelection = codeUtils::split(userSelection, '_');
    if (splitUserSelection.size() < 2)
    {
        throw std::runtime_error("userSelection (" + userSelection +
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
    for (auto& residue : this->getGlycoprotein()->getResidues())
    {
        pdb::PdbResidue* pdbResidue = static_cast<pdb::PdbResidue*>(residue);
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

void GlycoproteinBuilder::PrintDihedralAnglesAndOverlapOfGlycosites()
{
    std::stringstream logss;
    for (size_t n = 0; n < glycosites_.size(); n++)
    {
        auto& glycosite       = glycosites_[n];
        auto glycan           = glycosite.GetGlycan();
        auto glycanResidues   = glycan->getResidues();
        auto& overlapResidues = overlapResidues_[n];
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
    return;
}
