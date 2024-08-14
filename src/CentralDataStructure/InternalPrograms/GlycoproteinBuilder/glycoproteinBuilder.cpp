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
{
    try
    {
        this->SetIsDeterministic(inputStruct.isDeterministic_);
        this->SetStructuresToOutput(inputStruct.number3DStructures_);
        this->SetPersistCycles(inputStruct.persistCycles_);
        this->SetOverlapTolerance(inputStruct.overlapTolerance_);
        pdb::PdbFile pdbFile(inputStruct.substrateFileName_);
        if (!inputStruct.skipMDPrep_)
        {
            gmml::log(__LINE__, __FILE__, gmml::INF, "Performing MDPrep aka preprocessing.");
            pdbFile.PreProcess(preprocessingOptions);
        }
        glycoprotein_                                = std::move(*(pdbFile.getAssemblies().front()));
        std::vector<cds::Residue*> gpInitialResidues = glycoprotein_.getResidues();
        cds::setIntraConnectivity(gpInitialResidues);
        gmml::log(__LINE__, __FILE__, gmml::INF, "Attaching Glycans To Glycosites.");
        proteinResidues_ = getGlycoprotein()->getResidues();
        this->CreateGlycosites(inputStruct.glycositesInputVector_);
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

    double proteinOverlapWeight = 5.0;
    double glycanOverlapWeight  = 1.0;

    overlapResidues_.clear();
    for (size_t n = 0; n < glycosites_.size(); n++)
    {
        std::vector<Residue*> protein = proteinResidues_;
        protein.erase(std::remove(protein.begin(), protein.end(), glycosites_[n].GetResidue()), protein.end());

        std::vector<Residue*> glycan;
        for (auto& otherSite : codeUtils::withoutNth(n, glycositeResidues))
        {
            codeUtils::insertInto(glycan, otherSite);
        }

        std::vector<double> proteinWeight(protein.size(), proteinOverlapWeight);
        std::vector<double> glycanWeight(glycan.size(), glycanOverlapWeight);

        overlapResidues_.push_back(
            {codeUtils::vectorAppend(protein, glycan), codeUtils::vectorAppend(proteinWeight, glycanWeight)});
    }

    int persistCycles    = GetPersistCycles();
    int overlapTolerance = GetOverlapTolerance();

    uint64_t seed = GetIsDeterministic() ? 0 : codeUtils::generateRandomSeed();
    pcg32 rng(seed);

    auto randomMetadata = [&rng](gmml::MolecularMetadata::GLYCAM::DihedralAngleDataVector metadataVector)
    {
        auto weights = gmml::MolecularMetadata::GLYCAM::dihedralAngleDataWeights(metadataVector);
        return codeUtils::weightedRandomOrder(rng, weights);
    };
    double angleStandardDeviation = 2.0;
    double angleIncrement         = 1.0;
    auto searchAngles =
        [&angleStandardDeviation, &angleIncrement](const gmml::MolecularMetadata::GLYCAM::DihedralAngleData& metadata)
    {
        return cds::evenlySpacedAngles(angleStandardDeviation, angleIncrement, metadata);
    };
    auto randomAngle = [&rng, &angleStandardDeviation](gmml::MolecularMetadata::GLYCAM::DihedralAngleData metadata)
    {
        double stdCutoff = angleStandardDeviation;
        double num       = codeUtils::normalDistributionRandomDoubleWithCutoff(rng, -stdCutoff, stdCutoff);
        return metadata.default_angle_value_ + num * (num < 0 ? metadata.lower_deviation_ : metadata.upper_deviation_);
    };
    auto randomizeShape = [&randomMetadata, &randomAngle](std::vector<cds::ResidueLinkage> linkages)
    {
        return cds::linkageShapePreference(randomMetadata, randomAngle, linkages);
    };

    auto resolveOverlapsWithWiggler = [&](std::vector<std::vector<cds::ResidueLinkage>>& glycosidicLinkages,
                                          const std::vector<cds::ResiduesWithOverlapWeight>& overlapResidues,
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
        bool useMonteCarlo          = true;
        bool wiggleFirstLinkageOnly = true;
        cds::Overlap initialOverlap = countTotalOverlaps(overlapResidues, glycositeResidues);
        auto initialState           = GlycoproteinState {initialOverlap, glycositePreferences, glycositeShape};
        GlycoproteinState currentState =
            wiggleSitesWithOverlaps(rng, searchAngles, overlapTolerance, persistCycles, wiggleFirstLinkageOnly,
                                    glycosidicLinkages, initialState, overlapResidues, glycositeResidues);
        gmml::log(__LINE__, __FILE__, gmml::INF, "1. Overlap: " + std::to_string(currentState.overlap.count));
        currentState = randomDescent(rng, randomizeShape, searchAngles, useMonteCarlo, persistCycles, overlapTolerance,
                                     glycosidicLinkages, currentState, overlapResidues, glycositeResidues);
        gmml::log(__LINE__, __FILE__, gmml::INF, "1r. Overlap: " + std::to_string(currentState.overlap.count));
        currentState = wiggleSitesWithOverlaps(rng, searchAngles, overlapTolerance, persistCycles, false,
                                               glycosidicLinkages, currentState, overlapResidues, glycositeResidues);
        gmml::log(__LINE__, __FILE__, gmml::INF, "2. Overlap: " + std::to_string(currentState.overlap.count));
        currentState =
            wiggleSitesWithOverlaps(rng, searchAngles, overlapTolerance, persistCycles, wiggleFirstLinkageOnly,
                                    glycosidicLinkages, currentState, overlapResidues, glycositeResidues);
        gmml::log(__LINE__, __FILE__, gmml::INF, "3. Overlap: " + std::to_string(currentState.overlap.count));
    };

    this->WritePdbFile("glycoprotein_initial");
    resolveOverlapsWithWiggler(glycosidicLinkages, overlapResidues_, glycositeResiduesWithWeights);
    this->PrintDihedralAnglesAndOverlapOfGlycosites();
    this->WritePdbFile("glycoprotein");
    this->WriteOffFile("glycoprotein");
    this->WritePdbFile("glycoprotein_serialized");

    while (this->GetNumberOfOuputtedStructures() < this->GetNumberOfStructuresToOutput())
    {
        resolveOverlapsWithWiggler(glycosidicLinkages, overlapResidues_, glycositeResiduesWithWeights);
        this->PrintDihedralAnglesAndOverlapOfGlycosites();
        std::stringstream prefix;
        prefix << this->GetNumberOfOuputtedStructures() << "_glycoprotein";
        this->WritePdbFile(prefix.str(), true);
        this->incrementOutputtedStructures();
    }
}

void GlycoproteinBuilder::CreateGlycosites(std::vector<glycoprotein::GlycositeInput> glycositesInputVector)
{
    for (auto& glycositeInput : glycositesInputVector)
    {
        gmml::log(__LINE__, __FILE__, gmml::INF,
                  "Creating glycosite on residue " + glycositeInput.proteinResidueId_ + " with glycan " +
                      glycositeInput.glycanInput_);
        Carbohydrate* carb = static_cast<Carbohydrate*>(
            glycoprotein_.addMolecule(std::make_unique<Carbohydrate>(glycositeInput.glycanInput_)));
        Residue* glycositeResidue = this->SelectResidueFromInput(glycositeInput.proteinResidueId_);
        if (glycositeResidue == nullptr)
        {
            throw std::runtime_error("Did not find a residue with id matching " + glycositeInput.proteinResidueId_);
        }
        gmml::log(__LINE__, __FILE__, gmml::INF,
                  "About to emplace_back to glycosites with: " + glycositeInput.proteinResidueId_ + " and glycan " +
                      glycositeInput.glycanInput_);
        unsigned int highestResidueNumber =
            cdsSelections::findHighestResidueNumber(this->getGlycoprotein()->getResidues());
        glycosites_.emplace_back(glycositeResidue, carb, highestResidueNumber);
        for (auto& linkage : carb->GetGlycosidicLinkages())
        {
            cds::determineResiduesForOverlapCheck(linkage); // Now that the protein residue is attached.
        }
        //	    std::cout << "Done with glycan" << std::endl;
        gmml::log(__LINE__, __FILE__, gmml::INF,
                  "Completed creating glycosite on residue " + glycositeInput.proteinResidueId_ + " with glycan " +
                      glycositeInput.glycanInput_);
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
        auto& glycosite     = glycosites_[n];
        auto glycanResidues = glycosite.GetGlycan()->getResidues();
        cds::ResiduesWithOverlapWeight glycanResiduesWithWeight {glycanResidues,
                                                                 std::vector<double>(glycanResidues.size(), 1.0)};
        logss << "Residue ID: " << glycosite.GetResidue()->getStringId()
              << ", overlap: " << cds::CountOverlappingAtoms(overlapResidues_[n], glycanResiduesWithWeight).count;
        gmml::log(__LINE__, __FILE__, gmml::INF, logss.str());
    }
    return;
}
