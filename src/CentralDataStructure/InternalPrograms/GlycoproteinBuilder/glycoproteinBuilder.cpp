#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycoproteinBuilder.hpp"
#include "includes/CentralDataStructure/cdsFunctions/bondByDistance.hpp"
#include "includes/CentralDataStructure/cdsFunctions/atomicConnectivity.hpp"
#include "includes/CentralDataStructure/Writers/offWriter.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/gpInputStructs.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/files.hpp"
#include "includes/CodeUtils/strings.hpp" // split
#include "includes/CodeUtils/containers.hpp"
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

    ResiduesByType overlapResidues;
    for (size_t n = 0; n < glycosites_.size(); n++)
    {
        std::vector<Residue*> protein = proteinResidues_;
        protein.erase(std::remove(protein.begin(), protein.end(), glycosites_[n].GetResidue()), protein.end());

        std::vector<Residue*> glycan;
        for (auto& otherSite : codeUtils::withoutNth(n, glycositeResidues))
        {
            codeUtils::insertInto(glycan, otherSite);
        }

        overlapResidues.protein.push_back(protein);
        overlapResidues.glycan.push_back(glycan);
        overlapResidues.all.push_back(codeUtils::vectorAppend(protein, glycan));
    }
    overlapResidues_ = overlapResidues;

    bool random          = !GetIsDeterministic();
    int persistCycles    = GetPersistCycles();
    int overlapTolerance = GetOverlapTolerance();

    this->WritePdbFile("glycoprotein_initial");
    resolveOverlapsWithWiggler(random, persistCycles, overlapTolerance, glycosidicLinkages, overlapResidues,
                               glycositeResidues);
    this->PrintDihedralAnglesAndOverlapOfGlycosites();
    this->WritePdbFile("glycoprotein");
    this->WriteOffFile("glycoprotein");
    this->WritePdbFile("glycoprotein_serialized");
    while (this->GetNumberOfOuputtedStructures() < this->GetNumberOfStructuresToOutput())
    {
        for (auto& linkages : glycosidicLinkages)
        {
            cds::setRandomShapeUsingMetadata(linkages);
        }
        resolveOverlapsWithWiggler(random, persistCycles, overlapTolerance, glycosidicLinkages, overlapResidues,
                                   glycositeResidues);
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
        auto& glycosite = glycosites_[n];
        logss << "Residue ID: " << glycosite.GetResidue()->getStringId() << ", overlap: "
              << countGlycositeOverlaps(overlapResidues_.all[n], glycosite.GetGlycan()->getResidues()).count;
        gmml::log(__LINE__, __FILE__, gmml::INF, logss.str());
    }
    return;
}
