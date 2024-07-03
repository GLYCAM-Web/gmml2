#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycoproteinBuilder.hpp"
#include "includes/CentralDataStructure/cdsFunctions/bondByDistance.hpp"
#include "includes/CentralDataStructure/cdsFunctions/atomicConnectivity.hpp"
#include "includes/CentralDataStructure/Writers/offWriter.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/gpInputStructs.hpp"
#include "includes/CodeUtils/metropolisCriterion.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/files.hpp"
#include "includes/CodeUtils/strings.hpp" // split
#include "includes/CodeUtils/containers.hpp"
#include "includes/CentralDataStructure/Writers/pdbWriter.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbFile.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbResidue.hpp"
#include "includes/CentralDataStructure/Selections/residueSelections.hpp" // selectResiduesByType
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
    this->WritePdbFile("glycoprotein_initial");
    this->ResolveOverlapsWithWiggler();
    this->WritePdbFile("glycoprotein");
    this->WriteOffFile("glycoprotein");
    this->WritePdbFile("glycoprotein_serialized");
    while (this->GetNumberOfOuputtedStructures() < this->GetNumberOfStructuresToOutput())
    {
        for (auto& glycosite : this->GetGlycosites())
        {
            // gmml::log(__LINE__, __FILE__, gmml::INF, "Setting random dihedrals for " + glycosite.GetResidueId());
            glycosite.SetRandomDihedralAnglesUsingMetadata();
        }
        this->ResolveOverlapsWithWiggler();
        std::stringstream prefix;
        prefix << this->GetNumberOfOuputtedStructures() << "_glycoprotein";
        this->WritePdbFile(prefix.str(), true);
        this->incrementOutputtedStructures();
    }
}

void GlycoproteinBuilder::ResolveOverlapsWithWiggler()
{
    bool randomize = !this->GetIsDeterministic();
    if (randomize)
    {                               // First try a very fast/cheap approach
        if (this->DumbRandomWalk()) // returns true if it fully resolves overlaps.
        {
            return;
        }
    }
    bool useMonteCarlo            = true;
    bool wiggleFirstLinkageOnly   = true; // Only happens when passed to Wiggle function.
    bool useAllResiduesForOverlap = true; // Only happens when passed to Wiggle function.
    unsigned int currentOverlap   = this->Wiggle(this->GetPersistCycles(), wiggleFirstLinkageOnly);
    gmml::log(__LINE__, __FILE__, gmml::INF, "1. Overlap: " + std::to_string(currentOverlap));
    if (randomize)
    {
        currentOverlap = this->RandomDescent(this->GetPersistCycles(), useMonteCarlo);
        gmml::log(__LINE__, __FILE__, gmml::INF, "1r. Overlap: " + std::to_string(currentOverlap));
    }
    currentOverlap = this->Wiggle(this->GetPersistCycles());
    gmml::log(__LINE__, __FILE__, gmml::INF, "2. Overlap: " + std::to_string(currentOverlap));
    currentOverlap = this->Wiggle(this->GetPersistCycles(), wiggleFirstLinkageOnly, 5, useAllResiduesForOverlap);
    gmml::log(__LINE__, __FILE__, gmml::INF, "3. Overlap: " + std::to_string(currentOverlap));
    this->PrintDihedralAnglesAndOverlapOfGlycosites();
    return;
}

unsigned int GlycoproteinBuilder::RandomDescent(int persistCycles, bool monte_carlo)
{
    std::stringstream logss;
    logss << "Random Decent, persisting for " << persistCycles << " cycles and monte carlo is set as " << std::boolalpha
          << monte_carlo << ".\n";
    gmml::log(__LINE__, __FILE__, gmml::INF, logss.str());
    int cycle = 1;
    unsigned int previous_glycan_overlap, new_glycan_overlap, previous_protein_overlap, new_protein_overlap;
    unsigned int lowest_global_overlap = this->CountOverlaps(ALL);
    unsigned int new_global_overlap;
    int overlap_difference                              = 0;
    std::vector<GlycosylationSite*> sites_with_overlaps = this->DetermineSitesWithOverlap(ALL);
    if (sites_with_overlaps.size() == 0)
    {
        gmml::log(__LINE__, __FILE__, gmml::INF, "Stopping RandomDesent with all overlaps resolved.");
        return lowest_global_overlap;
    }
    // Make a random number engine for std::shuffle;
    pcg_extras::seed_seq_from<std::random_device> metropolis_seed_source;
    pcg32 rng_engine(metropolis_seed_source);
    while (cycle < persistCycles)
    {
        gmml::log(__LINE__, __FILE__, gmml::INF,
                  "Cycle " + std::to_string(cycle) + "/" + std::to_string(persistCycles));
        ++cycle;
        std::shuffle(sites_with_overlaps.begin(), sites_with_overlaps.end(), rng_engine);
        for (auto& current_glycosite : sites_with_overlaps)
        {
            previous_glycan_overlap  = current_glycosite->CountOverlaps(GLYCAN);
            previous_protein_overlap = current_glycosite->CountOverlaps(PROTEIN);
            current_glycosite->SetRandomDihedralAnglesUsingMetadata();
            // logss << "Site: " << current_glycosite->GetResidueNumber() << "\n";
            new_glycan_overlap  = current_glycosite->CountOverlaps(GLYCAN);
            new_protein_overlap = current_glycosite->CountOverlaps(PROTEIN);
            overlap_difference  = (new_glycan_overlap + (new_protein_overlap * 5)) -
                                 (previous_glycan_overlap + (previous_protein_overlap * 5));
            if (overlap_difference >= 0) // if the change made it worse
            {
                current_glycosite->ResetDihedralAngles();
            }
            else if ((monte_carlo) && (!monte_carlo::accept_via_metropolis_criterion(overlap_difference)))
            {
                current_glycosite->ResetDihedralAngles();
            }
            else
            {
                gmml::log(__LINE__, __FILE__, gmml::INF,
                          "RandomDescent accepted a change of " + std::to_string(overlap_difference));
            }
        }
        new_global_overlap = this->CountOverlaps(ALL);
        sites_with_overlaps =
            this->DetermineSitesWithOverlap(ALL); // Moved glycans may clash with other glycans. Need to check.
        if (sites_with_overlaps.size() == 0)
        {
            gmml::log(__LINE__, __FILE__, gmml::INF, "Stopping RandomDesent with all overlaps resolved.");
            return new_global_overlap;
        }
        if (lowest_global_overlap > new_global_overlap + 1)
        {
            lowest_global_overlap = new_global_overlap;
            cycle                 = 1;
        }
    }
    return lowest_global_overlap;
}

unsigned int GlycoproteinBuilder::Wiggle(int persistCycles, bool firstLinkageOnly, int interval,
                                         bool useAllResiduesForOverlap)
{
    std::vector<GlycosylationSite*> sites_with_overlaps;
    int cycle            = 0;
    unsigned int overlap = this->CountOverlaps(ALL);
    while (cycle < persistCycles)
    {
        ++cycle;
        gmml::log(__LINE__, __FILE__, gmml::INF,
                  "Cycle " + std::to_string(cycle) + "/" + std::to_string(persistCycles) +
                      ". Overlap: " + std::to_string(overlap));
        sites_with_overlaps = this->DetermineSitesWithOverlap(ALL);
        if (sites_with_overlaps.size() == 0)
        {
            return overlap;
        }
        std::random_shuffle(sites_with_overlaps.begin(), sites_with_overlaps.end());
        for (auto& glycosite : sites_with_overlaps)
        {
            glycosite->Wiggle(firstLinkageOnly, interval, useAllResiduesForOverlap);
        }
        unsigned int newOverlap = this->CountOverlaps(ALL);
        if (overlap > newOverlap)
        {
            overlap = newOverlap;
            cycle   = 1;
        }
    }
    return overlap;
}

bool GlycoproteinBuilder::DumbRandomWalk(int maxCycles)
{
    gmml::log(__LINE__, __FILE__, gmml::INF, "Starting DumbRandomWalk.");
    int cycle                                           = 1;
    std::vector<GlycosylationSite*> sites_with_overlaps = DetermineSitesWithOverlap();
    while (cycle < maxCycles)
    {
        ++cycle;
        for (auto& currentGlycosite : sites_with_overlaps)
        {
            currentGlycosite->SetRandomDihedralAnglesUsingMetadata();
        }
        gmml::log(__LINE__, __FILE__, gmml::INF, "Updating list of sites with overlaps.");
        sites_with_overlaps = this->DetermineSitesWithOverlap();
        if (sites_with_overlaps.empty())
        {
            gmml::log(__LINE__, __FILE__, gmml::INF, "DumbRandomWalk resolved the overlaps. Stopping.");
            return true;
        }
    }
    gmml::log(__LINE__, __FILE__, gmml::INF, "DumbRandomWalk did not resolve the overlaps.");
    return false;
}

void GlycoproteinBuilder::CreateGlycosites(std::vector<glycoprotein::GlycositeInput> glycositesInputVector)
{
    std::vector<Residue*> proteinResidues = this->getGlycoprotein()->getResidues(); // Before any glycans are added.
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
        std::vector<Residue*> otherResidues = proteinResidues;
        otherResidues.erase(std::remove(otherResidues.begin(), otherResidues.end(), glycositeResidue),
                            otherResidues.end());
        gmml::log(__LINE__, __FILE__, gmml::INF,
                  "About to emplace_back to glycosites with: " + glycositeInput.proteinResidueId_ + " and glycan " +
                      glycositeInput.glycanInput_);
        unsigned int highestResidueNumber =
            cdsSelections::findHighestResidueNumber(this->getGlycoprotein()->getResidues());
        glycosites_.emplace_back(glycositeResidue, carb, otherResidues, highestResidueNumber);
        //	    std::cout << "Done with glycan" << std::endl;
        gmml::log(__LINE__, __FILE__, gmml::INF,
                  "Completed creating glycosite on residue " + glycositeInput.proteinResidueId_ + " with glycan " +
                      glycositeInput.glycanInput_);
    }
    //    std::cout << "Done attaching all glycans" << std::endl;
    this->SetOtherGlycosites();
    this->AddOtherGlycositesToLinkageOverlapAtoms();
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

void GlycoproteinBuilder::SetOtherGlycosites()
{
    std::vector<GlycosylationSite>& glycosites    = this->GetGlycosites();
    std::vector<GlycosylationSite*> glycositePtrs = codeUtils::pointerVector(glycosites);
    for (size_t n = 0; n < glycosites.size(); n++)
    {
        glycosites[n].SetOtherGlycosites(codeUtils::withoutNth(n, glycositePtrs));
    }
}

void GlycoproteinBuilder::AddOtherGlycositesToLinkageOverlapAtoms()
{
    for (auto& glycosite : this->GetGlycosites())
    {
        glycosite.AddOtherGlycositesToLinkageOverlapAtoms();
    }
    return;
}

void GlycoproteinBuilder::PrintDihedralAnglesAndOverlapOfGlycosites()
{
    for (auto& glycosite : this->GetGlycosites())
    {
        glycosite.Print("All");
    }
    return;
}

void GlycoproteinBuilder::SetRandomDihedralAnglesUsingMetadata()
{
    for (auto& glycosite : this->GetGlycosites())
    {
        glycosite.SetRandomDihedralAnglesUsingMetadata();
    }
    return;
}

unsigned int GlycoproteinBuilder::CountOverlaps(MoleculeType)
{
    unsigned int overlap = 0;
    for (auto& glycosite : this->GetGlycosites())
    {
        overlap += glycosite.CountOverlapsFast();
    }
    return overlap;
}

std::vector<GlycosylationSite*> GlycoproteinBuilder::DetermineSitesWithOverlap(MoleculeType)
{
    std::vector<GlycosylationSite*> sites_to_return;
    unsigned int overlap = 0;
    for (auto& glycosite : this->GetGlycosites())
    {
        overlap = glycosite.CountOverlapsFast();
        if (overlap > this->GetOverlapTolerance())
        {
            sites_to_return.push_back(&glycosite);
        }
    }
    return sites_to_return;
}
