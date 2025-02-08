#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycoproteinCreation.hpp"
#include "includes/CentralDataStructure/Shapers/residueLinkageCreation.hpp"
#include "includes/CentralDataStructure/Selections/residueSelections.hpp"
#include "includes/CentralDataStructure/assembly.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbResidue.hpp"
#include "includes/CodeUtils/casting.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/strings.hpp"

#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <stdexcept>

namespace glycoproteinBuilder
{
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
            if (residue->GetType() != cds::ResidueType::Protein)
            {
                gmml::log(__LINE__, __FILE__, gmml::INF,
                          "Glycosite selector skipping non-protein residue: " + residue->getStringId());
                continue;
            }
            pdb::PdbResidue* pdbResidue = codeUtils::erratic_cast<pdb::PdbResidue*>(residue);
            //            std::cout << pdbResidue->getChainId() << "_";
            //            std::cout << pdbResidue->getNumberAndInsertionCode() << "\n";
            if ((pdbResidue->getChainId() == userSelectedChain) &&
                (pdbResidue->getNumberAndInsertionCode() == userSelectedResidue))
            {
                gmml::log(__LINE__, __FILE__, gmml::INF, "Id of selected glycosite: " + pdbResidue->printId());
                return residue;
            }
        }
        return nullptr;
    }

    std::vector<GlycosylationSite>
    createGlycosites(cds::Assembly* glycoprotein,
                     std::vector<glycoproteinBuilder::GlycositeInput> glycositesInputVector)
    {
        std::vector<GlycosylationSite> glycosites;
        for (auto& glycositeInput : glycositesInputVector)
        {
            gmml::log(__LINE__, __FILE__, gmml::INF,
                      "Creating glycosite on residue " + glycositeInput.proteinResidueId + " with glycan " +
                          glycositeInput.glycanInput);
            Residue* glycositeResidue = selectResidueFromInput(glycoprotein, glycositeInput.proteinResidueId);
            if (glycositeResidue == nullptr)
            {
                throw std::runtime_error("Error: Did not find a residue with id matching " +
                                         glycositeInput.proteinResidueId + "\n");
            }
            Carbohydrate* carb = codeUtils::erratic_cast<Carbohydrate*>(
                glycoprotein->addMolecule(std::make_unique<Carbohydrate>(glycositeInput.glycanInput)));
            gmml::log(__LINE__, __FILE__, gmml::INF,
                      "About to emplace_back to glycosites with: " + glycositeInput.proteinResidueId + " and glycan " +
                          glycositeInput.glycanInput);
            unsigned int highestResidueNumber = cdsSelections::findHighestResidueNumber(glycoprotein->getResidues());
            glycosites.emplace_back(glycositeResidue, carb, glycositeInput, highestResidueNumber);
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
} // namespace glycoproteinBuilder