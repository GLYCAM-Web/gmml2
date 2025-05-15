#include "includes/CentralDataStructure/Readers/Pdb/pdbFile.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbResidue.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbSelections.hpp"
#include "includes/CentralDataStructure/Selections/residueSelections.hpp"
#include "includes/CentralDataStructure/Selections/templatedSelections.hpp"
#include "includes/CentralDataStructure/Writers/pdbWriter.hpp"
#include "includes/CodeUtils/casting.hpp"
#include "includes/CodeUtils/files.hpp"

#include <vector>
#include <string>
#include <iostream>
#include <ostream>

int main(int argc, char* argv[])
{
    if (argc != 2)
    {
        std::cout << "Usage: " << argv[0] << " inputFile.pdb\n";
        std::cout << "Example: " << argv[0] << " tests/inputs/4mbz.pdb\n";
        std::exit(EXIT_FAILURE);
    }
    // requirement: a chain ID for every single ATOM entry, and all ligand atoms should be put in a single residue.
    pdb::PdbFile pdbFile(argv[1]); // PdbFile is an "Ensemble" (made up of "Assemblies"), but if you want to just set
                                   // every molecule to have any chain ID you can do:
    pdbFile.data.residues.chainIds = std::vector<std::string>(pdbFile.data.indices.residueCount, "Y");
    // ResidueTypes are guessed upon input. Using that guess to find the ligand, can improve this if you need:
    std::vector<size_t> ligandResidues =
        codeUtils::indicesOfElement(pdbFile.data.residues.types, cds::ResidueType::Undefined);
    if (ligandResidues.empty())
    {
        std::cout << "No ligand residues found in input file\n";
        return 0;
    }
    size_t firstLigandResidue    = ligandResidues.front();
    const std::string& firstName = pdbFile.data.residues.names[firstLigandResidue];
    uint firstNumber             = pdbFile.data.residues.numbers[firstLigandResidue];
    for (size_t ligandResidue : ligandResidues) // Each MODEL in PdbFile is converted into an "Assembly"
    {                                           // Every ligand residue gets the same residue number as the first one.
        // std::cout << "Renumbering and renaming " << ligandResidue->getStringId() << "\n";
        cds::Residue* residue = pdbFile.data.objects.residues[ligandResidue];
        residue->setName(firstName);
        residue->setNumber(firstNumber);
    }
    pdbFile.Write("./026.outputPdbFile.pdb");
    // ************************************************************************ //
    // Separate thing showing how to read/write PDB files as "trajectories/frames"
    pdb::PdbFile pdbFileTraj(argv[1], pdb::InputType::modelsAsCoordinates);
    std::vector<cds::Residue*> myResidues = pdb::getResidues(pdbFileTraj.getAssemblies());
    // somehow you specify number in inputs. e.g. A_405 chain A, residue 405.
    cds::Residue* queryResidue            = codeUtils::findElementWithNumber(myResidues, 5);
    double distance                       = 12.345; // inputs
    std::vector<cds::Residue*> selectedResidues =
        cdsSelections::selectResiduesWithinDistanceN(myResidues, queryResidue, distance);
    std::cout << "Found " << selectedResidues.size() << " residues\n";
    cds::Assembly newAssembly(selectedResidues);
    const std::string outName = "026.outputSelection.pdb";
    codeUtils::writeToFile(outName,
                           [&](std::ostream& stream)
                           {
                               cds::writeTrajectoryToPdb(stream, newAssembly.getMolecules());
                           });
    return 0;
}
