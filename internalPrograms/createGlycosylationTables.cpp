#include "includes/CentralDataStructure/InternalPrograms/glycosylationSiteFinder.hpp"
#include "includes/CentralDataStructure/cdsFunctions/bondByDistance.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbFile.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbSelections.hpp"

#include <string>
#include <vector>
#include <iostream>

int main(int argc, char* argv[])
{
    if (argc != 2)
    {
        std::cout << "Usage: createGlycosylationTables.exe inputFile.pdb\n";
        std::cout << "Example: pdb2glycam 1RVX.pdb\n";
        std::exit(EXIT_FAILURE);
    }
    std::string inputFileName = argv[1];
    pdb::PdbFile inputFile(inputFileName);
    auto residues = pdb::getResidues(inputFile.getAssemblies());
    cds::bondAtomsAndResiduesByDistance(residues);
    glycoproteinBuilder::GlycosylationSiteFinder siteFinder(residues);
    std::cout << siteFinder.PrintTable();
    // std::vector<GlycosylationSiteInfo> tableInfo = siteFinder.GetTable();
    //    for (auto &tableElement : tableInfo)
    //    {
    //        std::cout << tableElement.Print() << "\n";
    //    }
    return 0;
}
