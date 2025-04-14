#include "includes/CentralDataStructure/Readers/Prep/prepFile.hpp"
#include "includes/CodeUtils/files.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CentralDataStructure/FileFormats/offFileWriter.hpp"
#include "includes/CentralDataStructure/Writers/pdbWriter.hpp"
#include "includes/CentralDataStructure/Writers/offWriter.hpp"
#include "includes/CentralDataStructure/Editors/glycamResidueCombinator.hpp"
#include "includes/CentralDataStructure/cdsFunctions/cdsFunctions.hpp"
#include "includes/CentralDataStructure/cdsFunctions/graphInterface.hpp"
#include "includes/Assembly/assemblyGraph.hpp"

#include <ostream>
#include <stdexcept>

int main()
{
    std::string prepFilePath                        = "../dat/prep/GLYCAM_06j-1_GAGS.prep";
    std::vector<std::string> residuesToLoadFromPrep = {"0GA", "4YB", "4uA", "Cake", "4YA"};
    cds::Molecule molecule                          = cds::Molecule();
    prep::readPrepFile(&molecule, prepFilePath, residuesToLoadFromPrep);
    // PREP residues
    prep::Write(&molecule, "./prepAsPrepFile.prep");
    // PDB
    std::string fileName = "./prepAsPdbFile.pdb";
    try
    {
        cds::GraphIndexData indices = cds::toIndexData({&molecule});
        codeUtils::writeToFile(fileName,
                               [&](std::ostream& stream)
                               {
                                   cds::WritePdb(stream, indices, {});
                               });
    }
    catch (...)
    {
        gmml::log(__LINE__, __FILE__, gmml::ERR, "Error when writing pdbFile class to file:\n" + fileName);
        throw std::runtime_error("Error when writing pdbFile class to file:\n" + fileName);
    }
    // OFF molecule
    fileName = "./prepAsOffFile.off";
    try
    {
        cds::GraphIndexData indices = cds::toIndexData(molecule.getResidues());
        assembly::Graph graph       = cds::createVisibleAssemblyGraph(indices);
        cds::serializeResiduesIndividually(indices.residues);
        codeUtils::writeToFile(fileName,
                               [&](std::ostream& stream)
                               {
                                   cds::WriteResiduesIndividuallyToOffFile(stream, graph,
                                                                           cds::toOffFileData(indices.residues));
                               });
    }
    catch (...)
    {
        gmml::log(__LINE__, __FILE__, gmml::ERR, "Error when writing to file:\n" + fileName);
        throw std::runtime_error("Error when writing to file:\n" + fileName);
    }
    // OFF separate residues
    molecule.setName("LIBRARY");
    fileName = "./prepAsLibFile.lib";
    try
    {
        cds::GraphIndexData indices = cds::toIndexData(molecule.getResidues());
        assembly::Graph graph       = cds::createVisibleAssemblyGraph(indices);
        cds::serializeResiduesIndividually(indices.residues);
        codeUtils::writeToFile(fileName,
                               [&](std::ostream& stream)
                               {
                                   cds::WriteResiduesIndividuallyToOffFile(stream, graph,
                                                                           cds::toOffFileData(indices.residues));
                               });
    }
    catch (...)
    {
        gmml::log(__LINE__, __FILE__, gmml::ERR, "Error when writing to file:\n" + fileName);
        throw std::runtime_error("Error when writing to file:\n" + fileName);
    }
}
