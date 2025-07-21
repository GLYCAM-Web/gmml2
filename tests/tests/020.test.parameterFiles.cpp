#include "include/CentralDataStructure/Editors/glycamResidueCombinator.hpp"
#include "include/CentralDataStructure/FileFormats/offFileWriter.hpp"
#include "include/CentralDataStructure/cdsFunctions.hpp"
#include "include/CentralDataStructure/graphInterface.hpp"
#include "include/CentralDataStructure/molecule.hpp"
#include "include/assembly/assemblyGraph.hpp"
#include "include/readers/Prep/prepFile.hpp"
#include "include/util/files.hpp"
#include "include/util/logging.hpp"
#include "include/writers/offWriter.hpp"
#include "include/writers/pdbWriter.hpp"

#include <ostream>
#include <stdexcept>

int main()
{
    using namespace gmml;
    std::string prepFilePath = "../dat/prep/GLYCAM_06j-1_GAGS.prep";
    std::vector<std::string> residuesToLoadFromPrep = {"0GA", "4YB", "4uA", "Cake", "4YA"};
    Molecule molecule = Molecule();
    prep::readPrepFile(&molecule, prepFilePath, residuesToLoadFromPrep);
    // PREP residues
    prep::Write(&molecule, "./prepAsPrepFile.prep");
    // PDB
    std::string fileName = "./prepAsPdbFile.pdb";
    try
    {
        GraphIndexData indices = toIndexData({&molecule});
        util::writeToFile(fileName, [&](std::ostream& stream) { WritePdb(stream, indices, {}); });
    }
    catch (...)
    {
        util::log(__LINE__, __FILE__, util::ERR, "Error when writing pdbFile class to file:\n" + fileName);
        throw std::runtime_error("Error when writing pdbFile class to file:\n" + fileName);
    }
    // OFF molecule
    fileName = "./prepAsOffFile.off";
    try
    {
        GraphIndexData graphData = toIndexData(molecule.getResidues());
        assembly::Graph graph = createVisibleAssemblyGraph(graphData);
        serializeResiduesIndividually(graphData.objects.residues);
        util::writeToFile(
            fileName,
            [&](std::ostream& stream)
            { WriteResiduesIndividuallyToOffFile(stream, graph, toOffFileData(graphData.objects.residues)); });
    }
    catch (...)
    {
        util::log(__LINE__, __FILE__, util::ERR, "Error when writing to file:\n" + fileName);
        throw std::runtime_error("Error when writing to file:\n" + fileName);
    }
    // OFF separate residues
    molecule.setName("LIBRARY");
    fileName = "./prepAsLibFile.lib";
    try
    {
        GraphIndexData graphData = toIndexData(molecule.getResidues());
        assembly::Graph graph = createVisibleAssemblyGraph(graphData);
        serializeResiduesIndividually(graphData.objects.residues);
        util::writeToFile(
            fileName,
            [&](std::ostream& stream)
            { WriteResiduesIndividuallyToOffFile(stream, graph, toOffFileData(graphData.objects.residues)); });
    }
    catch (...)
    {
        util::log(__LINE__, __FILE__, util::ERR, "Error when writing to file:\n" + fileName);
        throw std::runtime_error("Error when writing to file:\n" + fileName);
    }
}
