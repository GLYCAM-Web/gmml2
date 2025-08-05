#include "include/CentralDataStructure/cdsFunctions.hpp"
#include "include/CentralDataStructure/graphInterface.hpp"
#include "include/CentralDataStructure/molecule.hpp"
#include "include/CentralDataStructure/offWriter.hpp"
#include "include/CentralDataStructure/pdbWriter.hpp"
#include "include/assembly/assemblyGraph.hpp"
#include "include/fileType/lib/libraryFile.hpp"
#include "include/fileType/off/offFileWriter.hpp"
#include "include/readers/parameterManager.hpp"
#include "include/util/files.hpp"
#include "include/util/logging.hpp"

#include <iostream>
#include <ostream>
#include <stdexcept>

int main()
{
    using namespace gmml;
    std::string libFilePath = "../dat/CurrentParams/leaprc.ff12SB_2014-04-24/amino12.lib";
    lib::LibraryData data = lib::loadLibraryData(libFilePath);
    ParameterManager params {data};
    std::cout << "Finished loading libfile" << std::endl;
    Molecule molecule = Molecule();
    for (size_t n = 0; n < data.residues.size(); n++)
    {
        const std::string& name = data.residueNames[n];
        Residue* residue = molecule.addResidue(std::make_unique<Residue>());
        residue->setName(name);
        residue->determineType(name);
        createAtomsForResidue(params, residue, name);
    }
    GraphIndexData moleculeIndices = toIndexData({&molecule});
    // Need a central place for this:
    std::string fileName = "./libAsPdbFile.pdb";
    try
    {
        util::writeToFile(fileName, [&](std::ostream& stream) { WritePdb(stream, moleculeIndices, {}); });
    }
    catch (...)
    {
        util::log(__LINE__, __FILE__, util::ERR, "Error when writing pdbFile class to file:\n" + fileName);
        throw std::runtime_error("Error when writing pdbFile class to file:\n" + fileName);
    }
    // OFF molecule
    fileName = "./libAsOffFile.off";
    try
    {
        util::writeToFile(fileName, [&](std::ostream& stream) { WriteOff(stream, "MOLECULE", moleculeIndices); });
    }
    catch (...)
    {
        util::log(__LINE__, __FILE__, util::ERR, "Error when writing to file:\n" + fileName);
        throw std::runtime_error("Error when writing to file:\n" + fileName);
    }
    // OFF separate residues
    fileName = "./libAsLibFile.lib";
    try
    {
        GraphIndexData graphData = toIndexData(molecule.getResidues());
        assembly::Graph graph = createVisibleAssemblyGraph(graphData);
        serializeResiduesIndividually(graphData.objects.residues);
        util::writeToFile(
            fileName,
            [&](std::ostream& stream)
            { off::writeResiduesIndividually(stream, graph, toOffFileData(graphData.objects.residues)); });
    }
    catch (...)
    {
        util::log(__LINE__, __FILE__, util::ERR, "Error when writing to file:\n" + fileName);
        throw std::runtime_error("Error when writing to file:\n" + fileName);
    }
}
