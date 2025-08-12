#include "include/CentralDataStructure/cdsFunctions.hpp"
#include "include/CentralDataStructure/graphInterface.hpp"
#include "include/CentralDataStructure/molecule.hpp"
#include "include/assembly/assemblyGraph.hpp"
#include "include/carbohydrate/carbohydrate.hpp"
#include "include/carbohydrate/offWriter.hpp"
#include "include/carbohydrate/pdbWriter.hpp"
#include "include/fileType/lib/libraryFile.hpp"
#include "include/fileType/off/offFileWriter.hpp"
#include "include/preprocess/parameterManager.hpp"
#include "include/util/containers.hpp"
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
    preprocess::ParameterManager params {data};
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
    carbohydrate::CarbohydrateData carbData = structured({&molecule});
    assembly::Graph graph = createVisibleAssemblyGraph(carbData);
    // Need a central place for this:
    std::string fileName = "./libAsPdbFile.pdb";
    try
    {
        util::writeToFile(fileName, [&](std::ostream& stream) { writePdb(stream, carbData, graph, {}); });
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
        off::OffFileData offData = toOffFileData(carbData, graph);
        util::writeToFile(
            fileName,
            [&](std::ostream& stream)
            { off::writeResiduesTogether(stream, offData, graph, indices(graph.residues.nodes), "MOLECULE"); });
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
        off::OffFileData offData = toOffFileData(carbData, graph);
        offData.residues.numbers = std::vector<uint>(graph.indices.residueCount, 1);
        for (size_t n : indices(graph.residues.nodes))
        {
            std::vector<size_t> residueAtoms = assembly::residueAtoms(graph, n);
            std::vector<bool> alive = util::indicesToValues(graph.atoms.source.nodes.alive, residueAtoms);
            std::vector<size_t> aliveAtoms = util::boolsToValues(residueAtoms, alive);
            util::setIndicesTo(offData.atoms.numbers, aliveAtoms, util::serializedNumberVector(aliveAtoms.size()));
        }
        util::writeToFile(
            fileName, [&](std::ostream& stream) { off::writeResiduesIndividually(stream, offData, graph); });
    }
    catch (...)
    {
        util::log(__LINE__, __FILE__, util::ERR, "Error when writing to file:\n" + fileName);
        throw std::runtime_error("Error when writing to file:\n" + fileName);
    }
}
