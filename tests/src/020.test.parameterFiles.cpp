#include "include/CentralDataStructure/cdsFunctions.hpp"
#include "include/CentralDataStructure/graphInterface.hpp"
#include "include/CentralDataStructure/molecule.hpp"
#include "include/CentralDataStructure/offWriter.hpp"
#include "include/CentralDataStructure/pdbWriter.hpp"
#include "include/assembly/assemblyGraph.hpp"
#include "include/fileType/off/offFileWriter.hpp"
#include "include/fileType/pdb/pdbFileWriter.hpp"
#include "include/fileType/prep/prepFile.hpp"
#include "include/fileType/prep/prepFunctions.hpp"
#include "include/util/containers.hpp"
#include "include/util/files.hpp"
#include "include/util/logging.hpp"

#include <ostream>
#include <stdexcept>

int main()
{
    using namespace gmml;
    std::string prepFilePath = "../dat/prep/GLYCAM_06j-1_GAGS.prep";
    std::vector<std::string> residuesToLoadFromPrep = {"0GA", "4YB", "4uA", "Cake", "4YA"};
    prep::PrepData data = prep::readPrepFile(prepFilePath, residuesToLoadFromPrep);
    // PREP residues
    prep::write(data, "./prepAsPrepFile.prep");
    // PDB
    assembly::Indices indices {
        data.atomCount,
        data.residueCount,
        1,
        1,
        data.atomGraph.nodeAlive,
        data.atomResidue,
        std::vector<size_t>(data.residueCount, 0),
        {0}};
    assembly::Graph graph = assembly::createAssemblyGraph(indices, data.atomGraph);
    std::vector<bool> residueTER(data.residueCount, true);

    std::function<std::string(const size_t&)> toElement = [&](size_t n)
    { return (data.atomGraph.nodeAlive[n] && !data.atoms.name[n].empty()) ? std::string {data.atoms.name[n][0]} : ""; };

    std::vector<size_t> atomIndices = util::indexVector(data.atomCount);
    std::vector<std::string> elementStrings = util::vectorMap(toElement, atomIndices);

    std::function<uint(const size_t&)> atomicNumber = [&](size_t n)
    { return data.atomGraph.nodeAlive[n] ? findElementAtomicNumber(elementStrings[n]) : Element::Unknown; };

    std::vector<uint> atomicNumbers = util::vectorMap(atomicNumber, atomIndices);
    std::vector<uint> residueNumbers(data.residueCount, 1);

    pdb::PdbFileAtomData pdbAtomData {
        data.atoms.coordinate,
        data.atoms.number,
        data.atoms.name,
        elementStrings,
        std::vector<std::string>(data.atomCount, "ATOM"),
        std::vector<double>(data.atomCount, 1.0),
        std::vector<double>(data.atomCount, 0.0)};

    pdb::PdbFileResidueData pdbResidueData {
        residueNumbers,
        data.residues.name,
        std::vector<std::string>(data.residueCount, ""),
        std::vector<std::string>(data.residueCount, "")};

    pdb::PdbFileData pdbData {pdb::PdbFileFormat(), {}, pdbResidueData, pdbAtomData};
    std::string fileName = "./prepAsPdbFile.pdb";
    try
    {
        util::writeToFile(
            fileName,
            [&](std::ostream& stream)
            {
                pdb::writeMoleculeToPdb(stream, graph, util::indexVector(data.residueCount), residueTER, pdbData);
                pdb::theEnd(stream);
            });
    }
    catch (...)
    {
        util::log(__LINE__, __FILE__, util::ERR, "Error when writing pdbFile class to file:\n" + fileName);
        throw std::runtime_error("Error when writing pdbFile class to file:\n" + fileName);
    }
    // OFF molecule

    off::OffFileAtomData offAtomData {
        data.atoms.number, data.atoms.name, data.atoms.type, atomicNumbers, data.atoms.charge, data.atoms.coordinate};

    for (size_t n = 0; n < data.residueCount; n++)
    {
        std::vector<size_t> atoms = prep::residueAtoms(data, n);
        util::setIndicesTo(offAtomData.numbers, atoms, serializedNumberVector(atoms.size()));
    }

    off::OffFileResidueData offResidueData {
        residueNumbers,
        data.residues.name,
        std::vector<ResidueType>(data.residueCount, Undefined),
        std::vector<std::vector<size_t>>(data.residueCount)};

    off::OffFileData offData {off::OffFileFormat(), offResidueData, offAtomData};

    fileName = "./prepAsOffFile.off";
    try
    {
        util::writeToFile(
            fileName, [&](std::ostream& stream) { off::writeResiduesIndividually(stream, graph, offData); });
    }
    catch (...)
    {
        util::log(__LINE__, __FILE__, util::ERR, "Error when writing to file:\n" + fileName);
        throw std::runtime_error("Error when writing to file:\n" + fileName);
    }
    // OFF separate residues
    fileName = "./prepAsLibFile.lib";
    try
    {
        util::writeToFile(
            fileName, [&](std::ostream& stream) { off::writeResiduesIndividually(stream, graph, offData); });
    }
    catch (...)
    {
        util::log(__LINE__, __FILE__, util::ERR, "Error when writing to file:\n" + fileName);
        throw std::runtime_error("Error when writing to file:\n" + fileName);
    }
}
