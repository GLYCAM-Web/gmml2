#include "include/CentralDataStructure/cdsFunctions.hpp"
#include "include/CentralDataStructure/graphInterface.hpp"
#include "include/CentralDataStructure/offWriter.hpp"
#include "include/assembly/assemblyGraph.hpp"
#include "include/assembly/assemblyIndices.hpp"
#include "include/assembly/assemblyTypes.hpp"
#include "include/fileType/off/offFileWriter.hpp"
#include "include/fileType/pdb/bondByDistance.hpp"
#include "include/fileType/pdb/pdbFile.hpp"
#include "include/fileType/pdb/pdbFunctions.hpp"
#include "include/preprocess/parameterManager.hpp"
#include "include/preprocess/pdbPreProcess.hpp"
#include "include/preprocess/pdbPreprocessorInputs.hpp"
#include "include/util/containers.hpp"
#include "include/util/files.hpp"
#include "include/util/filesystem.hpp"
#include "include/util/logging.hpp"

#include <fstream>
#include <iostream>
#include <string>

int main(int argc, char* argv[])
{
    using namespace gmml;

    if (argc != 3)
    {

        std::cout << "Usage: " << argv[0] << " inputFile.pdb outputFile.pdb\n";
        std::cout << "Example: " << argv[0] << " inputs/4mbz.pdb output/4mbz.pdb\n";
        std::exit(EXIT_FAILURE);
    }
    std::string inputFile = argv[1];
    std::string outputFile = argv[2];
    pdb::PdbFile pdbFile = pdb::toPdbFile(argv[1], pdb::modelsAsMolecules);
    std::string baseDir = util::toString(util::pathAboveCurrentExecutableDir());
    preprocess::PreprocessorOptions options = preprocess::defaultPreprocessorOptions; // Default values are good.
    std::cout << "Preprocessing\n";
    const preprocess::ParameterManager parameterManager = preprocess::loadParameters(baseDir);
    preprocess::PreprocessorInformation ppInfo = preProcess(pdbFile, parameterManager, options);
    std::vector<Assembly*> assemblies = getAssemblies(pdbFile);
    for (size_t assemblyId = 0; assemblyId < pdbFile.data.indices.assemblyCount; assemblyId++)
    {
        std::vector<size_t> atomIds = assemblyAtoms(pdbFile.data.indices, assemblyId);
        std::cout << "Bonding atoms by distance for assembly" << std::endl;
        pdb::bondAtomsByDistance(pdbFile.data, atomIds);
        // OFF molecule
        try
        {
            GraphIndexData graphData = toIndexData(assemblies[assemblyId]->getMolecules());
            assembly::Graph graph = createVisibleAssemblyGraph(graphData);
            off::OffFileData data = toOffFileData(graphData.objects.residues);
            serializeNumbers(graphData.objects.atoms);
            serializeNumbers(graphData.objects.residues);
            util::writeToFile(
                baseDir + "/outputOffFile.off",
                [&](std::ostream& stream) {
                    off::writeResiduesTogether(
                        stream, graph, data, util::indexVector(graphData.objects.residues), "Assembly");
                });
        }
        catch (std::runtime_error& error)
        {
            std::stringstream ss;
            ss << "Runtime error thrown when writing to off file:\n" << error.what() << "\n";
            std::cout << ss.str();
            util::log(__LINE__, __FILE__, util::ERR, ss.str());
        }
        catch (...)
        {
            std::cout << "Unknown error when writing to off file.\n";
            util::log(__LINE__, __FILE__, util::ERR, "Unknown error when writing to off file.\n");
        }
    }

    std::cout << "Finished bonding atoms by distance" << std::endl;
    pdbFile.data.atoms.numbers = atomNumbers(pdbFile.data.objects.atoms);
    pdbFile.data.residues.numbers = residueNumbers(pdbFile.data.objects.residues);
    write(pdbFile, outputFile);

    auto residueIdString = [](const pdb::ResidueId& id)
    { return id.residueName + " | " + id.chainId + " | " + id.sequenceNumber + id.insertionCode; };

    auto residueIdStringChainFirst = [](const pdb::ResidueId& id)
    { return id.chainId + " | " + id.residueName + " | " + id.sequenceNumber + id.insertionCode; };

    auto atomInfoString = [&](const preprocess::AtomInfo& atom)
    { return atom.name + " | " + residueIdString(atom.residue); };

    // Just showing what's in the ppInfo and how to access it
    std::cout << "Unrecognized atoms:\n";
    for (auto& unrecognized : ppInfo.unrecognizedAtoms)
    {
        std::cout << atomInfoString(unrecognized) << "\n";
    }
    std::cout << "Missing heavy atoms:\n";
    for (auto& missing : ppInfo.missingHeavyAtoms)
    {
        std::cout << atomInfoString(missing) << "\n";
    }
    std::cout << "Unrecognized residues:\n";
    for (auto& unrecognized : ppInfo.unrecognizedResidues)
    {
        std::cout << residueIdString(unrecognized) << "\n";
    }
    std::cout << "Gaps in amino acid chain:\n";
    for (auto& gap : ppInfo.missingResidues)
    {
        std::cout << gap.chainId << " | " << gap.residueBeforeGap << " | " << gap.residueAfterGap << " | "
                  << gap.terminationBeforeGap << " | " << gap.terminationAfterGap << "\n";
    }
    std::cout << "Histidine Protonation:\n";
    for (auto& his : ppInfo.hisResidues)
    {
        std::cout << residueIdString(his) << "\n";
    }
    std::cout << "Disulphide bonds:\n";
    for (auto& cysBond : ppInfo.cysBondResidues)
    {
        std::cout << residueIdStringChainFirst(cysBond.residue1) << " | " << cysBond.distance << " | "
                  << residueIdStringChainFirst(cysBond.residue2) << "\n";
    }
    std::cout << "Chain terminations:\n";
    for (auto& chainT : ppInfo.chainTerminals)
    {
        std::cout << chainT.chainId << " | " << chainT.startIndex << " | " << chainT.nTermination << " | "
                  << chainT.endIndex << " | " << chainT.cTermination << "\n";
    }
    std::cout << "NonNatural Protein Residues:\n";
    for (auto& nonNaturalResidue : ppInfo.nonNaturalProteinResidues)
    {
        std::cout << residueIdStringChainFirst(nonNaturalResidue) << "\n";
    }
    return 0;
}
