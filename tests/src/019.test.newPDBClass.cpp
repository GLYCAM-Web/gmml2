#include "include/CentralDataStructure/FileFormats/offFileWriter.hpp"
#include "include/CentralDataStructure/cdsFunctions.hpp"
#include "include/CentralDataStructure/graphInterface.hpp"
#include "include/assembly/assemblyGraph.hpp"
#include "include/assembly/assemblyIndices.hpp"
#include "include/assembly/assemblyTypes.hpp"
#include "include/readers/Pdb/bondByDistance.hpp"
#include "include/readers/Pdb/pdbFile.hpp"
#include "include/readers/Pdb/pdbFunctions.hpp"
#include "include/readers/Pdb/pdbPreprocessorInputs.hpp"
#include "include/readers/parameterManager.hpp"
#include "include/util/files.hpp"
#include "include/util/filesystem.hpp"
#include "include/util/logging.hpp"
#include "include/writers/offWriter.hpp"

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
    pdb::PdbFile pdbFile(argv[1], {pdb::InputType::modelsAsMolecules, false});
    std::string baseDir = util::toString(util::pathAboveCurrentExecutableDir());
    pdb::PreprocessorOptions options; // Default values are good.
    std::cout << "Preprocessing\n";
    const ParameterManager parameterManager = loadParameters(baseDir);
    pdb::PreprocessorInformation ppInfo = pdbFile.PreProcess(parameterManager, options);
    std::vector<Assembly*> assemblies = pdbFile.getAssemblies();
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
            OffFileData data = toOffFileData(graphData.objects.residues);
            serializeNumbers(graphData.objects.atoms);
            serializeNumbers(graphData.objects.residues);
            util::writeToFile(
                baseDir + "/outputOffFile.off",
                [&](std::ostream& stream) {
                    WriteResiduesTogetherToOffFile(
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
    pdbFile.Write(outputFile);

    // Just showing what's in the ppInfo and how to access it
    std::cout << "Unrecognized atoms:\n";
    for (auto& unrecognized : ppInfo.unrecognizedAtoms_)
    {
        std::cout << unrecognized.name_ << " | " << unrecognized.residue_.getName() << " | "
                  << unrecognized.residue_.getChainId() << " | " << unrecognized.residue_.getNumberAndInsertionCode()
                  << "\n";
    }
    std::cout << "Missing heavy atoms:\n";
    for (auto& missing : ppInfo.missingHeavyAtoms_)
    {
        std::cout << missing.name_ << " | " << missing.residue_.getName() << " | " << missing.residue_.getChainId()
                  << " | " << missing.residue_.getNumberAndInsertionCode() << "\n";
    }
    std::cout << "Unrecognized residues:\n";
    for (auto& unrecognized : ppInfo.unrecognizedResidues_)
    {
        std::cout << unrecognized.getName() << " | " << unrecognized.getChainId() << " | "
                  << unrecognized.getNumberAndInsertionCode() << "\n";
    }
    std::cout << "Gaps in amino acid chain:\n";
    for (auto& gap : ppInfo.missingResidues_)
    {
        std::cout << gap.chainId_ << " | " << gap.residueBeforeGap_ << " | " << gap.residueAfterGap_ << " | "
                  << gap.terminationBeforeGap_ << " | " << gap.terminationAfterGap_ << "\n";
    }
    std::cout << "Histidine Protonation:\n";
    for (auto& his : ppInfo.hisResidues_)
    {
        std::cout << his.getName() << " | " << his.getChainId() << " | " << his.getNumberAndInsertionCode() << "\n";
    }
    std::cout << "Disulphide bonds:\n";
    for (auto& cysBond : ppInfo.cysBondResidues_)
    {
        std::cout << cysBond.residue1_.getChainId() << " | " << cysBond.residue1_.getName() << " | "
                  << cysBond.residue1_.getNumberAndInsertionCode() << " | " << cysBond.distance_ << " | "
                  << cysBond.residue2_.getChainId() << " | " << cysBond.residue2_.getName() << " | "
                  << cysBond.residue2_.getNumberAndInsertionCode() << "\n";
    }
    std::cout << "Chain terminations:\n";
    for (auto& chainT : ppInfo.chainTerminals_)
    {
        std::cout << chainT.chainId_ << " | " << chainT.startIndex_ << " | " << chainT.nTermination_ << " | "
                  << chainT.endIndex_ << " | " << chainT.cTermination_ << "\n";
    }
    std::cout << "NonNatural Protein Residues:\n";
    for (auto& nonNaturalResidue : ppInfo.nonNaturalProteinResidues_)
    {
        std::cout << nonNaturalResidue.residue_.getChainId() << " | " << nonNaturalResidue.residue_.getName() << " | "
                  << nonNaturalResidue.residue_.getNumberAndInsertionCode() << "\n";
    }
    return 0;
}
