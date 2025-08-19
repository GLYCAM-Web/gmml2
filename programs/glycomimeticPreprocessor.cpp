#include "include/CentralDataStructure/graphInterface.hpp"
#include "include/assembly/assemblyGraph.hpp"
#include "include/assembly/assemblyIndices.hpp"
#include "include/assembly/assemblyTypes.hpp"
#include "include/fileType/pdb/pdbAtom.hpp"
#include "include/fileType/pdb/pdbFile.hpp"
#include "include/fileType/pdb/pdbModel.hpp"
#include "include/fileType/pdb/pdbResidue.hpp"
#include "include/preprocess/parameterManager.hpp"
#include "include/preprocess/pdbPreProcess.hpp"
#include "include/preprocess/pdbPreprocessorInputs.hpp"
#include "include/util/containers.hpp"
#include "include/util/filesystem.hpp"

#include <fstream>
#include <iostream>
#include <string>

// Adapted the preprocessor to prepare files for Yao's glycomimetic program.
// In addition to regular preprocessing it:
// Removes water
// Renames HETATM entries to ATOM. // I think this will be automatic anyway.
// Only writes out the first MODEL if there are many.

int main(int argc, char* argv[])
{
    using namespace gmml;
    if (argc != 3)
    {
        std::cout << "Usage: " << argv[0] << " inputFile.pdb outfileName\n";
        std::cout << "Example: " << argv[0] << " inputs/4mbz.pdb 030.outputPdbFile.pdb\n";
        std::exit(EXIT_FAILURE);
    }
    pdb::PdbFile pdbFile = pdb::toPdbFile(argv[1], pdb::modelsAsMolecules);
    auto panic = [&]()
    {
        for (size_t n = 0; n < residueCount(pdbFile.data.assembly); n++)
        {
            auto atomIds = residueAtoms(pdbFile.data.assembly.indices, n);
            auto atoms = pdbFile.data.objects.residues[n]->getAtoms();
            if (atomIds.size() != atoms.size())
            {
                std::cout << "residue " << n << "\n";
                std::cout << atomIds.size() << " vs " << atoms.size() << "\n";
                throw std::runtime_error("panic");
            }
        }
    };
    std::string baseDir = util::toString(util::pathAboveCurrentExecutableDir());
    preprocess::PreprocessorOptions options = preprocess::defaultPreprocessorOptions; // Default values are good.
    std::cout << "Preprocessing\n";
    preprocess::ParameterManager parameterManager = preprocess::loadParameters(baseDir);
    parameterManager.lib.residueNames = util::reverse(parameterManager.lib.residueNames);
    parameterManager.lib.residues = util::reverse(parameterManager.lib.residues);
    preprocess::PreprocessorInformation ppInfo = preProcess(pdbFile, parameterManager, options);
    std::vector<Assembly*> assemblies = getAssemblies(pdbFile);
    assembly::Indices& graphIndices = pdbFile.data.assembly.indices;
    std::vector<size_t> moleculeIds = util::indicesOfElement(graphIndices.moleculeAssembly, size_t(0));
    size_t residueCount = assembly::residueCount(pdbFile.data.assembly);
    std::vector<bool> residueAlive(residueCount, true);
    for (size_t residueId = 0; residueId < residueCount; residueId++)
    {
        Residue* residue = pdbFile.data.objects.residues[residueId];
        if (residue->GetType() != ResidueType::Protein)
        {
            std::vector<size_t> atomIds = residueAtoms(graphIndices, residueId);
            for (size_t atomId : atomIds)
            {
                pdbFile.data.atoms.recordNames[atomId] = "ATOM";
            }
        }
        if (residue->getName() == "HOH" || residue->getName() == "WAT")
        {
            std::cout << "Deleting " << residue->getName() << "\n";
            residueAlive[residueId] = false;
        }
    }
    std::ofstream outFileStream;
    outFileStream.open(argv[2]);
    std::vector<size_t> firstAssemblyMoleculeIds = util::indicesOfElement(graphIndices.moleculeAssembly, size_t(0));
    std::function<std::vector<size_t>(const std::vector<size_t>&)> onlyAliveResidues =
        [&](const std::vector<size_t>& residueIds)
    {
        std::function<bool(const size_t&)> isAlive = [&](size_t id) { return residueAlive[id]; };
        return util::vectorFilter(isAlive, residueIds);
    };
    std::vector<std::vector<size_t>> residueOrder = util::vectorMap(
        onlyAliveResidues, util::indicesToValues(pdbFile.data.molecules.residueOrder, firstAssemblyMoleculeIds));
    pdb::write(pdbFile.data, residueOrder, outFileStream);
    outFileStream << "END\n"; // Original GMML needs this.
    outFileStream.close();

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
