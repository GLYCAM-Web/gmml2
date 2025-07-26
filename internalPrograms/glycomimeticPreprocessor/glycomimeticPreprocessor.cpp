#include "include/CentralDataStructure/graphInterface.hpp"
#include "include/assembly/assemblyGraph.hpp"
#include "include/assembly/assemblyIndices.hpp"
#include "include/assembly/assemblyTypes.hpp"
#include "include/pdb/pdbAtom.hpp"
#include "include/pdb/pdbFile.hpp"
#include "include/pdb/pdbPreprocessorInputs.hpp"
#include "include/pdb/pdbResidue.hpp"
#include "include/readers/parameterManager.hpp"
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
    pdb::PdbFile pdbFile(argv[1], {pdb::InputType::modelsAsMolecules, false});
    auto panic = [&]()
    {
        for (size_t n = 0; n < pdbFile.data.indices.residueCount; n++)
        {
            auto atomIds = residueAtoms(pdbFile.data.indices, n);
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
    pdb::PreprocessorOptions options; // Default values are good.
    std::cout << "Preprocessing\n";
    ParameterManager parameterManager = loadParameters(baseDir);
    parameterManager.lib.residueNames = util::reverse(parameterManager.lib.residueNames);
    parameterManager.lib.residues = util::reverse(parameterManager.lib.residues);
    pdb::PreprocessorInformation ppInfo = pdbFile.PreProcess(parameterManager, options);
    std::vector<Assembly*> assemblies = pdbFile.getAssemblies();
    assembly::Indices& graphIndices = pdbFile.data.indices;
    std::vector<size_t> moleculeIds = util::indicesOfElement(graphIndices.moleculeAssembly, size_t(0));
    size_t residueCount = pdbFile.data.indices.residueCount;
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
    std::vector<size_t> firstAssemblyMoleculeIds =
        util::indicesOfElement(pdbFile.data.indices.moleculeAssembly, size_t(0));
    std::function<std::vector<size_t>(const std::vector<size_t>&)> onlyAliveResidues =
        [&](const std::vector<size_t>& residueIds)
    {
        std::function<bool(const size_t&)> isAlive = [&](size_t id) { return residueAlive[id]; };
        return util::vectorFilter(isAlive, residueIds);
    };
    std::vector<std::vector<size_t>> residueOrder = util::vectorMap(
        onlyAliveResidues, util::indicesToValues(pdbFile.data.molecules.residueOrder, firstAssemblyMoleculeIds));
    pdb::Write(pdbFile.data, residueOrder, outFileStream);
    outFileStream << "END\n"; // Original GMML needs this.
    outFileStream.close();

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
