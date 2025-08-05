#include "include/CentralDataStructure/pdbWriter.hpp"
#include "include/assembly/assemblyIndices.hpp"
#include "include/assembly/assemblyTypes.hpp"
#include "include/fileType/pdb/pdbFile.hpp"
#include "include/fileType/pdb/pdbFunctions.hpp"
#include "include/fileType/pdb/pdbResidue.hpp"
#include "include/geometry/geometryFunctions.hpp"
#include "include/util/containers.hpp"
#include "include/util/files.hpp"

#include <iostream>
#include <ostream>
#include <string>
#include <vector>

int main(int argc, char* argv[])
{
    using namespace gmml;
    if (argc != 2)
    {
        std::cout << "Usage: " << argv[0] << " inputFile.pdb\n";
        std::cout << "Example: " << argv[0] << " inputs/4mbz.pdb\n";
        std::exit(EXIT_FAILURE);
    }
    // requirement: a chain ID for every single ATOM entry, and all ligand atoms should be put in a single residue.
    pdb::PdbFile pdbFile = pdb::toPdbFile(argv[1], pdb::modelsAsMolecules);
    // PdbFile is an "Ensemble" (made up of "Assemblies"), but if you want to just set
    // every molecule to have any chain ID you can do:
    pdbFile.data.residues.chainIds = std::vector<std::string>(pdbFile.data.indices.residueCount, "Y");
    // ResidueTypes are guessed upon input. Using that guess to find the ligand, can improve this if you need:
    std::vector<size_t> ligandResidues = util::indicesOfElement(pdbFile.data.residues.types, ResidueType::Undefined);
    if (ligandResidues.empty())
    {
        std::cout << "No ligand residues found in input file\n";
        return 0;
    }
    size_t firstLigandResidue = ligandResidues.front();
    const std::string& firstName = pdbFile.data.residues.names[firstLigandResidue];
    uint firstNumber = pdbFile.data.residues.numbers[firstLigandResidue];
    for (size_t ligandResidue : ligandResidues) // Each MODEL in PdbFile is converted into an "Assembly"
    {                                           // Every ligand residue gets the same residue number as the first one.
        pdbFile.data.residues.names[ligandResidue] = firstName;
        pdbFile.data.residues.numbers[ligandResidue] = firstNumber;
    }
    write(pdbFile, "026.outputPdbFile.pdb");
    // ************************************************************************ //
    // Separate thing showing how to read/write PDB files as "trajectories/frames"
    pdb::PdbFile pdbFileTraj = pdb::toPdbFile(argv[1], {pdb::InputType::modelsAsCoordinates, false});

    std::function<Coordinate(const size_t&)> residueMean = [&](size_t n)
    {
        std::vector<size_t> atomIds = residueAtoms(pdbFileTraj.data.indices, n);
        return coordinateMean(util::indicesToValues(pdbFileTraj.data.atoms.coordinates, atomIds));
    };
    std::vector<Coordinate> residueMeans =
        util::vectorMap(residueMean, util::indexVector(pdbFileTraj.data.indices.residueCount));
    // somehow you specify number in inputs. e.g. A_405 chain A, residue 405.
    size_t residueIndex = util::indexOf(pdbFileTraj.data.residues.numbers, uint(5));
    double distance = 12.345; // inputs

    std::vector<size_t> selectedResidues;
    for (size_t n = 0; n < pdbFileTraj.data.indices.residueCount; n++)
    {
        if (withinDistance(distance, residueMeans[residueIndex], residueMeans[n]))
        {
            selectedResidues.push_back(n);
        }
    }

    std::cout << "Found " << selectedResidues.size() << " residues\n";
    pdbFileTraj.data.indices.residueMolecule = std::vector<size_t>(pdbFileTraj.data.indices.residueCount, 0);
    util::writeToFile(
        "026.outputSelection.pdb",
        [&](std::ostream& stream) { writeTrajectoryToPdb(stream, pdbFileTraj.data, selectedResidues); });
    return 0;
}
