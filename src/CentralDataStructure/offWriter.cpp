#include "include/CentralDataStructure/offWriter.hpp"

#include "include/CentralDataStructure/atom.hpp"
#include "include/CentralDataStructure/cdsFunctions.hpp"
#include "include/CentralDataStructure/graphInterface.hpp"
#include "include/CentralDataStructure/residue.hpp"
#include "include/assembly/assemblyGraph.hpp"
#include "include/fileType/off/offFileData.hpp"
#include "include/fileType/off/offFileWriter.hpp"
#include "include/util/containers.hpp"

#include <string>
#include <vector>

namespace gmml
{
    namespace
    {
        std::vector<size_t> atomIndices(const std::vector<Atom*>& atoms, const std::vector<Atom*>& reindex)
        {
            std::vector<size_t> result;
            result.reserve(reindex.size());
            for (auto& find : reindex)
            {
                result.push_back(atomVectorIndex(atoms, find));
            }
            return result;
        }

        std::vector<Atom*> atomsConnectedToOtherResidues(std::vector<Atom*> atoms)
        {
            std::vector<Atom*> foundAtoms;
            for (auto& atom : atoms)
            {
                bool connected = false;
                for (auto& neighbor : atom->getNeighbors())
                { // check if neighbor is not one of the atoms in this residue.
                    connected = connected || !util::contains(atoms, neighbor);
                }
                if (connected)
                {
                    foundAtoms.push_back(atom);
                }
            }
            return foundAtoms;
        }
    } // namespace

    off::OffFileData toOffFileData(const std::vector<Residue*>& residues)
    {
        std::vector<Atom*> atoms;
        std::vector<std::vector<size_t>> connections;
        for (auto& residue : residues)
        {
            std::vector<Atom*> residueAtoms = residue->getAtoms();
            std::vector<Atom*> atomsConnected = atomsConnectedToOtherResidues(residueAtoms);
            util::insertInto(atoms, residueAtoms);
            connections.push_back(atomIndices(atoms, atomsConnected));
        }
        off::OffFileAtomData atomData {
            atomNumbers(atoms),
            atomNames(atoms),
            atomTypes(atoms),
            atomAtomicNumbers(atoms),
            atomCharges(atoms),
            atomCoordinates(atoms)};
        off::OffFileResidueData residueData {
            residueNumbers(residues), residueNames(residues), residueTypes(residues), connections};
        off::OffFileFormat format;
        return off::OffFileData {format, residueData, atomData};
    }

    void serializeResiduesIndividually(std::vector<Residue*>& residues)
    {
        for (auto& residue : residues)
        {
            residue->setNumber(1);
            serializeNumbers(residue->getAtoms());
        }
    }

    void WriteOff(std::ostream& stream, const std::string& name, const GraphIndexData& data)
    {
        assembly::Graph graph = createVisibleAssemblyGraph(data);
        off::OffFileData offData = toOffFileData(data.objects.residues);
        offData.atoms.numbers = serializedNumberVector(graph.atoms.source.nodes.alive);
        offData.residues.numbers = serializedNumberVector(data.indices.residueCount);
        off::writeResiduesTogether(stream, graph, offData, util::indexVector(data.objects.residues), name);
    }
} // namespace gmml
