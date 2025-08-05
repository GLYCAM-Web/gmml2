#include "include/fileType/prep/prepFunctions.hpp"

#include "include/graph/graphFunctions.hpp"
#include "include/graph/graphManipulation.hpp"
#include "include/util/containers.hpp"

#include <functional>
#include <string>
#include <vector>

namespace gmml
{
    namespace prep
    {
        std::vector<size_t> residueAtoms(const PrepData& data, size_t residueId)
        {
            std::vector<size_t> atoms = util::indicesOfElement(data.atomResidue, residueId);
            std::vector<bool> alive = util::indicesToValues(data.atomGraph.nodeAlive, atoms);
            return util::boolsToValues(atoms, alive);
        }

        std::vector<size_t> residueEdges(const PrepData& data, size_t residueId)
        {
            auto partOfResidue = [&](size_t n) { return data.atomResidue[n] == residueId; };
            std::function<bool(const size_t&)> condition = [&](size_t edgeId)
            {
                const std::array<size_t, 2>& nodes = data.atomGraph.edgeNodes[edgeId];
                return partOfResidue(nodes[0]) && partOfResidue(nodes[1]) && graph::edgeAlive(data.atomGraph, edgeId);
            };
            return util::vectorFilter(condition, data.atomGraph.edges);
        }

        std::vector<size_t> atomNeighbors(const PrepData& data, size_t atomId)
        {
            std::vector<size_t> result;
            size_t residue = data.atomResidue[atomId];
            std::vector<size_t> residueEdges = prep::residueEdges(data, residue);
            for (size_t edgeId : prep::residueEdges(data, residue))
            {
                const std::array<size_t, 2>& nodes = data.atomGraph.edgeNodes[edgeId];
                if (nodes[0] == atomId)
                {
                    result.push_back(nodes[1]);
                }
                else if (nodes[1] == atomId)
                {
                    result.push_back(nodes[0]);
                }
            }
            return result;
        }

        size_t findAtom(const PrepData& data, size_t residueId, const std::string& name)
        {
            std::function<bool(const size_t&)> condition = [&](size_t id) {
                return data.atomResidue[id] == residueId && data.atoms.name[id] == name && data.atomGraph.nodeAlive[id];
            };
            return util::findIf(data.atomGraph.nodes, condition);
        }

        size_t addAtom(prep::PrepData& data, size_t residueId)
        {
            prep::AtomData& atoms = data.atoms;
            size_t atomId = graph::addNode(data.atomGraph);
            data.atomResidue.push_back(residueId);
            data.atomCount++;
            atoms.number.push_back(0);
            atoms.name.push_back("");
            atoms.type.push_back("");
            atoms.topologicalType.push_back(prep::kTopTypeE);
            atoms.bondIndex.push_back(0);
            atoms.angleIndex.push_back(0);
            atoms.dihedralIndex.push_back(0);
            atoms.bondLength.push_back(0);
            atoms.angle.push_back(0);
            atoms.dihedral.push_back(0);
            atoms.charge.push_back(0);
            atoms.coordinate.push_back({0, 0, 0});
            return atomId;
        }

        size_t copyResidue(prep::PrepData& data, const prep::PrepData& reference, size_t residue)
        {
            prep::AtomData& atoms = data.atoms;
            prep::ResidueData& residues = data.residues;
            size_t newResidue = data.residueCount;
            data.residueCount++;
            residues.name.push_back(reference.residues.name[residue]);
            residues.title.push_back(reference.residues.title[residue]);
            residues.coordinateType.push_back(reference.residues.coordinateType[residue]);
            residues.outputFormat.push_back(reference.residues.outputFormat[residue]);
            residues.geometryType.push_back(reference.residues.geometryType[residue]);
            residues.dummyAtomOmission.push_back(reference.residues.dummyAtomOmission[residue]);
            residues.dummyAtomType.push_back(reference.residues.dummyAtomType[residue]);
            residues.dummyAtomPosition.push_back(reference.residues.dummyAtomPosition[residue]);
            residues.charge.push_back(reference.residues.charge[residue]);
            residues.improperDihedrals.push_back(reference.residues.improperDihedrals[residue]);
            residues.loops.push_back(reference.residues.loops[residue]);

            std::vector<size_t> newAtomIndex = reference.atomGraph.nodes;

            for (size_t atom : prep::residueAtoms(reference, residue))
            {
                size_t newAtom = graph::addNode(data.atomGraph);
                newAtomIndex[atom] = newAtom;
                data.atomResidue.push_back(newResidue);
                data.atomCount++;
                atoms.number.push_back(reference.atoms.number[atom]);
                atoms.name.push_back(reference.atoms.name[atom]);
                atoms.type.push_back(reference.atoms.type[atom]);
                atoms.topologicalType.push_back(reference.atoms.topologicalType[atom]);
                atoms.bondIndex.push_back(reference.atoms.bondIndex[atom]);
                atoms.angleIndex.push_back(reference.atoms.angleIndex[atom]);
                atoms.dihedralIndex.push_back(reference.atoms.dihedralIndex[atom]);
                atoms.bondLength.push_back(reference.atoms.bondLength[atom]);
                atoms.angle.push_back(reference.atoms.angle[atom]);
                atoms.dihedral.push_back(reference.atoms.dihedral[atom]);
                atoms.charge.push_back(reference.atoms.charge[atom]);
                atoms.coordinate.push_back(reference.atoms.coordinate[atom]);
            }

            for (size_t edgeId : prep::residueEdges(reference, residue))
            {
                const std::array<size_t, 2>& nodes = reference.atomGraph.edgeNodes[edgeId];
                graph::addEdge(data.atomGraph, {newAtomIndex[nodes[0]], newAtomIndex[nodes[1]]});
            }
            return newResidue;
        }
    } // namespace prep
} // namespace gmml
