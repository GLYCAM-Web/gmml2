#include "include/carbohydrate/offWriter.hpp"

#include "include/assembly/assemblyGraph.hpp"
#include "include/fileType/off/offFileData.hpp"
#include "include/fileType/off/offFileWriter.hpp"
#include "include/graph/graphFunctions.hpp"
#include "include/util/containers.hpp"

#include <string>
#include <vector>

namespace gmml
{
    std::vector<std::vector<size_t>> atomsConnectedToOtherResidues(const assembly::Graph& graph)
    {
        const graph::Database& atomDB = graph.atoms.source;
        const std::vector<size_t>& atomResidue = graph.indices.atomResidue;
        std::vector<std::vector<size_t>> connections(graph.indices.residueCount, std::vector<size_t> {});
        for (size_t n : atomDB.edges.indices)
        {
            if (graph::edgeAlive(atomDB, n))
            {
                std::array<size_t, 2> nodes = atomDB.edges.nodes[n];
                std::array<size_t, 2> res = {atomResidue[nodes[0]], atomResidue[nodes[1]]};
                if (res[0] != res[1])
                {
                    for (size_t k = 0; k < 2; k++)
                    {
                        connections[res[k]].push_back(nodes[k]);
                    }
                }
            }
        }
        return connections;
    }

    off::OffFileData toOffFileData(const carbohydrate::CarbohydrateData& data, const assembly::Graph& graph)
    {
        std::vector<std::vector<size_t>> connections = atomsConnectedToOtherResidues(graph);

        off::OffFileAtomData atomData {
            util::serializedNumberVector(data.atoms.visible),
            data.atoms.names,
            data.atoms.types,
            data.atoms.atomicNumbers,
            data.atoms.charges,
            data.atoms.coordinates};
        off::OffFileResidueData residueData {
            util::serializedNumberVector(graph.indices.residueCount),
            data.residues.names,
            data.residues.types,
            connections};
        off::OffFileFormat format;
        return off::OffFileData {format, residueData, atomData};
    }
} // namespace gmml
