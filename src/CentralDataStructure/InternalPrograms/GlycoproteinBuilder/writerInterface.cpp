#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/writerInterface.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycoproteinStructs.hpp"
#include "includes/CentralDataStructure/FileFormats/offFileData.hpp"
#include "includes/CentralDataStructure/FileFormats/pdbFileData.hpp"
#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
#include "includes/Assembly/assemblyGraph.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/version.h"

#include <vector>

namespace glycoproteinBuilder
{
    cds::OffFileData toOffFileData(const assembly::Graph& graph, const AssemblyData& data,
                                   const std::vector<cds::Coordinate>& atomCoordinates)
    {
        cds::OffFileAtomData atomData {data.atoms.serializedNumbers, data.atoms.names,   data.atoms.types,
                                       data.atoms.atomicNumbers,     data.atoms.charges, atomCoordinates};

        size_t residueCount = graph.indices.residueCount;
        std::vector<std::vector<size_t>> atomsConnectedToOtherResidues;
        atomsConnectedToOtherResidues.resize(residueCount);
        for (size_t n = 0; n < edgeCount(graph.residues); n++)
        {
            std::array<size_t, 2> residueAdj = graph.residues.edges.nodeAdjacencies[n];
            size_t atomEdge                  = residueEdgeToAtomEdgeIndex(graph, n);
            std::array<size_t, 2> atomAdj    = graph.atoms.edges.nodeAdjacencies[atomEdge];
            for (size_t k = 0; k < 2; k++)
            {
                atomsConnectedToOtherResidues[residueAdj[k]].push_back(atomAdj[k]);
            }
        }
        // sorting more closely aligns with cds behavior in print, but perhaps unnecessary
        for (auto& vec : atomsConnectedToOtherResidues)
        {
            std::sort(vec.begin(), vec.end());
        }
        cds::OffFileResidueData residueData {data.residues.serializedNumbers, data.residues.names, data.residues.types,
                                             atomsConnectedToOtherResidues};
        cds::OffFileFormat format;
        format.coordinate.decimals.precision = 5;
        return cds::OffFileData {format, residueData, atomData};
    }

    cds::PdbFileData toPdbFileData(const assembly::Indices& indices, const AssemblyData& data,
                                   const std::vector<cds::Coordinate>& atomCoordinates,
                                   const std::vector<uint>& atomNumbers, const std::vector<uint>& residueNumbers,
                                   const std::vector<std::string>& chainIds,
                                   const std::vector<std::string>& headerLines)
    {
        size_t atomCount    = indices.atomCount;
        size_t residueCount = indices.residueCount;
        std::vector<std::string> recordNames(atomCount, "ATOM");
        std::vector<std::string> insertionCodes(residueCount, "");

        cds::PdbFileResidueData residuePdbData {residueNumbers, data.residues.names, chainIds, insertionCodes};
        cds::PdbFileAtomData atomPdbData {atomCoordinates,
                                          atomNumbers,
                                          data.atoms.names,
                                          data.atoms.elementStrings,
                                          recordNames,
                                          std::vector<double>(atomCount, 1.0),
                                          std::vector<double>(atomCount, 0.0)};
        cds::PdbFileFormat format;
        return cds::PdbFileData {format, headerLines, residuePdbData, atomPdbData};
    }
} // namespace glycoproteinBuilder
