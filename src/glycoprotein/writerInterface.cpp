#include "include/glycoprotein/writerInterface.hpp"

#include "include/assembly/assemblyGraph.hpp"
#include "include/assembly/assemblyIndices.hpp"
#include "include/fileType/off/offFileData.hpp"
#include "include/fileType/pdb/pdbFileData.hpp"
#include "include/geometry/geometryTypes.hpp"
#include "include/glycoprotein/glycoproteinStructs.hpp"
#include "include/util/containers.hpp"
#include "include/version.h"

#include <vector>

namespace gmml
{
    namespace gpbuilder
    {
        off::OffFileData toOffFileData(
            const assembly::Graph& graph, const AssemblyData& data, const std::vector<Coordinate>& atomCoordinates)
        {
            off::OffFileAtomData atomData {
                data.atoms.serializedNumbers,
                data.atoms.names,
                data.atoms.types,
                data.atoms.atomicNumbers,
                data.atoms.charges,
                atomCoordinates};

            size_t residueCount = assembly::residueCount(graph.source);
            std::vector<std::vector<size_t>> atomsConnectedToOtherResidues;
            atomsConnectedToOtherResidues.resize(residueCount);
            for (size_t n = 0; n < edgeCount(graph.residues); n++)
            {
                std::array<size_t, 2> residueAdj = graph.residues.edges.nodeAdjacencies[n];
                size_t atomEdge = residueEdgeToAtomEdgeIndex(graph, n);
                std::array<size_t, 2> atomAdj = graph.atoms.edges.nodeAdjacencies[atomEdge];
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
            off::OffFileResidueData residueData {
                data.residues.serializedNumbers,
                data.residues.names,
                data.residues.types,
                atomsConnectedToOtherResidues};
            off::OffFileFormat format;
            format.coordinate.decimals.precision = 5;
            return off::OffFileData {format, residueData, atomData};
        }

        pdb::PdbFileData toPdbFileData(
            const assembly::Indices& indices,
            const AssemblyData& data,
            const std::vector<Coordinate>& atomCoordinates,
            const std::vector<uint>& atomNumbers,
            const std::vector<uint>& residueNumbers,
            const std::vector<std::string>& chainIds,
            const std::vector<std::string>& headerLines)
        {
            size_t atomCount = indices.atomCount;
            size_t residueCount = indices.residueCount;
            std::vector<std::string> recordNames(atomCount, "ATOM");
            std::vector<std::string> insertionCodes(residueCount, "");

            pdb::PdbFileResidueData residuePdbData {residueNumbers, data.residues.names, chainIds, insertionCodes};
            pdb::PdbFileAtomData atomPdbData {
                atomCoordinates,
                atomNumbers,
                data.atoms.names,
                data.atoms.elementStrings,
                recordNames,
                std::vector<double>(atomCount, 1.0),
                std::vector<double>(atomCount, 0.0)};
            pdb::PdbFileFormat format;
            return pdb::PdbFileData {format, headerLines, residuePdbData, atomPdbData};
        }
    } // namespace gpbuilder
} // namespace gmml
