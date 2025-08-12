#include "include/carbohydrate/pdbWriter.hpp"

#include "include/assembly/assemblyIndices.hpp"
#include "include/assembly/assemblyTypes.hpp"
#include "include/fileType/pdb/pdbFileData.hpp"
#include "include/fileType/pdb/pdbFileWriter.hpp"
#include "include/fileType/pdb/pdbFunctions.hpp"
#include "include/graph/graphManipulation.hpp"
#include "include/graph/graphTypes.hpp"
#include "include/util/containers.hpp"
#include "include/util/strings.hpp"

#include <ostream>
#include <string>
#include <vector>

namespace gmml
{
    void writePdb(
        std::ostream& stream,
        const carbohydrate::CarbohydrateData& data,
        const assembly::Graph& graph,
        const std::vector<std::string>& headerLines)
    {
        size_t atomCount = data.indices.atomCount;
        size_t residueCount = data.indices.residueCount;
        std::function<std::string(const std::string&)> truncate = [](const std::string& str)
        { return util::truncate(3, str); };

        std::vector<bool> ter = pdb::residueTER(data.residues.types);

        pdb::PdbFileAtomData atomData {
            data.atoms.coordinates,
            data.atoms.numbers,
            data.atoms.names,
            elementNames(data.atoms.elements),
            std::vector<std::string>(atomCount, "ATOM"),
            std::vector<double>(atomCount, 1.0),
            std::vector<double>(atomCount, 0.0)};

        std::vector<std::string> chainIds(residueCount, "");
        std::vector<std::string> insertionCodes(residueCount, "");
        pdb::PdbFileResidueData residueData {
            data.residues.numbers, util::vectorMap(truncate, data.residues.names), chainIds, insertionCodes};
        pdb::PdbFileFormat format;
        pdb::PdbFileData pdbData {format, {}, residueData, atomData};

        for (auto& line : headerLines)
        {
            stream << "HEADER    " << line << "\n";
        }
        pdb::writeMoleculeToPdb(stream, pdbData, graph, util::indexVector(residueCount), ter);
        std::vector<std::array<size_t, 2>> connectionIndices =
            util::indicesToValues(graph.residues.source.edges.nodes, graph.residues.edges.sourceIndices);
        pdb::writeConectCards(stream, data.atoms.numbers, connectionIndices);
        pdb::theEnd(stream);
    }
} // namespace gmml
