#include "include/CentralDataStructure/pdbWriter.hpp"

#include "include/CentralDataStructure/cdsFunctions.hpp"
#include "include/CentralDataStructure/graphInterface.hpp"
#include "include/assembly/assemblyGraph.hpp"
#include "include/assembly/assemblyIndices.hpp"
#include "include/assembly/assemblyTypes.hpp"
#include "include/fileType/pdb/pdbFileData.hpp"
#include "include/fileType/pdb/pdbFileWriter.hpp"
#include "include/fileType/pdb/pdbFunctions.hpp"
#include "include/graph/graphManipulation.hpp"
#include "include/graph/graphTypes.hpp"
#include "include/util/containers.hpp"
#include "include/util/strings.hpp"

#include <iomanip>
#include <ostream>
#include <string>
#include <vector>

namespace gmml
{
    pdb::PdbFileAtomData toPdbFileAtomData(const std::vector<Atom*>& atoms, std::vector<std::string> recordNames)
    {
        return {
            atomCoordinates(atoms),
            atomNumbers(atoms),
            atomNames(atoms),
            atomElementStrings(atoms),
            recordNames,
            std::vector<double>(atoms.size(), 1.0),
            std::vector<double>(atoms.size(), 0.0)};
    }

    pdb::PdbFileData toPdbFileData(const GraphObjects& objects)
    {
        const std::vector<Atom*>& atoms = objects.atoms;
        const std::vector<Residue*>& residues = objects.residues;
        std::vector<std::string> recordNames(atoms.size(), "ATOM");
        std::vector<std::string> chainIds(residues.size(), "");
        std::vector<std::string> insertionCodes(residues.size(), "");
        pdb::PdbFileResidueData residueData {
            residueNumbers(residues), truncatedResidueNames(residues), chainIds, insertionCodes};
        pdb::PdbFileFormat format;
        return pdb::PdbFileData {format, {}, residueData, toPdbFileAtomData(atoms, recordNames)};
    }

    void writeTrajectoryToPdb(
        std::ostream& stream, const pdb::PdbData& data, const std::vector<size_t>& selectedResidues)
    {
        std::function<bool(const size_t&)> atomAlive = [&](size_t n) { return data.indices.atomAlive[n]; };
        std::function<std::string(const size_t&)> elementString = [&](size_t n)
        { return data.atoms.names[n].empty() ? "" : data.atoms.names[n].substr(0, 1); };
        std::function<std::string(const size_t&)> truncatedResidueNames = [&](size_t n)
        { return util::truncate(3, data.residues.names[n]); };

        std::vector<bool> residueIncluded = util::indicesToBools(data.indices.residueCount, selectedResidues);
        assembly::Graph graph = createAssemblyGraph(data.indices, data.atomGraph);
        pdb::PdbFileResidueData residueData {
            data.residues.numbers,
            util::vectorMap(truncatedResidueNames, util::indexVector(data.indices.residueCount)),
            std::vector<std::string>(data.indices.residueCount, ""),
            data.residues.insertionCodes};
        std::vector<std::string> elements = util::vectorMap(elementString, util::indexVector(data.indices.atomCount));
        pdb::PdbFileFormat format;

        size_t modelCount = data.trajectory.coordinates.size();
        for (size_t coordinateSet = 0; coordinateSet < modelCount; coordinateSet++)
        {
            pdb::PdbFileAtomData atomData {
                data.trajectory.coordinates[coordinateSet],
                data.atoms.numbers,
                data.atoms.names,
                elements,
                data.atoms.recordNames,
                data.atoms.occupancies,
                data.atoms.temperatureFactors};
            pdb::PdbFileData writerData {format, {}, residueData, atomData};
            stream << "MODEL " << std::right << std::setw(8) << (coordinateSet + 1) << "\n";
            for (size_t moleculeId = 0; moleculeId < data.indices.moleculeCount; moleculeId++)
            {
                std::vector<size_t> residueIds =
                    util::boolsToIndices(util::vectorAnd(residueIncluded, isMoleculeResidue(data.indices, moleculeId)));
                if (residueIds.size() > 0)
                {
                    std::vector<ResidueType> types = util::indicesToValues(data.residues.types, residueIds);
                    std::vector<bool> ter = pdb::residueTER(types);
                    pdb::writeMoleculeToPdb(stream, graph, residueIds, ter, writerData);
                }
            }
            stream << "ENDMDL\n";
        }
        pdb::theEnd(stream);
    }

    void WritePdb(std::ostream& stream, const GraphIndexData& graphData, const std::vector<std::string>& headerLines)
    {
        std::vector<bool> includedAtoms = atomVisibility(graphData.objects.atoms);
        graph::Database atomGraphData = createGraphData(graphData.objects);
        atomGraphData.nodeAlive = atomVisibility(graphData.objects.atoms);
        assembly::Graph graph = createAssemblyGraph(graphData.indices, atomGraphData);
        const std::vector<Residue*>& residues = graphData.objects.residues;
        std::vector<ResidueType> types = residueTypes(residues);
        std::vector<bool> ter = pdb::residueTER(types);
        pdb::PdbFileData data = toPdbFileData(graphData.objects);
        for (auto& line : headerLines)
        {
            stream << "HEADER    " << line << "\n";
        }
        pdb::writeMoleculeToPdb(stream, graph, util::indexVector(residues), ter, data);
        std::vector<bool> selectedResidueTypes(ResidueTypeCount, false);
        for (auto type : {Sugar, Derivative, Aglycone, Undefined})
        {
            selectedResidueTypes[type] = true;
        }
        std::vector<bool> residueSelected = util::indicesToValues(selectedResidueTypes, types);
        std::vector<bool> atomSelected = util::indicesToValues(residueSelected, graphData.indices.atomResidue);
        graph::Database subgraphData = createGraphData(graphData.objects);
        subgraphData.nodeAlive = util::vectorAnd(includedAtoms, atomSelected);
        assembly::Graph subgraph = createAssemblyGraph(graphData.indices, subgraphData);
        std::vector<std::array<size_t, 2>> connectionIndices =
            util::indicesToValues(subgraph.residues.source.edgeNodes, subgraph.residues.edges.indices);
        pdb::writeConectCards(stream, data.atoms.numbers, connectionIndices);
        pdb::theEnd(stream);
    }
} // namespace gmml
