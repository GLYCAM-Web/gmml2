#include "includes/CentralDataStructure/Writers/pdbWriter.hpp"
#include "includes/CentralDataStructure/FileFormats/pdbFileData.hpp"
#include "includes/CentralDataStructure/FileFormats/pdbFileWriter.hpp"
#include "includes/CentralDataStructure/cdsFunctions/cdsFunctions.hpp"
#include "includes/CentralDataStructure/cdsFunctions/graphInterface.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbFunctions.hpp"
#include "includes/Assembly/assemblyGraph.hpp"
#include "includes/Graph/graphTypes.hpp"
#include "includes/Graph/graphManipulation.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/strings.hpp"

#include <iomanip>
#include <ostream>
#include <string>
#include <vector>

std::vector<bool> cds::residueTER(const std::vector<ResidueType>& types)
{
    std::vector<bool> result;
    result.reserve(types.size());
    for (size_t n = 0; n < types.size(); n++)
    {
        const ResidueType& type = types[n];
        size_t next             = n + 1;
        bool isSugar            = type == cds::ResidueType::Undefined || type == cds::ResidueType::Sugar ||
                       type == cds::ResidueType::Derivative || type == cds::ResidueType::Aglycone;
        bool nextIsCapping           = next < types.size() && types[next] == cds::ResidueType::ProteinCappingGroup;
        bool betweenTwoCappingGroups = (type == cds::ResidueType::ProteinCappingGroup) && nextIsCapping;
        bool isLast                  = next == types.size();
        bool isProtein               = type == cds::ResidueType::Protein;
        result.push_back(isSugar || betweenTwoCappingGroups || (isLast && isProtein));
    }
    return result;
}

cds::PdbFileAtomData cds::toPdbFileAtomData(const std::vector<cds::Atom*>& atoms, std::vector<std::string> recordNames)
{
    return {atomCoordinates(atoms),
            atomNumbers(atoms),
            atomNames(atoms),
            atomElementStrings(atoms),
            recordNames,
            std::vector<double>(atoms.size(), 1.0),
            std::vector<double>(atoms.size(), 0.0)};
}

cds::PdbFileData cds::toPdbFileData(const cds::GraphObjects& objects)
{
    const std::vector<cds::Atom*>& atoms       = objects.atoms;
    const std::vector<cds::Residue*>& residues = objects.residues;
    std::vector<std::string> recordNames(atoms.size(), "ATOM");
    std::vector<std::string> chainIds(residues.size(), "");
    std::vector<std::string> insertionCodes(residues.size(), "");
    PdbFileResidueData residueData {residueNumbers(residues), truncatedResidueNames(residues), chainIds,
                                    insertionCodes};
    PdbFileFormat format;
    return PdbFileData {format, {}, residueData, toPdbFileAtomData(atoms, recordNames)};
}

void cds::writeTrajectoryToPdb(std::ostream& stream, const pdb::PdbData& data,
                               const std::vector<size_t>& selectedResidues)
{
    std::function<bool(const size_t&)> atomAlive = [&](size_t n)
    {
        return data.atomGraph.nodeAlive[n];
    };
    std::function<std::string(const size_t&)> elementString = [&](size_t n)
    {
        return data.atoms.names[n].empty() ? "" : data.atoms.names[n].substr(0, 1);
    };
    std::function<std::string(const size_t&)> truncatedResidueNames = [&](size_t n)
    {
        return codeUtils::truncate(3, data.residues.names[n]);
    };

    std::vector<bool> residueIncluded = codeUtils::indicesToBools(data.indices.residueCount, selectedResidues);
    assembly::Graph graph             = cds::createAssemblyGraph(data.indices, data.atomGraph);
    cds::PdbFileResidueData residueData {
        data.residues.numbers,
        codeUtils::vectorMap(truncatedResidueNames, codeUtils::indexVector(data.indices.residueCount)),
        std::vector<std::string>(data.indices.residueCount, ""), data.residues.insertionCodes};
    std::vector<std::string> elements =
        codeUtils::vectorMap(elementString, codeUtils::indexVector(data.indices.atomCount));
    cds::PdbFileFormat format;

    size_t modelCount = data.trajectory.coordinates.size();
    for (size_t coordinateSet = 0; coordinateSet < modelCount; coordinateSet++)
    {
        cds::PdbFileAtomData atomData {data.trajectory.coordinates[coordinateSet],
                                       data.atoms.numbers,
                                       data.atoms.names,
                                       elements,
                                       data.atoms.recordNames,
                                       data.atoms.occupancies,
                                       data.atoms.temperatureFactors};
        cds::PdbFileData writerData {format, {}, residueData, atomData};
        stream << "MODEL " << std::right << std::setw(8) << (coordinateSet + 1) << "\n";
        for (size_t moleculeId = 0; moleculeId < data.indices.moleculeCount; moleculeId++)
        {
            std::vector<size_t> residueIds = codeUtils::boolsToIndices(codeUtils::vectorAnd(
                residueIncluded, codeUtils::vectorEquals(data.indices.residueMolecule, moleculeId)));
            if (residueIds.size() > 0)
            {
                std::vector<cds::ResidueType> types = codeUtils::indicesToValues(data.residues.types, residueIds);
                std::vector<bool> ter               = residueTER(types);
                cds::writeMoleculeToPdb(stream, graph, residueIds, ter, writerData);
            }
        }
        stream << "ENDMDL\n";
    }
    theEnd(stream);
}

void cds::WritePdb(std::ostream& stream, const GraphIndexData& graphData, const std::vector<std::string>& headerLines)
{
    std::vector<bool> includedAtoms       = cds::atomVisibility(graphData.objects.atoms);
    graph::Database atomGraphData         = cds::createGraphData(graphData.objects);
    atomGraphData.nodeAlive               = cds::atomVisibility(graphData.objects.atoms);
    assembly::Graph graph                 = createAssemblyGraph(graphData.indices, atomGraphData);
    const std::vector<Residue*>& residues = graphData.objects.residues;
    std::vector<ResidueType> types        = residueTypes(residues);
    std::vector<bool> ter                 = residueTER(types);
    PdbFileData data                      = toPdbFileData(graphData.objects);
    for (auto& line : headerLines)
    {
        stream << "HEADER    " << line << "\n";
    }
    cds::writeMoleculeToPdb(stream, graph, codeUtils::indexVector(residues), ter, data);
    std::vector<bool> selectedResidueTypes(ResidueTypeCount, false);
    for (auto type : {Sugar, Derivative, Aglycone, Undefined})
    {
        selectedResidueTypes[type] = true;
    }
    std::vector<bool> residueSelected = codeUtils::indicesToValues(selectedResidueTypes, types);
    std::vector<bool> atomSelected    = codeUtils::indicesToValues(residueSelected, graphData.indices.atomResidue);
    graph::Database subgraphData      = cds::createGraphData(graphData.objects);
    subgraphData.nodeAlive            = codeUtils::vectorAnd(includedAtoms, atomSelected);
    assembly::Graph subgraph          = createAssemblyGraph(graphData.indices, subgraphData);
    std::vector<std::array<size_t, 2>> connectionIndices =
        codeUtils::indicesToValues(subgraph.residues.source.edgeNodes, subgraph.residues.edges.indices);
    cds::writeConectCards(stream, data.atoms.numbers, connectionIndices);
    cds::theEnd(stream);
}
