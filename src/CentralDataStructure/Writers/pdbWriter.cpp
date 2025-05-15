#include "includes/CentralDataStructure/Writers/pdbWriter.hpp"
#include "includes/CentralDataStructure/FileFormats/pdbFileData.hpp"
#include "includes/CentralDataStructure/FileFormats/pdbFileWriter.hpp"
#include "includes/CentralDataStructure/cdsFunctions/cdsFunctions.hpp"
#include "includes/CentralDataStructure/cdsFunctions/graphInterface.hpp"
#include "includes/Assembly/assemblyGraph.hpp"
#include "includes/Graph/graphTypes.hpp"
#include "includes/Graph/graphManipulation.hpp"
#include "includes/CodeUtils/containers.hpp"

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

void cds::writeTrajectoryToPdb(std::ostream& stream, const std::vector<cds::Molecule*>& molecules)
{
    size_t modelCount = molecules.at(0)->getAtoms().at(0)->getNumberOfCoordinateSets();
    for (size_t coordinateSet = 0; coordinateSet < modelCount; coordinateSet++)
    {
        stream << "MODEL " << std::right << std::setw(8) << (coordinateSet + 1) << "\n";
        for (auto& molecule : molecules)
        {
            GraphIndexData graphData            = toIndexData({molecule});
            assembly::Graph graph               = createCompleteAssemblyGraph(graphData);
            std::vector<cds::Residue*> residues = molecule->getResidues();
            for (auto& atom : molecule->getAtoms())
            {
                atom->setCurrentCoordinate(coordinateSet);
            }
            std::vector<cds::ResidueType> types = residueTypes(residues);
            std::vector<bool> ter               = residueTER(types);
            PdbFileData pdbData                 = toPdbFileData(graphData.objects);
            cds::writeMoleculeToPdb(stream, graph, codeUtils::indexVector(residues), ter, pdbData);
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
