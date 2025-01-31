#include "includes/CentralDataStructure/FileFormats/offFileWriter.hpp"
#include "includes/CentralDataStructure/FileFormats/offFileData.hpp"
#include "includes/CentralDataStructure/residueTypes.hpp"
#include "includes/Assembly/assemblyGraph.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/formatting.hpp"

#include <vector>
#include <string>
#include <iomanip>

namespace
{
    std::string residueOffType(const cds::ResidueType queryType)
    {
        if (queryType == cds::ResidueType::Protein)
        {
            return "p";
        }
        else if (queryType == cds::ResidueType::Solvent)
        {
            return "w";
        }
        return "?";
    }
} // namespace

void cds::WriteOffFileUnit(std::ostream& stream, const OffFileFormat& format, const assembly::Graph& graph,
                           const OffFileResidueData& residues, const OffFileAtomData& atoms,
                           const std::vector<size_t>& residueIndices, const std::string& unitName)
{
    // WriteAtomSection
    const std::string FLAG = "131072";
    stream << "!entry." << unitName
           << ".unit.atoms table  str name  str type  int typex  int resx  int flags  int seq  int elmnt  dbl chg"
           << "\n";
    for (size_t residueIndex : residueIndices)
    {
        unsigned int atomNumberInResidue = 1;
        for (size_t residueAtom : assembly::residueAtoms(graph, residueIndex))
        {
            size_t atomIndex = graph.residues.source.nodes[residueAtom];
            stream << " \"" << atoms.names[atomIndex] << "\" "
                   << "\"" << atoms.types[atomIndex] << "\" "
                   << "0 " << residues.numbers[residueIndex] << " " << FLAG << " " << atomNumberInResidue << " "
                   << atoms.atomicNumbers[atomIndex] << " ";
            codeUtils::writeFloat(stream, format.charge, atoms.charges[atomIndex]);
            stream << "\n";
            atomNumberInResidue++;
        }
    }
    // WriteAtomPertInfoSection
    stream << "!entry." << unitName
           << ".unit.atomspertinfo table  str pname  str ptype  int ptypex  int pelmnt  dbl pchg"
           << "\n";
    for (size_t residueIndex : residueIndices)
    {
        for (size_t residueAtom : assembly::residueAtoms(graph, residueIndex))
        {
            size_t atomIndex = graph.residues.source.nodes[residueAtom];
            stream << " \"" << atoms.names[atomIndex] << "\""
                   << " \"" << atoms.types[atomIndex] << "\" 0 -1 0.0"
                   << "\n";
        }
    }
    // WriteBoundBoxSection
    stream << "!entry." << unitName << ".unit.boundbox array dbl"
           << "\n";
    stream << " -1.000000"
           << "\n";
    stream << " 0.0"
           << "\n";
    stream << " 0.0"
           << "\n";
    stream << " 0.0"
           << "\n";
    stream << " 0.0"
           << "\n";
    // WriteChildSequenceSection
    stream << "!entry." << unitName << ".unit.childsequence single int"
           << "\n";
    stream << " " << residueIndices.size() + 1 << "\n";
    // WriteConnectSection
    //  Note: this is silly but fine for most cases. If you're reading this it's because it mattered and you need to
    //  make it better.
    stream << "!entry." << unitName << ".unit.connect array int"
           << "\n";
    stream << " 1"
           << "\n";
    size_t lastAtom = assembly::residueAtoms(graph, residueIndices.back()).back();
    stream << " " << atoms.numbers[graph.residues.source.nodes[lastAtom]] << "\n";
    // WriteConnectivitySection
    stream << "!entry." << unitName << ".unit.connectivity table  int atom1x  int atom2x  int flags"
           << "\n";
    std::vector<bool> residueIncluded = codeUtils::indicesToBools(residues.names.size(), residueIndices);
    for (auto& bond : graph.atoms.edges.nodeAdjacencies)
    {
        size_t first  = graph.atoms.nodes.indices[bond[0]];
        size_t second = graph.atoms.nodes.indices[bond[1]];
        if (residueIncluded[graph.atomResidue[first]])
        {
            int number         = atoms.numbers[first];
            int neighborNumber = atoms.numbers[second];
            int min            = std::min(number, neighborNumber);
            int max            = std::max(number, neighborNumber);
            // According to docs: (the *second* atom is the one with the larger index). So ordering
            stream << " " << min << " " << max << " 1"
                   << "\n";
        }
    }
    // WriteHierarchySection
    stream << "!entry." << unitName << ".unit.hierarchy table  str abovetype  int abovex  str belowtype  int belowx"
           << "\n";
    for (size_t residueIndex : residueIndices)
    {
        stream << " \"U\" 0 \"R\" " << residues.numbers[residueIndex] << "\n";
        for (size_t residueAtom : assembly::residueAtoms(graph, residueIndex))
        {
            size_t atomIndex = graph.residues.source.nodes[residueAtom];
            stream << " \"R\" " << residues.numbers[residueIndex] << " \"A\" " << atoms.numbers[atomIndex] << "\n";
        }
    }
    // WriteNameSection
    stream << "!entry." << unitName << ".unit.name single str"
           << "\n";
    stream << " \"" << unitName << "\""
           << "\n";
    // WritePositionSection
    stream << "!entry." << unitName << ".unit.positions table  dbl x  dbl y  dbl z"
           << "\n";
    for (size_t residueIndex : residueIndices)
    {
        for (size_t residueAtom : assembly::residueAtoms(graph, residueIndex))
        {
            size_t atomIndex = graph.residues.source.nodes[residueAtom];
            Coordinate coord = atoms.coordinates[atomIndex];
            for (size_t n = 0; n < 3; n++)
            {
                stream << " ";
                codeUtils::writeFloat(stream, format.coordinate, coord.nth(n));
            }
            stream << "\n";
        }
    }
    // WriteResidueConnectSection // Every residue needs a head/tail regardless of reality. tleap uses this info.
    stream << "!entry." << unitName
           << ".unit.residueconnect table  int c1x  int c2x  int c3x  int c4x  int c5x  int c6x"
           << "\n";
    for (size_t residueIndex : residueIndices)
    {
        std::vector<size_t> connectedAtoms = residues.atomsConnectedToOtherResidues[residueIndex];
        // Deal with residues that don't have a tail/head in reality:
        if (connectedAtoms.size() == 1)
        { // Repeating the same atom changes the tree structure in the parm7 file. Not sure anything uses that. Old gmml
          // code repeats so doing that.
            // For reducing terminal old code puts a 1 2 0 0 0 0. So not repeat. Changing to 2 2 0 0 0 0 causes first
            // atom not to be M (main). Might be an issue let's see.
            connectedAtoms.push_back(connectedAtoms.front());
        }
        for (auto& atomIndex : connectedAtoms)
        {
            stream << " " << atoms.numbers[atomIndex];
        }
        int columnsWithZero = 6 - connectedAtoms.size();
        for (int i = 0; i < columnsWithZero; ++i)
        {
            stream << " 0";
        }
        stream << "\n";
    }
    // WriteResiduesSection
    stream << "!entry." << unitName
           << ".unit.residues table  str name  int seq  int childseq  int startatomx  str restype  int imagingx"
           << "\n";
    for (size_t residueIndex : residueIndices)
    {
        const std::vector<size_t>& atomIndices = assembly::residueAtoms(graph, residueIndex);
        unsigned int childseq                  = atomIndices.size() + 1;
        unsigned int startatomx                = atoms.numbers[graph.residues.source.nodes[atomIndices.front()]];
        std::string restype                    = residueOffType(residues.types[residueIndex]);
        unsigned int imagingx                  = 0;
        stream << " \"" << residues.names[residueIndex] << "\""
               << " " << residues.numbers[residueIndex] << " " << childseq << " " << startatomx << " "
               << "\"" << restype << "\""
               << " " << imagingx << "\n";
    }
    // WriteSolventCapSection
    stream << "!entry." << unitName << ".unit.solventcap array dbl"
           << "\n";
    stream << " -1.000000"
           << "\n";
    stream << " 0.0"
           << "\n";
    stream << " 0.0"
           << "\n";
    stream << " 0.0"
           << "\n";
    stream << " 0.0"
           << "\n";
    // WriteVelocitiesSection
    stream << "!entry." << unitName << ".unit.velocities table  dbl x  dbl y  dbl z"
           << "\n";
    for (size_t residueIndex : residueIndices)
    {
        for (size_t n = 0; n < assembly::residueAtoms(graph, residueIndex).size(); n++)
        { // Maybe later we'll want to deal with atom velocities...
            stream << " 0.0 0.0 0.0"
                   << "\n";
        }
    }
    return;
}

void cds::WriteResiduesIndividuallyToOffFile(std::ostream& stream, const assembly::Graph& graph,
                                             const OffFileData& data)
{ // For writing each residue separately
    size_t residueCount = data.residues.names.size();
    stream << "!!index array str"
           << "\n";
    for (size_t n = 0; n < residueCount; n++)
    {
        stream << " \"" << data.residues.names[n] << "\""
               << "\n";
    }
    for (size_t n = 0; n < residueCount; n++)
    {
        cds::WriteOffFileUnit(stream, data.format, graph, data.residues, data.atoms, {n}, data.residues.names[n]);
    }
}

void cds::WriteResiduesTogetherToOffFile(std::ostream& stream, const assembly::Graph& graph, const OffFileData& data,
                                         const std::string& unitName)
{ // For writing residues together as a molecule
    stream << "!!index array str"
           << "\n";
    stream << " \"" << unitName << "\""
           << "\n";
    cds::WriteOffFileUnit(stream, data.format, graph, data.residues, data.atoms,
                          codeUtils::indexVector(data.residues.names), unitName);
}
