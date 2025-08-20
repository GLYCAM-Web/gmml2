#include "include/fileType/pdb/pdbFileWriter.hpp"

#include "include/assembly/assemblyGraph.hpp"
#include "include/assembly/assemblyIndices.hpp"
#include "include/fileType/pdb/pdbData.hpp"
#include "include/fileType/pdb/pdbFileData.hpp"
#include "include/metadata/residueTypes.hpp"
#include "include/util/containers.hpp"
#include "include/util/formatting.hpp"
#include "include/util/strings.hpp"

#include <array>
#include <functional>
#include <iomanip>
#include <ostream>
#include <string>
#include <vector>

namespace gmml
{
    namespace pdb
    {
        std::vector<bool> residueTER(const std::vector<ResidueType>& types)
        {
            std::vector<bool> result;
            result.reserve(types.size());
            for (size_t n = 0; n < types.size(); n++)
            {
                const ResidueType& type = types[n];
                size_t next = n + 1;
                bool isSugar = type == ResidueType::Undefined || type == ResidueType::Sugar ||
                               type == ResidueType::Derivative || type == ResidueType::Aglycone;
                bool nextIsCapping = next < types.size() && types[next] == ResidueType::ProteinCappingGroup;
                bool betweenTwoCappingGroups = (type == ResidueType::ProteinCappingGroup) && nextIsCapping;
                bool isLast = next == types.size();
                bool isProtein = type == ResidueType::Protein;
                result.push_back(isSugar || betweenTwoCappingGroups || (isLast && isProtein));
            }
            return result;
        }

        void writeAssemblyToPdb(
            std::ostream& stream,
            const PdbFileData& data,
            const assembly::Graph& graph,
            const std::vector<std::vector<size_t>>& residueIndices,
            const std::vector<std::vector<bool>>& residueTER,
            const std::vector<std::array<size_t, 2>>& connectionIndices)
        {
            for (auto& line : data.headerLines)
            {
                stream << "HEADER    " << line << "\n";
            }
            for (size_t n = 0; n < residueIndices.size(); n++)
            {
                writeMoleculeToPdb(stream, data, graph, residueIndices[n], residueTER[n]);
            }
            writeConectCards(stream, data.atoms.numbers, connectionIndices);
        }

        void writeMoleculeToPdb(
            std::ostream& stream,
            const PdbFileData& data,
            const assembly::Graph& graph,
            const std::vector<size_t>& residueIndices,
            const std::vector<bool>& residueTER)
        {
            const PdbFileResidueData& residues = data.residues;
            const PdbFileAtomData& atoms = data.atoms;
            for (size_t n = 0; n < residueIndices.size(); n++)
            {
                size_t residueIndex = residueIndices[n];
                for (size_t residueAtom : residueAtoms(graph, residueIndex))
                {
                    size_t atomIndex = graph.residues.source.nodes.indices[residueAtom];
                    writeAtomToPdb(stream, data.format, residues, residueIndex, atoms, atomIndex);
                }
                if (residueTER[n])
                {
                    stream << "TER\n";
                }
            }
        }

        void writeTrajectoryToPdb(
            std::ostream& stream,
            const PdbData& data,
            const assembly::Graph& graph,
            const std::vector<size_t>& selectedResidues)
        {
            std::function<bool(const size_t&)> atomAlive = [&](size_t n) { return data.assembly.indices.atomAlive[n]; };
            std::function<std::string(const size_t&)> elementString = [&](size_t n)
            { return data.atoms.names[n].empty() ? "" : data.atoms.names[n].substr(0, 1); };
            std::function<std::string(const size_t&)> truncatedResidueNames = [&](size_t n)
            { return util::truncate(3, data.residues.names[n]); };

            std::vector<bool> residueIncluded =
                util::indicesToBools(data.assembly.indices.residueCount, selectedResidues);
            PdbFileResidueData residueData {
                data.residues.numbers,
                util::vectorMap(truncatedResidueNames, util::indexVector(data.assembly.indices.residueCount)),
                std::vector<std::string>(data.assembly.indices.residueCount, ""),
                data.residues.insertionCodes};
            PdbFileFormat format;

            size_t modelCount = data.trajectory.coordinates.size();
            for (size_t coordinateSet = 0; coordinateSet < modelCount; coordinateSet++)
            {
                PdbFileAtomData atomData {
                    data.trajectory.coordinates[coordinateSet],
                    data.atoms.numbers,
                    data.atoms.names,
                    elementNames(data.atoms.elements),
                    data.atoms.recordNames,
                    data.atoms.occupancies,
                    data.atoms.temperatureFactors};
                PdbFileData pdbData {format, {}, residueData, atomData};
                stream << "MODEL " << std::right << std::setw(8) << (coordinateSet + 1) << "\n";
                for (size_t moleculeId = 0; moleculeId < data.assembly.indices.moleculeCount; moleculeId++)
                {
                    std::vector<size_t> residueIds = util::boolsToIndices(
                        util::vectorAnd(residueIncluded, isMoleculeResidue(data.assembly.indices, moleculeId)));
                    if (residueIds.size() > 0)
                    {
                        std::vector<ResidueType> types = util::indicesToValues(data.residues.types, residueIds);
                        std::vector<bool> ter = residueTER(types);
                        writeMoleculeToPdb(stream, pdbData, graph, residueIds, ter);
                    }
                }
                stream << "ENDMDL\n";
            }
            theEnd(stream);
        }

        // Used by the PdbAtom class and cds classes. Thus what it takes as input is a bit odd.
        void writeAtomToPdb(
            std::ostream& stream,
            const PdbFileFormat& format,
            const PdbFileResidueData& residues,
            size_t residueIndex,
            const PdbFileAtomData& atoms,
            size_t atomIndex)
        {
            std::string residueAlternativeLocation = ""; // don't know if this should live in residue or atom..
            const Coordinate& coord = atoms.coordinates[atomIndex];
            stream << std::left << std::setw(6) << atoms.recordNames[atomIndex];
            stream << std::right << std::setw(5) << atoms.numbers[atomIndex] << std::left << std::setw(1) << " ";
            stream << std::left << std::setw(4) << atoms.names[atomIndex];
            stream << std::left << std::setw(1) << residueAlternativeLocation;
            stream << std::right << std::setw(3) << residues.names[residueIndex] << std::left << std::setw(1) << " ";
            stream << std::left << std::setw(1) << residues.chainIds[residueIndex];
            stream << std::right << std::setw(4) << residues.numbers[residueIndex];
            stream << std::left << std::setw(1) << residues.insertionCodes[residueIndex] << std::left << std::setw(3)
                   << " ";
            for (size_t n = 0; n < 3; n++)
            {
                util::writeFloat(stream, format.coordinate, coord.nth(n));
            }
            util::writeFloat(stream, format.occupancy, atoms.occupancies[atomIndex]);
            util::writeFloat(stream, format.temperatureFactor, atoms.temperatureFactors[atomIndex]);
            stream << std::left << std::setw(10) << " ";
            stream << std::right << std::setw(2) << atoms.elements[atomIndex];
            //    We probably don't want to write charges into pdb file. Width allowed is 2
            stream << std::endl;
        }

        void writeConectCards(
            std::ostream& stream,
            const std::vector<uint>& atomNumbers,
            const std::vector<std::array<size_t, 2>>& connectionIndices)
        { // These are only written for atoms connecting residues. The numbers overflow/truncate when longer than 5, but
          // the format is what the format is.
            std::function<bool(const std::array<size_t, 2>&, const std::array<size_t, 2>&)> ascending =
                [&](const std::array<size_t, 2> a, const std::array<size_t, 2> b)
            { return atomNumbers[a[0]] < atomNumbers[b[0]]; };
            std::vector<std::array<size_t, 2>> doubled = connectionIndices;
            doubled.reserve(connectionIndices.size() * 2);
            for (auto& a : connectionIndices)
            {
                doubled.push_back({a[1], a[0]});
            }
            for (auto& atomPair : util::sortedBy(ascending, doubled))
            {
                auto writeLine = [&](int a, int b)
                { stream << "CONECT" << std::right << std::setw(5) << a << std::right << std::setw(5) << b << "\n"; };
                writeLine(atomNumbers[atomPair[0]], atomNumbers[atomPair[1]]);
            }
        }

        void theEnd(std::ostream& stream) { stream << "END\n"; }
    } // namespace pdb
} // namespace gmml
