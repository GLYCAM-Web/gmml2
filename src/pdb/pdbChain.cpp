#include "include/pdb/pdbChain.hpp"

#include "include/CentralDataStructure/residue.hpp"
#include "include/geometry/geometryTypes.hpp"
#include "include/metadata/aminoAcids.hpp"
#include "include/pdb/pdbFunctions.hpp"
#include "include/pdb/pdbResidue.hpp"
#include "include/util/casting.hpp"
#include "include/util/containers.hpp"
#include "include/util/logging.hpp"
#include "include/util/strings.hpp"

#include <string>
#include <variant>
#include <vector>

namespace gmml
{
    namespace pdb
    {
        void readChain(PdbData& data, size_t moleculeId, std::stringstream& stream_block)
        {
            std::string line;
            while (getline(stream_block, line))
            {
                std::string recordName = util::RemoveWhiteSpace(line.substr(0, 6));
                if ((recordName == "ATOM") || (recordName == "HETATM"))
                {
                    std::stringstream singleResidueSection = extractSingleResidueFromRecordSection(stream_block, line);
                    readResidue(data, moleculeId, singleResidueSection, line);
                }
                else
                {
                    util::log(
                        __LINE__,
                        __FILE__,
                        util::WAR,
                        "In PdbChain Constructor with record that isn't cool: " + recordName);
                    break;
                }
            }
            tagTerminalResidues(data, moleculeId);
        }

        void tagTerminalResidues(PdbData& data, size_t moleculeId)
        {
            size_t nTer = getNTerminal(data, moleculeId);
            if (nTer < data.indices.residueCount)
            {
                data.residues.isNTerminal[nTer] = true;
                data.objects.residues[nTer]->setNTerminal(true);
            }
            size_t cTer = getCTerminal(data, moleculeId);
            if (cTer < data.indices.residueCount)
            {
                data.residues.isCTerminal[cTer] = true;
                data.objects.residues[cTer]->setCTerminal(true);
            }
        }

        std::stringstream extractSingleResidueFromRecordSection(std::stringstream& pdbFileStream, std::string line)
        {
            std::streampos previousLinePosition = pdbFileStream.tellg(); // Save current line position
            std::stringstream singleResidueSection;
            ResidueId residueId(line);
            ResidueId initialResidueId = residueId;
            while (residueId == initialResidueId)
            {
                singleResidueSection << line << std::endl;
                previousLinePosition = pdbFileStream.tellg(); // Save current line position.
                if (!std::getline(pdbFileStream, line))
                {
                    break; // // If we hit the end, time to leave.
                }
                residueId = ResidueId(line);
            }
            // Go back to previous line position. E.g. was reading HEADER and found TITLE:
            pdbFileStream.seekg(previousLinePosition);
            return singleResidueSection;
        }

        // Only makes sense for proteins.
        // Assumes vector is populated from N-terminal to C-terminal.
        size_t getNTerminal(const PdbData& data, size_t moleculeId)
        {
            std::function<bool(const size_t&)> isProtein = [&](size_t id)
            { return util::contains(proteinResidueNames, data.residues.names[id]); };
            std::vector<size_t> proteinResidues =
                util::vectorFilter(isProtein, data.molecules.residueOrder[moleculeId]);
            if (proteinResidues.empty())
            {
                util::log(__LINE__, __FILE__, util::WAR, "Looked for terminal residue of chain with protein residues.");
                return data.indices.residueCount;
            }
            return proteinResidues.front();
        }

        size_t getCTerminal(const PdbData& data, size_t moleculeId)
        {
            std::function<bool(const size_t&)> isProtein = [&](size_t id)
            { return util::contains(proteinResidueNames, data.residues.names[id]); };
            std::vector<size_t> proteinResidues =
                util::vectorFilter(isProtein, data.molecules.residueOrder[moleculeId]);
            if (proteinResidues.empty())
            {
                util::log(__LINE__, __FILE__, util::WAR, "Looked for terminal residue of chain with protein residues.");
                return data.indices.residueCount;
            }
            return proteinResidues.back();
        }
    } // namespace pdb
} // namespace gmml
