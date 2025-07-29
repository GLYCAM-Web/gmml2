#include "include/pdb/pdbResidue.hpp"

#include "include/metadata/aminoAcids.hpp"
#include "include/pdb/pdbAtom.hpp"
#include "include/pdb/pdbData.hpp"
#include "include/pdb/pdbFunctions.hpp"
#include "include/pdb/pdbResidueId.hpp"
#include "include/util/containers.hpp"
#include "include/util/logging.hpp"
#include "include/util/strings.hpp" //RemoveWhiteSpace

#include <string>

namespace gmml
{
    namespace pdb
    {
        namespace
        {
            ResidueType residueType(const std::string& name)
            {
                if (util::contains(proteinResidueNames, name))
                {
                    return ResidueType::Protein;
                }
                // ToDo we want to figure out solvent, aglycone etc here too?.
                return ResidueType::Undefined;
            }
        } // namespace

        ResidueId pdbResidueId(const PdbData& data, size_t residueId)
        {
            return {
                data.residues.names[residueId],
                std::to_string(data.residues.numbers[residueId]),
                data.residues.insertionCodes[residueId],
                data.residues.chainIds[residueId],
                ""};
        }

        std::string residueStringId(const PdbData& data, size_t residueId)
        {
            return toString(pdbResidueId(data, residueId));
        }

        size_t readResidue(
            PdbData& data, size_t moleculeId, std::stringstream& singleResidueSecion, std::string firstLine)
        {
            ResidueId resId = readResidueId(firstLine);
            size_t position = data.molecules.residueOrder[moleculeId].size();
            std::string name = resId.residueName;
            size_t residueId = addResidue(
                data,
                moleculeId,
                position,
                {name,
                 residueType(name),
                 uint(std::stoi(resId.sequenceNumber)),
                 resId.insertionCode,
                 resId.chainId,
                 false});
            std::string firstFoundAlternativeLocationIndicator = resId.alternativeLocation; // Normally empty
            std::string line;
            while (getline(singleResidueSecion, line))
            {
                std::string recordName = util::RemoveWhiteSpace(line.substr(0, 6));
                if ((recordName == "ATOM") || (recordName == "HETATM"))
                {
                    ResidueId id = readResidueId(
                        line); // Check alternativeLocation and ignore any that aren't the same as the first one.
                    if (id.alternativeLocation.empty() ||
                        id.alternativeLocation == firstFoundAlternativeLocationIndicator)
                    { // If no alternative location (normal case) or alternateLocation is the first one (normally "A").
                        addAtom(data, residueId, readAtom(line));
                    }
                    else if (firstFoundAlternativeLocationIndicator.empty())
                    { // Sometimes first atom has one location, but later atoms have alternatives. Just set and use the
                      // first one (normally "A")
                        firstFoundAlternativeLocationIndicator = id.alternativeLocation;
                        addAtom(data, residueId, readAtom(line));
                    }
                    else
                    {
                        util::log(__LINE__, __FILE__, util::INF, "Skipping line with alternative location: " + line);
                    }
                }
            }
            return residueId;
        }

        size_t readResidue(
            PdbData& data,
            size_t moleculeId,
            const std::string& name,
            ResidueType type,
            bool hasTerCard,
            size_t referenceResidue)
        {
            size_t position = util::indexOf(data.molecules.residueOrder[moleculeId], referenceResidue) + 1;
            return addResidue(
                data,
                moleculeId,
                position,
                {name,
                 type,
                 data.residues.numbers[referenceResidue] + 1,
                 data.residues.insertionCodes[referenceResidue],
                 data.residues.chainIds[referenceResidue],
                 hasTerCard});
        }

        //////////////////////////////////////////////////////////
        //                         ACCESSOR                     //
        //////////////////////////////////////////////////////////
        std::string getNumberAndInsertionCode(const PdbData& data, size_t residueId)
        {
            return std::to_string(data.residues.numbers[residueId]) + data.residues.insertionCodes[residueId];
        }
    } // namespace pdb
} // namespace gmml
