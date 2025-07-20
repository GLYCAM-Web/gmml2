#include "includes/CentralDataStructure/Readers/Pdb/pdbResidue.hpp"

#include "includes/CentralDataStructure/Readers/Pdb/pdbAtom.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbData.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbFunctions.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbResidueId.hpp"
#include "includes/CodeUtils/biology.hpp"
#include "includes/CodeUtils/casting.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/strings.hpp" //RemoveWhiteSpace

#include <string>

namespace
{
    using cds::ResidueType;

    ResidueType residueType(const std::string& name)
    {
        if (codeUtils::contains(biology::proteinResidueNames, name))
        {
            return ResidueType::Protein;
        }
        // ToDo we want to figure out solvent, aglycone etc here too?.
        return ResidueType::Undefined;
    }
} // namespace

pdb::ResidueId pdb::pdbResidueId(const PdbData& data, size_t residueId)
{
    return ResidueId(
        data.residues.names[residueId],
        std::to_string(data.residues.numbers[residueId]),
        data.residues.insertionCodes[residueId],
        data.residues.chainIds[residueId]);
}

std::string pdb::residueStringId(const PdbData& data, size_t residueId)
{
    return pdbResidueId(data, residueId).print();
}

size_t pdb::readResidue(PdbData& data, size_t moleculeId, std::stringstream& singleResidueSecion, std::string firstLine)
{
    ResidueId resId(firstLine);
    size_t position = data.molecules.residueOrder[moleculeId].size();
    std::string name = resId.getName();
    size_t residueId = addResidue(
        data,
        moleculeId,
        position,
        {name,
         residueType(name),
         uint(std::stoi(resId.getNumber())),
         resId.getInsertionCode(),
         resId.getChainId(),
         false});
    std::string firstFoundAlternativeLocationIndicator = resId.getAlternativeLocation(); // Normally empty
    std::string line;
    while (getline(singleResidueSecion, line))
    {
        std::string recordName = codeUtils::RemoveWhiteSpace(line.substr(0, 6));
        if ((recordName == "ATOM") || (recordName == "HETATM"))
        {
            ResidueId id(line); // Check alternativeLocation and ignore any that aren't the same as the first one.
            if (id.getAlternativeLocation().empty() ||
                id.getAlternativeLocation() == firstFoundAlternativeLocationIndicator)
            { // If no alternative location (normal case) or alternateLocation is the first one (normally "A").
                addAtom(data, residueId, readAtom(line));
            }
            else if (firstFoundAlternativeLocationIndicator.empty())
            { // Sometimes first atom has one location, but later atoms have alternatives. Just set and use the first
              // one (normally "A")
                firstFoundAlternativeLocationIndicator = id.getAlternativeLocation();
                addAtom(data, residueId, readAtom(line));
            }
            else
            {
                gmml::log(__LINE__, __FILE__, gmml::INF, "Skipping line with alternative location: " + line);
            }
        }
    }
    return residueId;
}

size_t pdb::readResidue(
    PdbData& data,
    size_t moleculeId,
    const std::string& name,
    cds::ResidueType type,
    bool hasTerCard,
    size_t referenceResidue)
{
    size_t position = codeUtils::indexOf(data.molecules.residueOrder[moleculeId], referenceResidue) + 1;
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
std::string pdb::getNumberAndInsertionCode(const PdbData& data, size_t residueId)
{
    return std::to_string(data.residues.numbers[residueId]) + data.residues.insertionCodes[residueId];
}
