#ifndef INCLUDE_PDB_PDBRESIDUEID_HPP
#define INCLUDE_PDB_PDBRESIDUEID_HPP

#include <string>
#include <vector>

namespace gmml
{
    namespace pdb
    {
        struct ResidueId
        {
            std::string residueName;
            std::string sequenceNumber;
            std::string insertionCode;
            std::string chainId;
            std::string alternativeLocation;

            bool operator==(const ResidueId& rhs)
            {
                return (
                    (residueName == rhs.residueName) && (sequenceNumber == rhs.sequenceNumber) &&
                    (insertionCode == rhs.insertionCode) && (chainId == rhs.chainId));
            }
        };

        ResidueId readResidueId(const std::string& line);
        ResidueId readResidueId(const std::vector<std::string>& inputVector);
        std::string toString(const ResidueId& id);
    } // namespace pdb
} // namespace gmml

#endif
