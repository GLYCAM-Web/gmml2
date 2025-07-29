#include "include/pdb/pdbResidueId.hpp"

#include "include/pdb/pdbFunctions.hpp"
#include "include/util/constants.hpp"
#include "include/util/containers.hpp"
#include "include/util/strings.hpp"

#include <functional>
#include <string>
#include <vector>

namespace gmml
{
    namespace pdb
    {
        ResidueId readResidueId(const std::string& line)
        { // atom number can overrun, shifting everything to the right.
            int shift = checkShiftFromSerialNumberOverrun(line);
            std::string alternativeLocation = util::RemoveWhiteSpace(line.substr(
                16 + shift, 1)); // I created this to use it only when reading, so I can discard the additional entries.
            std::string residueName = util::RemoveWhiteSpace(line.substr(17 + shift, 3));
            std::string chainId = util::RemoveWhiteSpace(line.substr(21 + shift, 1));
            int secondShift = util::GetSizeOfIntInString(line.substr(26 + shift));
            std::string sequenceNumber = util::RemoveWhiteSpace(line.substr(22 + shift, 4 + secondShift));
            // Insertion code gets shifted right by every overrun in residue number.
            std::string insertionCode = util::RemoveWhiteSpace(line.substr(26 + shift + secondShift, 1));
            return {residueName, sequenceNumber, insertionCode, chainId, alternativeLocation};
        }

        ResidueId readResidueId(const std::vector<std::string>& inputVector)
        {
            auto blankOrValue = [](const std::string& str) { return (str == "?") ? "" : str; };
            if (inputVector.size() < 4)
            {
                throw std::runtime_error("ResidueId cannot be constructed from inputs");
            }
            std::string residueName = blankOrValue(inputVector[0]);
            std::string sequenceNumber = blankOrValue(inputVector[1]);
            std::string insertionCode = blankOrValue(inputVector[2]);
            std::string chainId = blankOrValue(inputVector[3]);
            return {residueName, sequenceNumber, insertionCode, chainId, ""};
        }

        std::string toString(const ResidueId& id)
        {
            std::function<std::string(const std::string&)> strOrNotSet = [](const std::string& str)
            { return str.empty() ? constants::sNotSet : str; };
            return util::join(
                "_", util::vectorMap(strOrNotSet, {id.residueName, id.sequenceNumber, id.insertionCode, id.chainId}));
        }
    } // namespace pdb
} // namespace gmml
