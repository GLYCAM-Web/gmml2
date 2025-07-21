#include "include/readers/Pdb/pdbResidueId.hpp"

#include "include/readers/Pdb/pdbFunctions.hpp"
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
        //////////////////////////////////////////////////////////
        //                       CONSTRUCTOR                    //
        //////////////////////////////////////////////////////////
        ResidueId::ResidueId(const std::string& line)
        { // atom number can overrun, shifting everything to the right.
            int shift = checkShiftFromSerialNumberOverrun(line);
            alternativeLocation_ = util::RemoveWhiteSpace(line.substr(
                16 + shift, 1)); // I created this to use it only when reading, so I can discard the additional entries.
            residueName_ = util::RemoveWhiteSpace(line.substr(17 + shift, 3));
            chainId_ = util::RemoveWhiteSpace(line.substr(21 + shift, 1));
            int secondShift = util::GetSizeOfIntInString(line.substr(26 + shift));
            sequenceNumber_ = util::RemoveWhiteSpace(line.substr(22 + shift, 4 + secondShift));
            // Insertion code gets shifted right by every overrun in residue number.
            insertionCode_ = util::RemoveWhiteSpace(line.substr(26 + shift + secondShift, 1));
        }

        ResidueId::ResidueId(
            const std::string name,
            const std::string number,
            const std::string insertionCode,
            const std::string chainId)
            : residueName_(name), sequenceNumber_(number), insertionCode_(insertionCode), chainId_(chainId)
        {}

        ResidueId::ResidueId(std::vector<std::string> inputVector)
        {
            if (inputVector.size() < 4)
            {
                throw std::runtime_error("ResidueId cannot be constructed from inputs");
            }
            if (inputVector.at(0) != "?")
            {
                residueName_ = inputVector.at(0);
            }
            if (inputVector.at(1) != "?")
            {
                sequenceNumber_ = inputVector.at(1);
            }
            if (inputVector.at(2) != "?")
            {
                insertionCode_ = inputVector.at(2);
            }
            if (inputVector.at(3) != "?")
            {
                chainId_ = inputVector.at(3);
            }
            // Just leave alternativeLocation as empty.
        }

        //////////////////////////////////////////////////////////
        //                       FUNCTIONS                      //
        //////////////////////////////////////////////////////////
        const std::string ResidueId::getNumberAndInsertionCode() const
        {
            return this->getNumber() + this->getInsertionCode();
        }

        //////////////////////////////////////////////////////////
        //                       DISPLAY                        //
        //////////////////////////////////////////////////////////
        std::string ResidueId::print() const
        {
            std::function<std::string(const std::string&)> strOrNotSet = [](const std::string& str)
            { return str.empty() ? constants::sNotSet : str; };
            return util::join(
                "_",
                util::vectorMap(
                    strOrNotSet, {this->getName(), this->getNumber(), this->getInsertionCode(), this->getChainId()}));
        }
    } // namespace pdb
} // namespace gmml
