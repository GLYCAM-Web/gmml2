#include "includes/CentralDataStructure/CondensedSequence/sequenceManipulation.hpp"
#include "includes/CentralDataStructure/CondensedSequence/sequenceTypes.hpp"
#include "includes/Graph/graphManipulation.hpp"
#include "includes/CodeUtils/containers.hpp"

#include <string>
#include <vector>

using cds::ResidueType;

std::string cdsCondensedSequence::getLink(cds::ResidueType type, const std::string& linkage)
{
    switch (type)
    {
        case (ResidueType::Sugar):
            return linkage.substr(linkage.size() - 1, 1);
        case (ResidueType::Derivative):
            return linkage.substr(linkage.size() - 1, 1);
        case (ResidueType::Deoxy):
            return linkage.substr(0, 1);
        default:
            return "0";
    }
}
