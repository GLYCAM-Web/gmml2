#ifndef INCLUDES_MOLECULARMETADATA_PROTEINBONDING_HPP
#define INCLUDES_MOLECULARMETADATA_PROTEINBONDING_HPP

#include <string>
#include <vector>

namespace biology
{
    const std::vector<std::pair<std::string, std::string>>& getBackboneBonding();
    const std::vector<std::pair<std::string, std::string>>& getSidechainBonding(std::string queryResidueName);
    const std::string equivalentResidueNameLookup(const std::string queryResidueName);
} // namespace biology

#endif
