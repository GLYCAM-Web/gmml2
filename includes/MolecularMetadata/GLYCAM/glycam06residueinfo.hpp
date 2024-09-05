#ifndef INCLUDES_MOLECULARMETADATA_GLYCAM_GLYCAM06RESIDUEINFO_HPP
#define INCLUDES_MOLECULARMETADATA_GLYCAM_GLYCAM06RESIDUEINFO_HPP

#include <string>
#include <vector>

namespace gmml
{
    namespace MolecularMetadata
    {
        namespace GLYCAM
        {
            std::vector<std::string> getTypesForResidue(std::string query);
        } // namespace GLYCAM
    }     // namespace MolecularMetadata
} // namespace gmml

#endif
