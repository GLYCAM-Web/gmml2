#ifndef INCLUDES_MOLECULARMETADATA_AMINOACIDS_HPP
#define INCLUDES_MOLECULARMETADATA_AMINOACIDS_HPP

#include <string>
#include <vector>

namespace MolecularMetadata
{
    struct MoleculeDefinition
    {
        std::vector<std::string> names;
        std::vector<std::pair<std::string, std::string>> bonds;
    };

    const MoleculeDefinition& proteinBackbone();
    const std::vector<std::string>& aminoAcidNames();
    const std::vector<MoleculeDefinition>& aminoAcids();
    const MoleculeDefinition& aminoAcid(const std::string& name);
} // namespace MolecularMetadata

#endif
