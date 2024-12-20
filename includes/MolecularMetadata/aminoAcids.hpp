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

    struct AminoAcid
    {
        MoleculeDefinition standard;
        MoleculeDefinition zwitterion;
    };

    const std::vector<std::string>& aminoAcidNames();
    const std::vector<AminoAcid>& aminoAcids();
    const AminoAcid& aminoAcid(const std::string& name);
} // namespace MolecularMetadata

#endif
