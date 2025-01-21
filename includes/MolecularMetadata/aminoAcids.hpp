#ifndef INCLUDES_MOLECULARMETADATA_AMINOACIDS_HPP
#define INCLUDES_MOLECULARMETADATA_AMINOACIDS_HPP

#include <string>
#include <array>
#include <vector>

namespace MolecularMetadata
{
    typedef std::vector<std::pair<std::string, std::string>> BondVector;

    struct AminoAcid
    {
        std::vector<std::string> atomNames;
        BondVector bonds;
    };

    std::string originalResidueName(const std::string& str);
    const BondVector& carboxylBonds();
    const std::vector<std::string>& aminoAcidNames();
    const std::vector<AminoAcid>& aminoAcids();
    const AminoAcid& aminoAcid(size_t index);
    const AminoAcid& aminoAcid(const std::string& name);
    const std::vector<std::array<std::string, 4>>& aminoAcidDihedrals(size_t index);
} // namespace MolecularMetadata

#endif
