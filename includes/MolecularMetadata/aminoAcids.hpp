#ifndef INCLUDES_MOLECULARMETADATA_AMINOACIDS_HPP
#define INCLUDES_MOLECULARMETADATA_AMINOACIDS_HPP

#include <string>
#include <array>
#include <vector>

namespace MolecularMetadata
{
    typedef std::vector<std::pair<std::string, std::string>> BondVector;

    struct AminoAcidTable
    {
        std::vector<std::string> names;
        std::vector<std::string> originalName;
        std::vector<double> weights;
        std::vector<std::vector<std::string>> atomNames;
        std::vector<BondVector> bonds;
        std::vector<std::vector<std::array<std::string, 4>>> sidechainDihedralAtoms;
    };

    const AminoAcidTable& aminoAcidTable();
    size_t aminoAcidIndex(const AminoAcidTable& table, const std::string& name);
    std::string originalResidueName(const std::string& str);
    const BondVector& carboxylBonds();
} // namespace MolecularMetadata

#endif
