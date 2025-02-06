#ifndef INCLUDES_MOLECULARMETADATA_SIDECHAINROTAMERS_HPP
#define INCLUDES_MOLECULARMETADATA_SIDECHAINROTAMERS_HPP

#include <string>
#include <array>
#include <vector>

namespace MolecularMetadata
{
    struct SidechainRotamerBin
    {
        size_t residue;
        int phi;
        int psi;
    };

    struct SidechainRotation
    {
        size_t bin;
        double probability;
        std::array<double, 4> chi;
    };

    struct SidechainRotamerData
    {
        uint binSize;
        std::vector<std::string> residues;
        std::vector<size_t> residueDihedralCount;
        std::vector<std::pair<size_t, size_t>> residueBins;
        std::vector<SidechainRotamerBin> bins;
        std::vector<std::pair<size_t, size_t>> binRotations;
        std::vector<SidechainRotation> rotations;
    };

    std::vector<size_t> sidechainRotationIndices(const SidechainRotamerData& data, const std::string& residue,
                                                 double phi, double psi);
    SidechainRotamerData readSidechainRotamerData(const std::string& filename);
} // namespace MolecularMetadata

#endif
