#include "includes/MolecularMetadata/sidechainRotamers.hpp"

#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/files.hpp"

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

std::vector<size_t> MolecularMetadata::sidechainRotationIndices(
    const SidechainRotamerData& data, const std::string& residue, double phi, double psi)
{
    size_t residueId = codeUtils::indexOf(data.residues, residue);
    if (residueId == data.residues.size())
    {
        return {};
    }
    auto angleBinValue = [&](double angle) { return (int)(round(angle / data.binSize) * data.binSize); };
    const std::pair<size_t, size_t>& residueBin = data.residueBins[residueId];
    int phiBin = angleBinValue(phi);
    int psiBin = angleBinValue(psi);
    size_t firstBin = residueBin.first;
    size_t lastBin = firstBin + residueBin.second;
    auto binIt = std::find_if(
        data.bins.begin() + firstBin,
        data.bins.begin() + lastBin,
        [&](const SidechainRotamerBin& bin) { return bin.phi == phiBin && bin.psi == psiBin; });
    if (binIt == data.bins.begin() + lastBin)
    {
        throw std::runtime_error(
            "rotamer bin not found: " + residue + " " + std::to_string(phiBin) + " " + std::to_string(psiBin));
    }
    size_t binId = binIt - data.bins.begin();
    const std::pair<size_t, size_t>& binRotations = data.binRotations[binId];
    size_t firstRotation = binRotations.first;
    size_t rotationCount = binRotations.second;
    return codeUtils::indexVectorWithOffset(firstRotation, rotationCount);
}

MolecularMetadata::SidechainRotamerData MolecularMetadata::readSidechainRotamerData(const std::string& filename)
{
    std::vector<char> buffer = codeUtils::readEntireFile(filename);

    size_t start = 0;

    // skip comments
    while (buffer[start] == '#')
    {
        auto it = std::find(buffer.begin() + start, buffer.end(), '\n');
        start = it - buffer.begin() + 1;
    }

    char* p = &buffer[start];
    char* strEnd {};
    auto readUlong = [&p, &strEnd]()
    {
        ulong num = std::strtoul(p, &strEnd, 10);
        p = strEnd;
        return num;
    };
    auto readInt = [&p, &strEnd]()
    {
        int num = std::strtol(p, &strEnd, 10);
        p = strEnd;
        return num;
    };
    auto readDouble = [&p, &strEnd]()
    {
        double num = std::strtod(p, &strEnd);
        p = strEnd;
        return num;
    };
    auto readString = [&p](size_t length)
    {
        std::string str(p + 1, length);
        p += length + 1;
        return str;
    };
    size_t residueCount = readUlong();
    std::vector<std::string> residues;
    residues.reserve(residueCount);
    std::vector<size_t> dihedralCount;
    dihedralCount.reserve(residueCount);
    for (size_t n = 0; n < residueCount; n++)
    {
        residues.push_back(readString(3));
        size_t dihedrals = readUlong();
        dihedralCount.push_back(dihedrals);
    }
    size_t binCount = readUlong();
    std::vector<SidechainRotamerBin> bins;
    bins.reserve(binCount);
    std::vector<std::pair<size_t, size_t>> residueBins(residueCount, {0, 0});
    {
        size_t currentResidue = 0;
        size_t currentCount = 0;
        for (size_t n = 0; n < binCount; n++)
        {
            size_t residue = readUlong();
            if (residue != currentResidue)
            {
                residueBins[currentResidue].second = currentCount;
                residueBins[residue].first = n;
                currentResidue = residue;
                currentCount = 0;
            }
            currentCount += 1;
            int phi = readInt();
            int psi = readInt();
            bins.push_back({residue, phi, psi});
        }
        residueBins[currentResidue].second = currentCount;
    }
    size_t lineCount = readUlong();
    std::vector<SidechainRotation> lines;
    lines.reserve(lineCount);
    std::vector<std::pair<size_t, size_t>> binLines(binCount, {0, 0});
    {
        size_t currentBin = 0;
        size_t currentCount = 0;
        for (size_t n = 0; n < lineCount; n++)
        {
            size_t bin = readUlong();
            if (bin != currentBin)
            {
                binLines[currentBin].second = currentCount;
                binLines[bin].first = n;
                currentBin = bin;
                currentCount = 0;
            }
            currentCount += 1;
            double probability = readDouble();
            std::array<double, 4> dihedrals {0.0, 0.0, 0.0, 0.0};
            size_t dc = dihedralCount[bins[bin].residue];
            for (size_t k = 0; k < dc; k++)
            {
                dihedrals[k] = readDouble();
            }
            lines.push_back({bin, probability, dihedrals});
        }
        binLines[currentBin].second = currentCount;
    }

    uint binSize = 10;
    return SidechainRotamerData {binSize, residues, dihedralCount, residueBins, bins, binLines, lines};
}
