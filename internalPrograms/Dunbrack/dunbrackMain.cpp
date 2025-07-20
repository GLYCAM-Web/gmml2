#include "includes/CodeUtils/arguments.hpp"
#include "includes/CodeUtils/constants.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/files.hpp"
#include "includes/CodeUtils/filesystem.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/parsing.hpp"
#include "includes/CodeUtils/strings.hpp"
#include "includes/MolecularMetadata/aminoAcids.hpp"
#include "includes/MolecularMetadata/sidechainRotamers.hpp"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>

using MolecularMetadata::SidechainRotamerBin;
using MolecularMetadata::SidechainRotamerData;
using MolecularMetadata::SidechainRotation;

struct RotamerInput
{
    std::string residue;
    int phi;
    int psi;
    uint count;
    uint r1;
    uint r2;
    uint r3;
    uint r4;
    double probability;
    double chi1Val;
    double chi2Val;
    double chi3Val;
    double chi4Val;
    double chi1Sig;
    double chi2Sig;
    double chi3Sig;
    double chi4Sig;
};

std::vector<RotamerInput> readInputFile(const std::string& filename)
{
    std::vector<RotamerInput> input;

    auto processLine = [&](const std::string& line, size_t number)
    {
        if (!(line.empty() || line[0] == '#'))
        {
            std::vector<std::string> contents = codeUtils::split(line, ' ');
            if (contents.size() != 17)
            {
                throw std::runtime_error("incorrect line length: " + line);
            }
            const std::string& residue = contents[0];
            int phi = std::stoi(contents[1]);
            int psi = std::stoi(contents[2]);
            uint count = std::stoi(contents[3]);
            uint r1 = std::stoi(contents[4]);
            uint r2 = std::stoi(contents[5]);
            uint r3 = std::stoi(contents[6]);
            uint r4 = std::stoi(contents[7]);
            double probability = std::stod(contents[8]);
            double chi1Val = std::stod(contents[9]);
            double chi2Val = std::stod(contents[10]);
            double chi3Val = std::stod(contents[11]);
            double chi4Val = std::stod(contents[12]);
            double chi1Sig = std::stod(contents[13]);
            double chi2Sig = std::stod(contents[14]);
            double chi3Sig = std::stod(contents[15]);
            double chi4Sig = std::stod(contents[16]);

            input.push_back(
                {residue,
                 phi,
                 psi,
                 count,
                 r1,
                 r2,
                 r3,
                 r4,
                 probability,
                 chi1Val,
                 chi2Val,
                 chi3Val,
                 chi4Val,
                 chi1Sig,
                 chi2Sig,
                 chi3Sig,
                 chi4Sig});
        }
        return;
    };

    codeUtils::readFileLineByLine(filename, processLine);
    return input;
}

SidechainRotamerData formatOutput(const std::vector<RotamerInput>& input)
{
    std::vector<std::string> residues;
    std::vector<size_t> dihedralCount;
    std::vector<SidechainRotamerBin> bins;
    std::vector<SidechainRotation> lines;

    auto lineDihedrals = [](const RotamerInput& line)
    {
        const std::array<double, 4>& values = {line.chi1Val, line.chi2Val, line.chi3Val, line.chi4Val};
        size_t result = 0;
        for (auto& val : values)
        {
            result += (val != 0.0);
        }
        return result;
    };

    const RotamerInput& firstLine = input[0];
    residues.push_back(firstLine.residue);
    dihedralCount.push_back(0);
    bins.push_back({0, firstLine.phi, firstLine.psi});
    size_t currentBin = 0;
    for (auto& a : input)
    {
        const SidechainRotamerBin& bin = bins[currentBin];
        const std::string& residue = residues[bin.residue];
        if (!(a.psi == bin.psi && a.phi == bin.phi && a.residue == residue))
        {
            auto binIt = std::find_if(
                bins.begin(),
                bins.end(),
                [&](const SidechainRotamerBin& bin)
                { return a.psi == bin.psi && a.phi == bin.phi && a.residue == residues[bin.residue]; });
            if (binIt != bins.end())
            {
                throw std::runtime_error(
                    "input lines out of order, (residue, phi, psi) should be grouped together. Offender: (" +
                    a.residue + " " + std::to_string(a.phi) + " " + std::to_string(a.psi) + ")");
            }
            if (a.residue != residue)
            {
                size_t residueIndex = codeUtils::indexOf(residues, a.residue);
                if (residueIndex != residues.size())
                {
                    throw std::runtime_error("input line residues out of order. Offender: " + a.residue);
                }
                residues.push_back(a.residue);
                dihedralCount.push_back(0);
            }
            bins.push_back({residues.size() - 1, a.phi, a.psi});
            currentBin = bins.size() - 1;
        }
        size_t currentResidue = bins[currentBin].residue;
        dihedralCount[currentResidue] = std::max(dihedralCount[currentResidue], lineDihedrals(a));
        lines.push_back({
            currentBin, a.probability, {a.chi1Val, a.chi2Val, a.chi3Val, a.chi4Val}
        });
    }
    return SidechainRotamerData {10, residues, dihedralCount, {}, bins, {}, lines};
}

void writeOutput(const SidechainRotamerData& output, const std::string& filename)
{
    std::ofstream os;
    os.open(filename.c_str());

    os << "# This file has been generated from:\n";
    os << "# Backbone-dependent rotamer library, April, 2011\n";
    os << "# Copyright (c) 2007-2010\n";
    os << "# Maxim V. Shapovalov and Roland L. Dunbrack Jr.\n";
    os << "# Fox Chase Cancer Center\n";
    os << "# Philadelphia, PA, USA\n";

    os << output.residues.size() << "\n";
    for (size_t n = 0; n < output.residues.size(); n++)
    {
        os << output.residues[n] << " " << output.residueDihedralCount[n] << "\n";
    }
    os << output.bins.size() << "\n";
    for (auto& bin : output.bins)
    {
        os << bin.residue << " " << bin.phi << " " << bin.psi << "\n";
    }
    os << output.rotations.size() << "\n";
    for (size_t n = 0; n < output.rotations.size(); n++)
    {
        const SidechainRotation& line = output.rotations[n];
        size_t dihedralCount = output.residueDihedralCount[output.bins[line.bin].residue];
        os << line.bin << " " << line.probability;
        for (size_t n = 0; n < dihedralCount; n++)
        {
            os << " " << line.chi[n];
        }
        if (n < output.rotations.size() - 1)
        {
            os << "\n";
        }
    }
    os.close();
}

int main(int argc, char* argv[])
{
    enum ARGUMENTS
    {
        INPUT_FILE,
        OUTPUT_FILE,
        HELP
    };

    using codeUtils::ArgReq;
    using codeUtils::ArgType;
    std::vector<codeUtils::ArgDef> argumentDefinitions = {
        {ArgReq::required, ArgType::unnamed,  INPUT_FILE,     "", ' ',  "input-file"},
        {ArgReq::required, ArgType::unnamed, OUTPUT_FILE,     "", ' ', "output-file"},
        {ArgReq::optional,    ArgType::flag,        HELP, "help", 'h',            ""}
    };
    std::string programName = codeUtils::programName(argv);
    codeUtils::Arguments arguments;
    try
    {
        arguments = codeUtils::readArguments(argc, argv, argumentDefinitions);
        if (codeUtils::contains<int>(arguments.ids, HELP))
        {
            std::cout << codeUtils::helpString(programName, argumentDefinitions);
            std::exit(0);
        }
        codeUtils::validateArguments(arguments, argumentDefinitions);
    }
    catch (const std::runtime_error& error)
    {
        std::cout << "error in program arguments\n";
        std::cout << error.what() << "\n";
        std::cout << "\n";
        std::cout << codeUtils::helpString(programName, argumentDefinitions);
        std::exit(1);
    }

    std::string inputFile = "";
    std::string outputFile = "";
    for (const auto& arg : arguments.args)
    {
        switch (arg.id)
        {
            case ARGUMENTS::INPUT_FILE:
                {
                    inputFile = arg.value;
                    break;
                }
            case ARGUMENTS::OUTPUT_FILE:
                {
                    outputFile = arg.value;
                    break;
                }
            default:
                break;
        }
    }

    const std::vector<RotamerInput> input = readInputFile(inputFile);
    const std::vector<std::string>& aminoAcids = MolecularMetadata::aminoAcidNames();

    std::function<bool(const RotamerInput&)> includeLine = [&aminoAcids](const RotamerInput& line)
    {
        size_t index = codeUtils::indexOf(aminoAcids, line.residue);
        return (index < aminoAcids.size()) && (MolecularMetadata::aminoAcidDihedrals(index).size() > 0) &&
               (line.probability > 0.0);
    };

    const std::vector<RotamerInput> included = codeUtils::filter(includeLine, input);
    const MolecularMetadata::SidechainRotamerData output = formatOutput(included);
    writeOutput(output, outputFile);

    return 0;
}
