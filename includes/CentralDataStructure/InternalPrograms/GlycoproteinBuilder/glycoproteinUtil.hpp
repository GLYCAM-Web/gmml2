#ifndef INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GLYCOPROTEINUTIL_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GLYCOPROTEINUTIL_HPP

#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycoproteinStructs.hpp"
#include "includes/Assembly/assemblyTypes.hpp"
#include "includes/CodeUtils/containers.hpp"

#include <vector>

namespace glycoproteinBuilder
{
    inline std::vector<bool> glycanIncluded(const AssemblyData& data, const std::vector<bool>& includedMolecules)
    {
        return codeUtils::indicesToValues(includedMolecules, data.glycans.moleculeId);
    }

    inline std::vector<size_t> includedGlycanMoleculeIds(const AssemblyData& data,
                                                         const std::vector<bool>& includedMolecules)
    {
        return codeUtils::indicesToValues(data.glycans.moleculeId,
                                          codeUtils::boolsToIndices(glycanIncluded(data, includedMolecules)));
    }

    inline std::vector<size_t> includedGlycanIndices(const AssemblyData& data,
                                                     const std::vector<bool>& includedMolecules)
    {
        return codeUtils::boolsToIndices(glycanIncluded(data, includedMolecules));
    }

    inline const std::vector<uint>& atomNumbers(bool serialized, const AssemblyData& data)
    {
        return serialized ? data.atoms.serializedNumbers : data.atoms.numbers;
    }

    inline const std::vector<uint>& residueNumbers(bool serialized, const AssemblyData& data)
    {
        return serialized ? data.residues.serializedNumbers : data.residues.numbers;
    }

    inline void deleteMolecule(MutableData& mutableData, size_t moleculeId)
    {
        mutableData.moleculeIncluded[moleculeId] = false;
    };
} // namespace glycoproteinBuilder
#endif
