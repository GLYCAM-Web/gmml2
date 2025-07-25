#ifndef INCLUDE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GLYCOPROTEINUTIL_HPP
#define INCLUDE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GLYCOPROTEINUTIL_HPP

#include "include/assembly/assemblyTypes.hpp"
#include "include/internalPrograms/GlycoproteinBuilder/glycoproteinStructs.hpp"
#include "include/util/containers.hpp"

#include <vector>

namespace gmml
{
    namespace gpbuilder
    {
        inline std::vector<bool> glycanIncluded(const AssemblyData& data, const std::vector<bool>& includedMolecules)
        {
            return util::indicesToValues(includedMolecules, data.glycans.moleculeId);
        }

        inline std::vector<size_t> includedGlycanMoleculeIds(
            const AssemblyData& data, const std::vector<bool>& includedMolecules)
        {
            return util::indicesToValues(
                data.glycans.moleculeId, util::boolsToIndices(glycanIncluded(data, includedMolecules)));
        }

        inline std::vector<size_t> includedGlycanIndices(
            const AssemblyData& data, const std::vector<bool>& includedMolecules)
        {
            return util::boolsToIndices(glycanIncluded(data, includedMolecules));
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
    } // namespace gpbuilder
} // namespace gmml

#endif
