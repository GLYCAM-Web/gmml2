#ifndef INCLUDE_SEQUENCE_CARBOHYDRATE_HPP
#define INCLUDE_SEQUENCE_CARBOHYDRATE_HPP

#include "include/CentralDataStructure/Shapers/dihedralAngleSearchTypes.hpp"
#include "include/CentralDataStructure/Shapers/residueLinkageTypes.hpp"
#include "include/CentralDataStructure/molecule.hpp"
#include "include/metadata/dihedralangledata.hpp"
#include "include/readers/parameterManager.hpp"
#include "include/sequence/sequenceTypes.hpp"
#include "include/util/containerTypes.hpp"

#include <string>
#include <variant>
#include <vector>

namespace gmml
{
    void initializeCarbohydrate(
        Molecule& molecule,
        std::vector<ResidueLinkage>& glycosidicLinkages,
        const ParameterManager& parameterManager,
        const util::SparseVector<double>& elementRadii,
        const sequence::SequenceData& sequence);

    void resolveOverlaps(
        Molecule& molecule,
        std::vector<ResidueLinkage>& glycosidicLinkages,
        const util::SparseVector<double>& elementRadii,
        const DihedralAngleDataTable& metadataTable,
        const AngleSearchSettings& searchSettings);

    void generate3DStructureFiles(
        Molecule& molecule,
        const std::string& fileOutputDirectory,
        const std::string& outputFileNaming,
        const std::vector<std::string>& headerLines);

} // namespace gmml

#endif
