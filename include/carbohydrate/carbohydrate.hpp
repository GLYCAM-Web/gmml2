#ifndef INCLUDE_CARBOHYDRATE_CARBOHYDRATE_HPP
#define INCLUDE_CARBOHYDRATE_CARBOHYDRATE_HPP

#include "include/CentralDataStructure/molecule.hpp"
#include "include/CentralDataStructure/residueLinkage/residueLinkageTypes.hpp"
#include "include/carbohydrate/carbohydrateTypes.hpp"
#include "include/carbohydrate/dihedralAngleSearchTypes.hpp"
#include "include/metadata/dihedralangledata.hpp"
#include "include/preprocess/parameterManager.hpp"
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
        const preprocess::ParameterManager& parameterManager,
        const util::SparseVector<double>& elementRadii,
        const sequence::SequenceData& sequence);

    void resolveOverlaps(
        Molecule& molecule,
        std::vector<ResidueLinkage>& glycosidicLinkages,
        const util::SparseVector<double>& elementRadii,
        const DihedralAngleDataTable& metadataTable,
        const AngleSearchSettings& searchSettings);

    void generate3DStructureFiles(
        const carbohydrate::CarbohydrateData& data,
        const std::string& name,
        const std::string& fileOutputDirectory,
        const std::string& outputFileNaming,
        const std::vector<std::string>& headerLines);

    carbohydrate::CarbohydrateData structured(const std::vector<Molecule*>& molecules);
    assembly::Graph createVisibleAssemblyGraph(const carbohydrate::CarbohydrateData& data);

} // namespace gmml

#endif
