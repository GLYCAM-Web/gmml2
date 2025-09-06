#ifndef INCLUDE_GLYCOPROTEIN_CDSINTERFACE_HPP
#define INCLUDE_GLYCOPROTEIN_CDSINTERFACE_HPP

#include "include/CentralDataStructure/graphInterface.hpp"
#include "include/fileType/pdb/pdbData.hpp"
#include "include/glycoprotein/glycoproteinStructs.hpp"
#include "include/metadata/aminoAcids.hpp"
#include "include/util/containerTypes.hpp"

#include <vector>

namespace gmml
{
    namespace gpbuilder
    {
        struct PartialLinkageData
        {
            RotatableBondData bonds;
            std::vector<size_t> initialDihedralMetadata;
            ResidueLinkageData linkages;
            std::vector<std::vector<size_t>> glycanLinkageIds;
        };

        PartialLinkageData linkageData(
            const DihedralAngleDataTable& dihedralAngleDataTable,
            const GraphIndexData& graphData,
            const assembly::Graph& graph,
            const std::vector<std::vector<ResidueLinkage>>& glycosidicLinkages);

        GlycoproteinAssembly toGlycoproteinAssemblyStructs(
            const AminoAcidTable& aminoAcidTable,
            const util::SparseVector<double>& elementRadii,
            const std::vector<bool>& includedElements,
            const pdb::PdbData& pdbData,
            const GraphIndexData& graphData,
            const assembly::Graph& graph,
            const std::vector<size_t>& glycanMoleculeIds,
            const std::vector<size_t>& glycositeIds,
            const PartialLinkageData& linkageData);
    } // namespace gpbuilder
} // namespace gmml

#endif
