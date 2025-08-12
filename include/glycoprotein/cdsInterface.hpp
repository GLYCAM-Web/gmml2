#ifndef INCLUDE_GLYCOPROTEIN_CDSINTERFACE_HPP
#define INCLUDE_GLYCOPROTEIN_CDSINTERFACE_HPP

#include "include/CentralDataStructure/cdsTypes.hpp"
#include "include/fileType/pdb/pdbData.hpp"
#include "include/glycoprotein/glycoproteinCreation.hpp"
#include "include/glycoprotein/glycoproteinStructs.hpp"
#include "include/metadata/aminoAcids.hpp"
#include "include/util/containerTypes.hpp"

#include <vector>

namespace gmml
{
    namespace gpbuilder
    {
        GlycoproteinAssembly toGlycoproteinAssemblyStructs(
            const AminoAcidTable& aminoAcidTable,
            const DihedralAngleDataTable& dihedralAngleDataTable,
            const util::SparseVector<double>& elementRadii,
            const pdb::PdbData& pdbData,
            std::vector<Molecule*>& molecules,
            std::vector<GlycosylationSite>& glycosites,
            std::vector<Molecule*>& glycans,
            const std::vector<std::vector<ResidueLinkage>>& glycosidicLinkages,
            bool excludeHydrogen);
    }
} // namespace gmml

#endif
