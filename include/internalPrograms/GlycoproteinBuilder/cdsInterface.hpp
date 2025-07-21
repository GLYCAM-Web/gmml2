#ifndef INCLUDES_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_CDSINTERFACE_HPP
#define INCLUDES_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_CDSINTERFACE_HPP

#include "include/CentralDataStructure/cdsTypes.hpp"
#include "include/internalPrograms/GlycoproteinBuilder/glycoproteinCreation.hpp"
#include "include/internalPrograms/GlycoproteinBuilder/glycoproteinStructs.hpp"
#include "include/metadata/aminoAcids.hpp"
#include "include/readers/Pdb/pdbData.hpp"
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
            std::vector<Carbohydrate*>& glycans,
            double overlapTolerance,
            double overlapRejectionThreshold,
            bool excludeHydrogen);
    }
} // namespace gmml

#endif
