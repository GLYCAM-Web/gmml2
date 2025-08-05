#ifndef INCLUDE_PROGRAMS_GLYCOPROTEINBUILDER_GLYCOPROTEINCREATION_HPP
#define INCLUDE_PROGRAMS_GLYCOPROTEINBUILDER_GLYCOPROTEINCREATION_HPP

#include "include/CentralDataStructure/cdsTypes.hpp"
#include "include/fileType/pdb/pdbData.hpp"
#include "include/metadata/dihedralangledata.hpp"
#include "include/programs/GlycoproteinBuilder/gpInputStructs.hpp"
#include "include/readers/parameterManager.hpp"
#include "include/sequence/carbohydrate.hpp"
#include "include/util/containerTypes.hpp"

#include <string>
#include <vector>

namespace gmml
{
    namespace gpbuilder
    {
        struct GlycosylationSite
        {
            size_t residueId;
            GlycositeInput input;
        };

        std::vector<GlycosylationSite> createGlycosites(
            const pdb::PdbData& pdbData,
            size_t glycoproteinAssemblyId,
            const std::vector<GlycositeInput>& glycositesInputVector);

        void addGlycansToProtein(
            std::vector<Molecule*>& glycans,
            std::vector<std::vector<ResidueLinkage>>& glycosidicLinkages,
            const ParameterManager& parameterManager,
            const util::SparseVector<double>& elementRadii,
            const DihedralAngleDataTable& metadataTable,
            const pdb::PdbData& pdbData,
            Assembly* glycoprotein,
            const std::vector<GlycosylationSite>& glycosites);
    } // namespace gpbuilder
} // namespace gmml

#endif
