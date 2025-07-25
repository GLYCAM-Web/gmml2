#ifndef INCLUDE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GLYCOPROTEINCREATION_HPP
#define INCLUDE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GLYCOPROTEINCREATION_HPP

#include "include/CentralDataStructure/cdsTypes.hpp"
#include "include/internalPrograms/GlycoproteinBuilder/gpInputStructs.hpp"
#include "include/metadata/dihedralangledata.hpp"
#include "include/readers/Pdb/pdbData.hpp"
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
            Residue* residue;
            GlycositeInput input;
        };

        std::vector<GlycosylationSite> createGlycosites(
            const pdb::PdbData& pdbData,
            Assembly* glycoprotein,
            const std::vector<GlycositeInput>& glycositesInputVector);

        std::vector<Carbohydrate*> addGlycansToProtein(
            const ParameterManager& parameterManager,
            const util::SparseVector<double>& elementRadii,
            const DihedralAngleDataTable& metadataTable,
            const pdb::PdbData& pdbData,
            Assembly* glycoprotein,
            const std::vector<GlycosylationSite>& glycosites);
    } // namespace gpbuilder
} // namespace gmml

#endif
