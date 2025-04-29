#ifndef INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GLYCOPROTEINCREATION_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GLYCOPROTEINCREATION_HPP

#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/gpInputStructs.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/assembly.hpp"
#include "includes/CentralDataStructure/CondensedSequence/carbohydrate.hpp"
#include "includes/CentralDataStructure/Parameters/parameterManager.hpp"
#include "includes/MolecularMetadata/GLYCAM/dihedralangledata.hpp"
#include "includes/CodeUtils/containerTypes.hpp"

#include <string>
#include <vector>

namespace glycoproteinBuilder
{
    struct GlycosylationSite
    {
        cds::Residue* residue;
        GlycositeInput input;
    };

    std::vector<GlycosylationSite> createGlycosites(cds::Assembly* glycoprotein,
                                                    const std::vector<GlycositeInput>& glycositesInputVector);

    std::vector<cdsCondensedSequence::Carbohydrate*>
    addGlycansToProtein(const cdsParameters::ParameterManager& parameterManager, cds::Assembly* glycoprotein,
                        const codeUtils::SparseVector<double>& elementRadii,
                        const GlycamMetadata::DihedralAngleDataTable& metadataTable,
                        const std::vector<GlycosylationSite>& glycosites);
} // namespace glycoproteinBuilder

#endif
