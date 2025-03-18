#ifndef INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GLYCOPROTEINCREATION_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GLYCOPROTEINCREATION_HPP

#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/gpInputStructs.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/assembly.hpp"
#include "includes/CentralDataStructure/CondensedSequence/carbohydrate.hpp"

#include <string>
#include <vector>

namespace glycoproteinBuilder
{
    struct GlycosylationSite
    {
        cds::Residue* residue;
        cdsCondensedSequence::Carbohydrate* glycan;
        GlycositeInput input;
    };

    std::vector<GlycosylationSite> createGlycosites(cds::Assembly* glycoprotein,
                                                    std::vector<GlycositeInput> glycositesInputVector);
} // namespace glycoproteinBuilder

#endif
