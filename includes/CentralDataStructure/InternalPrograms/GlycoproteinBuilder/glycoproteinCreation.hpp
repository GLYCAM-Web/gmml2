#ifndef INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GLYCOPROTEINCREATION_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GLYCOPROTEINCREATION_HPP

#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycosylationSite.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/gpInputStructs.hpp"
#include "includes/CentralDataStructure/assembly.hpp"
#include "includes/CentralDataStructure/residue.hpp"

#include <string>
#include <vector>

namespace glycoproteinBuilder
{
    cds::Residue* selectResidueFromInput(cds::Assembly* glycoprotein, const std::string userSelection);
    std::vector<GlycosylationSite> createGlycosites(cds::Assembly* glycoprotein,
                                                    std::vector<GlycositeInput> glycositesInputVector);
} // namespace glycoproteinBuilder

#endif
