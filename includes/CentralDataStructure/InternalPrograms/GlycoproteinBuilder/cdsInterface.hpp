#ifndef INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_CDSINTERFACE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_CDSINTERFACE_HPP

#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycoproteinStructs.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycosylationSite.hpp"
#include "includes/CentralDataStructure/molecule.hpp"

#include <vector>

namespace glycoproteinBuilder
{
    GlycoproteinAssembly toGlycoproteinAssemblyStructs(std::vector<cds::Molecule*>& molecules,
                                                       std::vector<GlycosylationSite>& glycosites,
                                                       const cds::OverlapProperties overlapProperties,
                                                       const OverlapMultiplier overlapWeight, bool excludeHydrogen);
} // namespace glycoproteinBuilder
#endif
