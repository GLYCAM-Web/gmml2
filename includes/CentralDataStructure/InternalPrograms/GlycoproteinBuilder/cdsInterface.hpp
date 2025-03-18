#ifndef INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_CDSINTERFACE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_CDSINTERFACE_HPP

#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycoproteinStructs.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycoproteinCreation.hpp"
#include "includes/CentralDataStructure/molecule.hpp"

#include <vector>

namespace glycoproteinBuilder
{
    GlycoproteinAssembly toGlycoproteinAssemblyStructs(std::vector<cds::Molecule*>& molecules,
                                                       std::vector<GlycosylationSite>& glycosites,
                                                       const OverlapMultiplier overlapWeight, double overlapTolerance,
                                                       double overlapRejectionThreshold, bool excludeHydrogen);
} // namespace glycoproteinBuilder
#endif
