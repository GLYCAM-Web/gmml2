#ifndef INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_CDSINTERFACE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_CDSINTERFACE_HPP

#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycoproteinStructs.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycoproteinCreation.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbData.hpp"
#include "includes/CentralDataStructure/molecule.hpp"
#include "includes/MolecularMetadata/aminoAcids.hpp"
#include "includes/CodeUtils/containerTypes.hpp"

#include <vector>

namespace glycoproteinBuilder
{
    GlycoproteinAssembly
    toGlycoproteinAssemblyStructs(const MolecularMetadata::AminoAcidTable& aminoAcidTable,
                                  const GlycamMetadata::DihedralAngleDataTable& dihedralAngleDataTable,
                                  const codeUtils::SparseVector<double>& elementRadii, const pdb::PdbData& pdbData,
                                  std::vector<cds::Molecule*>& molecules, std::vector<GlycosylationSite>& glycosites,
                                  std::vector<cdsCondensedSequence::Carbohydrate*>& glycans, double overlapTolerance,
                                  double overlapRejectionThreshold, bool excludeHydrogen);
} // namespace glycoproteinBuilder
#endif
