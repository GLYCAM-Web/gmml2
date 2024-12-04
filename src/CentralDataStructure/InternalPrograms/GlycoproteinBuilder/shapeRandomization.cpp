#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/shapeRandomization.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycoproteinStructs.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycanShape.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralShape.hpp"
#include "includes/CodeUtils/constants.hpp"
#include "includes/External_Libraries/PCG/pcg_random.h"

#include <functional>
#include <vector>

namespace glycoproteinBuilder
{
    std::vector<cds::ResidueLinkageShapePreference> randomLinkageShapePreference(
        pcg32& rng, const AssemblyGraphs& graphs, const AssemblyData& data, const MutableData& mutableData,
        size_t glycanId,
        std::function<std::vector<size_t>(pcg32&, GlycamMetadata::DihedralAngleDataVector metadataVector)>
            randomMetadata,
        std::function<double(pcg32&, GlycamMetadata::DihedralAngleData metadata)> randomAngle,
        bool freezeGlycositeResidueConformation)
    {
        const GlycanIndices& glycan                                       = graphs.indices.glycans[glycanId];
        const std::vector<cds::DihedralAngleDataVector>& dihedralMetadata = data.rotatableDihedralData.metadata;
        std::vector<cds::ResidueLinkageShapePreference> result;
        for (size_t n = 0; n < glycan.linkages.size(); n++)
        {
            size_t linkageId = glycan.linkages[n];
            const std::vector<size_t>& rotatableDihedrals =
                graphs.indices.residueLinkages[linkageId].rotatableDihedrals;
            bool isGlycositeLinkage = data.residueLinkageData.isGlycositeLinkage[linkageId];
            std::vector<std::vector<double>> angles;
            angles.resize(rotatableDihedrals.size());
            for (size_t k = 0; k < rotatableDihedrals.size(); k++)
            {
                size_t dihedralId = rotatableDihedrals[k];
                for (auto& metadata : dihedralMetadata[dihedralId])
                {
                    angles[k].push_back(randomAngle(rng, metadata));
                }
            }
            if (data.residueLinkageData.rotamerTypes[linkageId] == GlycamMetadata::RotamerType::conformer)
            {
                size_t firstDihedralId             = rotatableDihedrals[0];
                auto order                         = randomMetadata(rng, dihedralMetadata[firstDihedralId]);
                auto isFrozen                      = std::vector<bool>(rotatableDihedrals.size(), false);
                cds::ConformerShapePreference pref = {isFrozen, angles, order};
                if (isGlycositeLinkage && freezeGlycositeResidueConformation)
                {
                    for (size_t n = 0; n < rotatableDihedrals.size(); n++)
                    {
                        size_t dihedralId = rotatableDihedrals[n];
                        auto& name        = data.rotatableDihedralData.metadata[dihedralId][0].dihedral_angle_name_;
                        if ((name == "Chi1") || (name == "Chi2"))
                        {
                            pref.isFrozen[n] = true;
                            size_t metadata  = mutableData.currentDihedralShape[dihedralId].metadataIndex;
                            pref.angles[n][metadata] =
                                constants::toDegrees(cds::angle(dihedralCoordinates(graphs, mutableData, dihedralId)));
                            pref.metadataOrder = {metadata};
                        }
                    }
                }
                result.push_back(pref);
            }
            else
            {
                std::vector<std::vector<size_t>> order;
                order.reserve(rotatableDihedrals.size());
                for (size_t dihedralId : rotatableDihedrals)
                {
                    order.push_back(randomMetadata(rng, dihedralMetadata[dihedralId]));
                }
                result.push_back(cds::PermutationShapePreference {angles, order});
            }
        }
        return result;
    };
} // namespace glycoproteinBuilder
