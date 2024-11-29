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
        const GlycanIndices& glycan = graphs.indices.glycans[glycanId];
        std::vector<cds::ResidueLinkageShapePreference> result;
        for (size_t n = 0; n < glycan.linkages.size(); n++)
        {
            size_t linkageId    = glycan.linkages[n];
            bool isFirstLinkage = n == 0;
            std::vector<std::vector<double>> angles;
            auto& linkageMetadata = data.residueLinkageData.metadata[linkageId];
            angles.resize(linkageMetadata.size());
            for (size_t k = 0; k < linkageMetadata.size(); k++)
            {
                for (auto& metadata : linkageMetadata[k])
                {
                    angles[k].push_back(randomAngle(rng, metadata));
                }
            }
            if (data.residueLinkageData.rotamerTypes[linkageId] == GlycamMetadata::RotamerType::conformer)
            {
                const std::vector<size_t>& rotatableDihedrals =
                    graphs.indices.residueLinkages[linkageId].rotatableDihedrals;
                auto order                         = randomMetadata(rng, linkageMetadata[0]);
                auto isFrozen                      = std::vector<bool>(rotatableDihedrals.size(), false);
                cds::ConformerShapePreference pref = {isFrozen, angles, order};
                if (isFirstLinkage && freezeGlycositeResidueConformation)
                {
                    for (size_t n = 0; n < rotatableDihedrals.size(); n++)
                    {
                        size_t dihedralId = rotatableDihedrals[n];
                        auto& name        = data.residueLinkageData.metadata[linkageId][n][0].dihedral_angle_name_;
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
                order.reserve(linkageMetadata.size());
                for (auto& metadataVector : linkageMetadata)
                {
                    order.push_back(randomMetadata(rng, metadataVector));
                }
                result.push_back(cds::PermutationShapePreference {angles, order});
            }
        }
        return result;
    };
} // namespace glycoproteinBuilder
