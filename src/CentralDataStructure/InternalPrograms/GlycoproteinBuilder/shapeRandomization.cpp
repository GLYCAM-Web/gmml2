#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/shapeRandomization.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycoproteinStructs.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycanShape.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralShape.hpp"
#include "includes/CodeUtils/constants.hpp"

#include <functional>
#include <vector>

namespace glycoproteinBuilder
{
    std::vector<cds::ResidueLinkageShapePreference> randomLinkageShapePreference(
        const AssemblyGraphs& graphs, const AssemblyData& data, size_t glycanId,
        std::function<std::vector<size_t>(GlycamMetadata::DihedralAngleDataVector metadataVector)> randomMetadata,
        std::function<double(GlycamMetadata::DihedralAngleData metadata)> randomAngle,
        bool freezeGlycositeResidueConformation)
    {
        const GlycanIndices& glycan = graphs.glycans[glycanId];
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
                    angles[k].push_back(randomAngle(metadata));
                }
            }
            if (data.residueLinkageData.rotamerTypes[linkageId] == GlycamMetadata::RotamerType::conformer)
            {
                const std::vector<size_t>& rotatableDihedrals = graphs.residueLinkages[linkageId].rotatableDihedrals;
                auto order                                    = randomMetadata(linkageMetadata[0]);
                auto isFrozen                                 = std::vector<bool>(rotatableDihedrals.size(), false);
                cds::ConformerShapePreference pref            = {isFrozen, angles, order};
                if (isFirstLinkage && freezeGlycositeResidueConformation)
                {
                    for (size_t n = 0; n < rotatableDihedrals.size(); n++)
                    {
                        size_t dihedralId = rotatableDihedrals[n];
                        auto& name        = data.residueLinkageData.metadata[linkageId][n][0].dihedral_angle_name_;
                        if ((name == "Chi1") || (name == "Chi2"))
                        {
                            pref.isFrozen[n] = true;
                            size_t metadata  = data.rotatableDihedralData.currentShape[dihedralId].metadataIndex;
                            pref.angles[n][metadata] =
                                constants::toDegrees(cds::angle(dihedralCoordinates(graphs, data, dihedralId)));
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
                    order.push_back(randomMetadata(metadataVector));
                }
                result.push_back(cds::PermutationShapePreference {angles, order});
            }
        }
        return result;
    };
} // namespace glycoproteinBuilder
