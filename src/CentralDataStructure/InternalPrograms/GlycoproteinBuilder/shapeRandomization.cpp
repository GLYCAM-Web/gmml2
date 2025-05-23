#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/shapeRandomization.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycoproteinStructs.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycanShape.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralShape.hpp"
#include "includes/CodeUtils/constants.hpp"
#include "includes/External_Libraries/PCG/pcg_random.h"

#include <functional>
#include <vector>

namespace
{
    bool isChiAngle(const std::string& str)
    {
        return (str == "Chi1") || (str == "Chi2");
    }
} // namespace

namespace glycoproteinBuilder
{
    cds::GlycanShapePreference randomLinkageShapePreference(
        pcg32& rng, const AngleSettings& settings, const AssemblyData& data, const assembly::Bounds& bounds,
        size_t glycanId,
        std::function<double(pcg32&, const AngleSettings&, const GlycamMetadata::DihedralAngleData& metadata)>
            randomAngle,
        bool freezeGlycositeResidueConformation)
    {
        const std::vector<size_t>& linkages                      = data.glycans.linkages[glycanId];
        const std::vector<std::vector<size_t>>& dihedralMetadata = data.rotatableDihedralData.metadata;
        cds::GlycanShapePreference result;
        for (size_t n = 0; n < linkages.size(); n++)
        {
            size_t linkageId                              = linkages[n];
            const std::vector<size_t>& rotatableDihedrals = data.indices.residueLinkages[linkageId].rotatableDihedrals;
            bool isGlycositeLinkage                       = data.residueLinkageData.isGlycositeLinkage[linkageId];
            std::vector<std::vector<double>> angles;
            angles.resize(rotatableDihedrals.size());
            for (size_t k = 0; k < rotatableDihedrals.size(); k++)
            {
                size_t dihedralId = rotatableDihedrals[k];
                for (auto& metadata : dihedralMetadata[dihedralId])
                {
                    angles[k].push_back(randomAngle(rng, settings, data.dihedralAngleTable.entries[metadata]));
                }
            }
            if (data.residueLinkageData.rotamerTypes[linkageId] == GlycamMetadata::RotamerType::conformer)
            {
                size_t firstDihedralId = rotatableDihedrals[0];
                std::vector<size_t> order =
                    settings.randomMetadata(rng, data.dihedralAngleTable, dihedralMetadata[firstDihedralId]);
                std::vector<bool> isFrozen(rotatableDihedrals.size(), false);
                cds::ConformerShapePreference pref {isFrozen, angles, order};
                if (isGlycositeLinkage && freezeGlycositeResidueConformation)
                {
                    for (size_t n = 0; n < rotatableDihedrals.size(); n++)
                    {
                        size_t dihedralId = rotatableDihedrals[n];
                        if (isChiAngle(
                                data.dihedralAngleTable.entries[dihedralMetadata[dihedralId][0]].dihedral_angle_name_))
                        {
                            pref.isFrozen[n] = true;
                            size_t metadata  = data.rotatableDihedralData.initialShape[dihedralId].metadataIndex;
                            pref.angles[n][metadata] =
                                constants::toDegrees(cds::angle(dihedralCoordinates(data, bounds, dihedralId)));
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
                    order.push_back(
                        settings.randomMetadata(rng, data.dihedralAngleTable, dihedralMetadata[dihedralId]));
                }
                result.push_back(cds::PermutationShapePreference {angles, order});
            }
        }
        return result;
    };
} // namespace glycoproteinBuilder
