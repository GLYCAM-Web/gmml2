#include "include/glycoprotein/shapeRandomization.hpp"

#include "include/carbohydrate/dihedralShape.hpp"
#include "include/external/pcg/pcg_random.h"
#include "include/glycoprotein/glycanShape.hpp"
#include "include/glycoprotein/glycoproteinStructs.hpp"
#include "include/util/constants.hpp"

#include <functional>
#include <vector>

namespace gmml
{
    namespace gpbuilder
    {
        namespace
        {
            bool isChiAngle(const std::string& str) { return (str == "Chi1") || (str == "Chi2"); }
        } // namespace

        GlycanShapePreference randomLinkageShapePreference(
            pcg32& rng,
            const DihedralAngleDataTable& dihedralAngleTable,
            const AngleSettings& settings,
            const AssemblyData& data,
            const assembly::Bounds& bounds,
            size_t glycanId,
            std::function<double(pcg32&, const AngleSettings&, const DihedralAngleData& metadata)> randomAngle,
            bool freezeGlycositeResidueConformation)
        {
            const std::vector<size_t>& linkages = data.glycans.linkages[glycanId];
            const std::vector<std::vector<size_t>>& dihedralMetadata = data.rotatableBonds.metadata;
            GlycanShapePreference result;
            for (size_t n = 0; n < linkages.size(); n++)
            {
                size_t linkageId = linkages[n];
                const std::vector<size_t>& rotatableBonds = data.residueLinkages.rotatableBonds[linkageId];
                bool isGlycositeLinkage = data.residueLinkages.isGlycositeLinkage[linkageId];
                std::vector<std::vector<double>> angles;
                angles.resize(rotatableBonds.size());
                for (size_t k = 0; k < rotatableBonds.size(); k++)
                {
                    size_t bondId = rotatableBonds[k];
                    for (auto& metadata : dihedralMetadata[bondId])
                    {
                        angles[k].push_back(randomAngle(rng, settings, dihedralAngleTable.entries[metadata]));
                    }
                }
                if (data.residueLinkages.rotamerTypes[linkageId] == RotamerType::conformer)
                {
                    size_t firstbondId = rotatableBonds[0];
                    std::vector<size_t> order =
                        settings.randomMetadata(rng, dihedralAngleTable, dihedralMetadata[firstbondId]);
                    std::vector<bool> isFrozen(rotatableBonds.size(), false);
                    ConformerShapePreference pref {isFrozen, angles, order};
                    if (isGlycositeLinkage && freezeGlycositeResidueConformation)
                    {
                        for (size_t n = 0; n < rotatableBonds.size(); n++)
                        {
                            size_t bondId = rotatableBonds[n];
                            if (isChiAngle(
                                    dihedralAngleTable.entries[dihedralMetadata[bondId][0]].dihedral_angle_name_))
                            {
                                pref.isFrozen[n] = true;
                                size_t metadata = data.rotatableBonds.initialShape[bondId].metadataIndex;
                                pref.angles[n][metadata] =
                                    constants::toDegrees(angle(dihedralCoordinates(data, bounds, bondId)));
                                pref.metadataOrder = {metadata};
                            }
                        }
                    }
                    result.push_back(pref);
                }
                else
                {
                    std::vector<std::vector<size_t>> order;
                    order.reserve(rotatableBonds.size());
                    for (size_t bondId : rotatableBonds)
                    {
                        order.push_back(settings.randomMetadata(rng, dihedralAngleTable, dihedralMetadata[bondId]));
                    }
                    result.push_back(PermutationShapePreference {angles, order});
                }
            }
            return result;
        };
    } // namespace gpbuilder
} // namespace gmml
