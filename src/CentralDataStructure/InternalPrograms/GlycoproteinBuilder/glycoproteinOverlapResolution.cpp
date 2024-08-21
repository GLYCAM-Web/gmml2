#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycoproteinOverlapResolution.hpp"

#include "includes/CentralDataStructure/Geometry/overlap.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralShape.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralAngleSearch.hpp"
#include "includes/CentralDataStructure/Overlaps/atomOverlaps.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CodeUtils/metropolisCriterion.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/random.hpp"
#include "includes/External_Libraries/PCG/pcg_random.h"
#include "includes/External_Libraries/PCG/pcg_extras.h"

namespace
{
    using gmml::MolecularMetadata::GLYCAM::DihedralAngleData;
    using gmml::MolecularMetadata::GLYCAM::DihedralAngleDataVector;
    using gmml::MolecularMetadata::GLYCAM::RotamerType;

    std::vector<cds::AngleWithMetadata>
    wigglePermutationLinkage(cds::SearchAngles searchAngles, cds::ResidueLinkage& linkage,
                             const cds::PermutationShapePreference& shapePreference,
                             std::array<cds::ResiduesWithOverlapWeight, 2>& overlapInput)
    {
        auto& dihedrals = linkage.rotatableDihedrals;
        auto& metadata  = linkage.dihedralMetadata;
        std::vector<cds::AngleWithMetadata> shape;
        shape.resize(dihedrals.size());
        //  Reverse as convention is Glc1-4Gal and I want to wiggle in opposite direction i.e. from first
        //  rotatable bond in Asn outwards
        for (size_t rn = 0; rn < dihedrals.size(); rn++)
        {
            size_t n         = dihedrals.size() - 1 - rn;
            auto& dihedral   = dihedrals[n];
            auto preference  = cds::AngleSearchPreference {shapePreference.angles[n], shapePreference.metadataOrder[n]};
            auto coordinates = cds::dihedralCoordinates(dihedral);
            auto input       = cds::dihedralRotationInputData(dihedral, overlapInput);
            auto& metadataVector = metadata[n];
            auto best = cds::wiggleUsingRotamers(searchAngles, coordinates, codeUtils::indexVector(metadataVector),
                                                 metadataVector, preference, input);
            cds::setDihedralAngle(dihedral, best.angle);
            shape[n] = best.angle;
        }
        return shape;
    }

    std::vector<cds::AngleWithMetadata>
    wiggleConformerLinkage(cds::SearchAngles searchAngles, cds::ResidueLinkage& linkage,
                           const cds::ConformerShapePreference& shapePreference,
                           std::array<cds::ResiduesWithOverlapWeight, 2>& overlapInput)
    {
        auto& dihedrals         = linkage.rotatableDihedrals;
        auto& metadata          = linkage.dihedralMetadata;
        size_t numberOfMetadata = metadata[0].size();
        std::vector<std::vector<cds::AngleWithMetadata>> results;
        results.resize(numberOfMetadata);
        std::vector<cds::AngleOverlap> bestOverlaps;
        bestOverlaps.resize(numberOfMetadata);
        auto index = codeUtils::indexVector(metadata[0]);
        for (size_t k = 0; k < numberOfMetadata; k++)
        {
            cds::setShapeToPreference(
                linkage, cds::ConformerShapePreference {shapePreference.angles, {shapePreference.metadataOrder[k]}});
            //  Reverse as convention is Glc1-4Gal and I want to wiggle in opposite direction i.e. from first
            //  rotatable bond in Asn outwards
            for (size_t rn = 0; rn < dihedrals.size(); rn++)
            {
                size_t n       = dihedrals.size() - 1 - rn;
                auto& dihedral = dihedrals[n];
                auto preference =
                    cds::AngleSearchPreference {shapePreference.angles[n], {shapePreference.metadataOrder[k]}};
                auto coordinates = cds::dihedralCoordinates(dihedral);
                auto input       = cds::dihedralRotationInputData(dihedral, overlapInput);
                auto best = cds::wiggleUsingRotamers(searchAngles, coordinates, index, metadata[n], preference, input);
                results[k].push_back(best.angle);
                bestOverlaps[k] = best;
                cds::setDihedralAngle(dihedral, best.angle);
            }
        }
        size_t bestIndex = cds::bestOverlapResultIndex(bestOverlaps);
        std::vector<cds::AngleWithMetadata> shape;
        shape.reserve(dihedrals.size());
        for (size_t n = 0; n < dihedrals.size(); n++)
        {
            auto& bestShape = results[bestIndex][n];
            cds::setDihedralAngle(dihedrals[n], bestShape);
            shape.push_back(bestShape);
        }
        return shape;
    }

    std::vector<cds::AngleWithMetadata> wiggleLinkage(cds::SearchAngles searchAngles, cds::ResidueLinkage& linkage,
                                                      const cds::ResidueLinkageShapePreference& shapePreference,
                                                      const cds::ResiduesWithOverlapWeight& overlapResidues)
    {
        auto& nonReducing = linkage.nonReducingOverlapResidues;
        std::vector<double> nonReducingWeight(nonReducing.size(), 1.0);
        auto& reducing = linkage.reducingOverlapResidues;
        std::vector<double> reducingWeight(reducing.size(), 1.0);
        std::array<cds::ResiduesWithOverlapWeight, 2> overlapInput = {
            cds::ResiduesWithOverlapWeight {
                                            codeUtils::vectorAppend(linkage.nonReducingOverlapResidues, overlapResidues.residues),
                                            codeUtils::vectorAppend(nonReducingWeight, overlapResidues.weights)},
            cds::ResiduesWithOverlapWeight {reducing, reducingWeight}
        };
        switch (linkage.rotamerType)
        {
            case RotamerType::permutation:
                {
                    auto preference = std::get<cds::PermutationShapePreference>(shapePreference);
                    return wigglePermutationLinkage(searchAngles, linkage, preference, overlapInput);
                }
            case RotamerType::conformer:
                {
                    auto preference = std::get<cds::ConformerShapePreference>(shapePreference);
                    return wiggleConformerLinkage(searchAngles, linkage, preference, overlapInput);
                }
        }
    }

    cds::Overlap countOverlaps(const cds::ResiduesWithOverlapWeight& overlapResidues,
                               const cds::ResiduesWithOverlapWeight& glycositeResidues)
    {
        return cds::CountOverlappingAtoms(true, overlapResidues, glycositeResidues);
    }
} // namespace

cds::Overlap countTotalOverlaps(const std::vector<cds::ResiduesWithOverlapWeight>& overlapResidues,
                                const std::vector<cds::ResiduesWithOverlapWeight>& glycositeResidues)
{
    cds::Overlap overlap {0, 0.0};
    for (size_t n = 0; n < overlapResidues.size(); n++)
    {
        overlap += countOverlaps(overlapResidues[n], glycositeResidues[n]);
    }
    return overlap;
}

std::vector<size_t> determineSitesWithOverlap(uint overlapTolerance,
                                              const std::vector<cds::ResiduesWithOverlapWeight>& overlapResidues,
                                              const std::vector<cds::ResiduesWithOverlapWeight>& glycositeResidues)
{
    std::vector<size_t> indices;
    cds::Overlap overlap {0, 0.0};
    for (size_t n = 0; n < overlapResidues.size(); n++)
    {
        overlap = countOverlaps(overlapResidues[n], glycositeResidues[n]);
        if (overlap.count > overlapTolerance)
        {
            indices.push_back(n);
        }
    }
    return indices;
}

GlycoproteinState wiggleSitesWithOverlaps(pcg32& rng, cds::SearchAngles searchAngles, uint overlapTolerance,
                                          int persistCycles, bool firstLinkageOnly,
                                          std::vector<std::vector<cds::ResidueLinkage>>& glycosidicLinkages,
                                          const GlycoproteinState& initialState,
                                          const std::vector<cds::ResiduesWithOverlapWeight>& overlapResidues,
                                          const std::vector<cds::ResiduesWithOverlapWeight>& glycositeResidues)
{
    int cycle                 = 0;
    cds::Overlap overlap      = countTotalOverlaps(overlapResidues, glycositeResidues);
    auto glycositePreferences = initialState.preferences;
    auto glycositeShape       = initialState.shape;
    auto wiggledShape         = glycositeShape;
    std::vector<size_t> sites_with_overlaps =
        determineSitesWithOverlap(overlapTolerance, overlapResidues, glycositeResidues);
    while ((cycle < persistCycles) && (sites_with_overlaps.size() > 0))
    {
        ++cycle;
        gmml::log(__LINE__, __FILE__, gmml::INF,
                  "Cycle " + std::to_string(cycle) + "/" + std::to_string(persistCycles) +
                      ". Overlap: " + std::to_string(overlap.count));
        for (auto& glycosite : codeUtils::shuffleVector(rng, sites_with_overlaps))
        {
            auto& linkages = glycosidicLinkages[glycosite];
            for (size_t n = 0; (n == 0 || !firstLinkageOnly) && (n < linkages.size()); n++)
            {
                wiggledShape[glycosite][n] = wiggleLinkage(
                    searchAngles, linkages[n], glycositePreferences[glycosite][n], overlapResidues[glycosite]);
            }
        }
        cds::Overlap newOverlap = countTotalOverlaps(overlapResidues, glycositeResidues);
        if (cds::compareOverlaps(overlap, newOverlap) > 0)
        {
            overlap = newOverlap;
            cycle   = 0;
            for (auto& glycosite : sites_with_overlaps)
            {
                glycositeShape[glycosite] = wiggledShape[glycosite];
            }
        }
        else
        {
            for (auto& glycosite : sites_with_overlaps)
            {
                cds::setShape(glycosidicLinkages[glycosite], glycositeShape[glycosite]);
                wiggledShape[glycosite] = glycositeShape[glycosite];
            }
        }
        sites_with_overlaps = determineSitesWithOverlap(overlapTolerance, overlapResidues, glycositeResidues);
    }
    return {overlap, glycositePreferences, glycositeShape};
}

GlycoproteinState randomDescent(pcg32& rng, LinkageShapeRandomizer randomizeShape, cds::SearchAngles searchAngles,
                                bool monte_carlo, int persistCycles, uint overlapTolerance,
                                std::vector<std::vector<cds::ResidueLinkage>>& glycosidicLinkages,
                                const GlycoproteinState& initialState,
                                const std::vector<cds::ResiduesWithOverlapWeight>& overlapResidues,
                                const std::vector<cds::ResiduesWithOverlapWeight>& glycositeResidues)
{
    std::stringstream logss;
    logss << "Random Decent, persisting for " << persistCycles << " cycles and monte carlo is set as " << std::boolalpha
          << monte_carlo << ".\n";
    gmml::log(__LINE__, __FILE__, gmml::INF, logss.str());
    auto acceptViaMetropolisCriterion = [&rng](int difference)
    {
        double acceptance = codeUtils::uniformRandomDoubleWithinRange(rng, 0, 1);
        return monte_carlo::accept_via_metropolis_criterion(acceptance, difference);
    };

    int cycle                          = 0;
    cds::Overlap lowest_global_overlap = countTotalOverlaps(overlapResidues, glycositeResidues);
    auto glycositePreferences          = initialState.preferences;
    auto glycositeShape                = initialState.shape;
    std::vector<size_t> sites_with_overlaps =
        determineSitesWithOverlap(overlapTolerance, overlapResidues, glycositeResidues);
    while ((cycle < persistCycles) && (sites_with_overlaps.size() > 0))
    {
        gmml::log(__LINE__, __FILE__, gmml::INF,
                  "Cycle " + std::to_string(cycle) + "/" + std::to_string(persistCycles));
        ++cycle;
        for (auto& currentGlycosite : codeUtils::shuffleVector(rng, sites_with_overlaps))
        {
            auto siteResidues            = glycositeResidues[currentGlycosite];
            auto& linkages               = glycosidicLinkages[currentGlycosite];
            cds::Overlap previousOverlap = countOverlaps(overlapResidues[currentGlycosite], siteResidues);
            auto preferences             = randomizeShape(linkages);
            auto recordedShape           = glycositeShape[currentGlycosite];
            auto currentShape            = glycositeShape[currentGlycosite];
            cds::setShapeToPreference(linkages, preferences);
            for (size_t k = 0; k < 2; k++)
            {
                for (size_t n = 0; n < linkages.size(); n++)
                {
                    currentShape[n] =
                        wiggleLinkage(searchAngles, linkages[n], preferences[n], overlapResidues[currentGlycosite]);
                }
            }
            cds::Overlap newOverlap = countOverlaps(overlapResidues[currentGlycosite], siteResidues);
            double diff             = newOverlap.count - previousOverlap.count;
            bool isWorse            = cds::compareOverlaps(newOverlap, previousOverlap) > 0;
            if (isWorse || (monte_carlo && !acceptViaMetropolisCriterion(diff)))
            {
                cds::setShape(linkages, recordedShape);
            }
            else
            {
                glycositePreferences[currentGlycosite] = preferences;
                glycositeShape[currentGlycosite]       = currentShape;
                gmml::log(__LINE__, __FILE__, gmml::INF, "RandomDescent accepted a change of " + std::to_string(diff));
            }
        }
        std::vector<size_t> new_sites_with_overlaps =
            determineSitesWithOverlap(overlapTolerance, overlapResidues, glycositeResidues);
        cds::Overlap new_global_overlap = countTotalOverlaps(overlapResidues, glycositeResidues);
        if (cds::compareOverlaps(lowest_global_overlap, new_global_overlap) > 0)
        {
            lowest_global_overlap = new_global_overlap;
            cycle                 = 0;
        }
        sites_with_overlaps = new_sites_with_overlaps;
    }
    return {lowest_global_overlap, glycositePreferences, glycositeShape};
}
