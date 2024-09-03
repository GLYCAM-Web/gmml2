#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycoproteinOverlapResolution.hpp"

#include "includes/CentralDataStructure/Geometry/overlap.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralShape.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralAngleSearch.hpp"
#include "includes/CentralDataStructure/Overlaps/atomOverlaps.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CodeUtils/containers.hpp"
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
                             const std::array<cds::ResiduesWithOverlapWeight, 2>& overlapInput)
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
                           const std::array<cds::ResiduesWithOverlapWeight, 2>& overlapInput)
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
            results[k].resize(dihedrals.size());
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
                results[k][n]   = best.angle;
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
} // namespace

cds::Overlap intraGlycanOverlaps(const std::vector<cds::ResidueLinkage>& linkages)
{
    auto withWeight = [](const std::vector<Residue*> residues)
    {
        return cds::ResiduesWithOverlapWeight {residues, std::vector<double>(residues.size(), 1.0)};
    };
    cds::Overlap overlap = {0.0, 0.0};
    // skip first linkage as it connects to protein. We're only counting glycan atoms here
    for (size_t n = 1; n < linkages.size(); n++)
    {
        auto& linkage = linkages[n];
        // only take first non-reducing residue to avoid any double-counting
        overlap       += cds::CountOverlappingAtoms({{linkage.nonReducingOverlapResidues[0]}, {1.0}},
                                                    withWeight(linkage.reducingOverlapResidues));
    }
    return overlap;
}

cds::Overlap countOverlaps(const std::vector<Residue*>& overlapResidues,
                           const cds::ResiduesWithOverlapWeight& glycositeResidues)
{
    auto weights = std::vector<double>(overlapResidues.size(), 1.0);
    return cds::CountOverlappingAtoms({overlapResidues, weights}, glycositeResidues);
}

cds::Overlap siteOverlaps(OverlapWeight weight, const OverlapResidues& overlapResidues,
                          const cds::ResiduesWithOverlapWeight& glycositeResidues,
                          const std::vector<cds::ResidueLinkage>& linkages)
{
    auto proteinOverlaps = countOverlaps(overlapResidues.protein, glycositeResidues);
    auto glycanOverlaps  = countOverlaps(overlapResidues.glycan, glycositeResidues);
    auto selfOverlaps    = intraGlycanOverlaps(linkages);
    // glycans will be double counted, so halve them
    return (proteinOverlaps * weight.protein) + (selfOverlaps * weight.self) + (glycanOverlaps * weight.glycan);
}

cds::Overlap totalOverlaps(OverlapWeight weight, const std::vector<OverlapResidues>& overlapResidues,
                           const std::vector<cds::ResiduesWithOverlapWeight>& glycositeResidues,
                           const std::vector<std::vector<cds::ResidueLinkage>>& glycosidicLinkages)
{
    // glycans will be double counted, so halve their weight
    weight.glycan        *= 0.5;
    cds::Overlap overlap = {0.0, 0.0};
    for (size_t n = 0; n < overlapResidues.size(); n++)
    {
        overlap += siteOverlaps(weight, overlapResidues[n], glycositeResidues[n], glycosidicLinkages[n]);
    }
    return overlap;
}

std::vector<size_t> determineSitesWithOverlap(const std::vector<size_t>& movedSites,
                                              const std::vector<std::vector<cds::ResidueLinkage>>& glycosidicLinkages,
                                              const std::vector<OverlapResidues>& overlapResidues,
                                              const std::vector<cds::ResiduesWithOverlapWeight>& glycositeResidues)
{
    auto hasProteinOverlap = [&](size_t n)
    {
        auto overlap = countOverlaps(overlapResidues[n].protein, glycositeResidues[n]);
        return overlap.count > 0;
    };
    auto hasSelfOverlap = [&](size_t n)
    {
        auto overlap = intraGlycanOverlaps(glycosidicLinkages[n]);
        return overlap.count > 0;
    };
    std::vector<bool> justMoved(glycositeResidues.size(), false);
    for (size_t n : movedSites)
    {
        justMoved[n] = true;
    }
    std::vector<bool> glycanOverlap(glycositeResidues.size(), false);
    for (size_t n : movedSites)
    {
        for (size_t k = n + 1; k < glycositeResidues.size(); k++)
        {
            if (!(glycanOverlap[n] && glycanOverlap[k]))
            {
                auto& glycanA   = glycositeResidues[n];
                auto& glycanB   = glycositeResidues[k];
                auto& glycosite = overlapResidues[n].protein[0];
                if (cds::CountOverlappingAtoms(glycanA, glycanB).count > 0 ||
                    cds::CountOverlappingAtoms(glycanB, {{glycosite}, {1.0}}).count > 0)
                {
                    glycanOverlap[n] = true;
                    glycanOverlap[k] = true;
                }
            }
        }
    }
    std::vector<size_t> indices;
    for (size_t n = 0; n < overlapResidues.size(); n++)
    {
        // glycans which haven't moved won't overlap with protein or themselves (at least not more than before)
        if (glycanOverlap[n] || (justMoved[n] && (hasProteinOverlap(n) || hasSelfOverlap(n))))
        {
            indices.push_back(n);
        }
    }
    return indices;
}

std::vector<cds::AngleWithMetadata> wiggleLinkage(cds::SearchAngles searchAngles, cds::ResidueLinkage& linkage,
                                                  const cds::ResidueLinkageShapePreference& shapePreference,
                                                  const std::array<cds::ResiduesWithOverlapWeight, 2> overlapInput)
{
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
    throw std::runtime_error("unhandled linkage shape preference in glycoproteinOverlapResolution wiggleLinkage");
}

std::vector<std::vector<cds::AngleWithMetadata>>
wiggleGlycosite(cds::SearchAngles searchAngles, OverlapWeight weight, std::vector<cds::ResidueLinkage>& linkages,
                const std::vector<cds::ResidueLinkageShapePreference>& preferences,
                const OverlapResidues& overlapResidues)
{
    auto toWeight = [](const std::vector<cds::Residue*>& residues, double a)
    {
        return std::vector<double>(residues.size(), a);
    };
    auto baseResidues = codeUtils::vectorAppend(overlapResidues.protein, overlapResidues.glycan);
    auto baseWeights  = codeUtils::vectorAppend(toWeight(overlapResidues.protein, weight.protein),
                                                toWeight(overlapResidues.glycan, weight.glycan));
    std::vector<std::vector<cds::AngleWithMetadata>> shape;
    shape.resize(linkages.size());
    // wiggling twice gives the first linkages a second chance to resolve in a better structure
    for (size_t k = 0; k < 2; k++)
    {
        for (size_t n = 0; n < linkages.size(); n++)
        {
            auto& linkage       = linkages[n];
            auto& reducing      = linkage.reducingOverlapResidues;
            auto reducingWeight = toWeight(reducing, 1.0);

            std::vector<cds::Residue*> nonReducing = linkage.nonReducingOverlapResidues;
            // last is the glycosite, which should be present as the first entry of the protein vector
            nonReducing.pop_back();
            auto nonReducingWeight = toWeight(nonReducing, weight.self);

            std::array<cds::ResiduesWithOverlapWeight, 2> overlapInput = {
                cds::ResiduesWithOverlapWeight {codeUtils::vectorAppend(nonReducing, baseResidues),
                                                codeUtils::vectorAppend(nonReducingWeight, baseWeights)},
                cds::ResiduesWithOverlapWeight {reducing, reducingWeight}
            };
            shape[n] = wiggleLinkage(searchAngles, linkage, preferences[n], overlapInput);
        }
    }
    return shape;
}

GlycoproteinState randomDescent(pcg32 rng, LinkageShapeRandomizer randomizeShape, cds::SearchAngles searchAngles,
                                std::function<bool(int)> acceptViaMetropolisCriterion, int persistCycles,
                                OverlapWeight overlapWeight,
                                std::vector<std::vector<cds::ResidueLinkage>>& glycosidicLinkages,
                                const GlycoproteinState& initialState,
                                const std::vector<OverlapResidues>& overlapResidues,
                                const std::vector<cds::ResiduesWithOverlapWeight>& glycositeResidues)
{
    std::stringstream logss;
    logss << "Random Decent, persisting for " << persistCycles << " cycles and monte carlo is set as true.\n";
    gmml::log(__LINE__, __FILE__, gmml::INF, logss.str());

    int cycle                  = 0;
    cds::Overlap globalOverlap = totalOverlaps(overlapWeight, overlapResidues, glycositeResidues, glycosidicLinkages);
    auto glycositePreferences  = initialState.preferences;
    auto glycositeShape        = initialState.shape;
    std::vector<size_t> overlapSites = determineSitesWithOverlap(
        codeUtils::indexVector(glycosidicLinkages), glycosidicLinkages, overlapResidues, glycositeResidues);
    while ((cycle < persistCycles) && (overlapSites.size() > 0))
    {
        gmml::log(__LINE__, __FILE__, gmml::INF,
                  "Cycle " + std::to_string(cycle) + "/" + std::to_string(persistCycles));
        ++cycle;
        cds::Overlap newGlobalOverlap = globalOverlap;
        for (auto& site : codeUtils::shuffleVector(rng, overlapSites))
        {
            auto& siteResidues        = glycositeResidues[site];
            auto& siteOverlapResidues = overlapResidues[site];
            auto& linkages            = glycosidicLinkages[site];
            auto localOverlaps        = [&]()
            {
                cds::Residue* glycosite = siteOverlapResidues.protein[0];
                auto glycan             = siteOverlaps(overlapWeight, siteOverlapResidues, siteResidues, linkages);
                auto protein = countOverlaps(siteOverlapResidues.glycan, {{glycosite}, {overlapWeight.protein}});
                return glycan + protein;
            };
            auto previousOverlap = localOverlaps();
            auto preferences     = randomizeShape(linkages);
            auto recordedShape   = glycositeShape[site];
            cds::setShapeToPreference(linkages, preferences);
            auto currentShape =
                wiggleGlycosite(searchAngles, overlapWeight, linkages, preferences, siteOverlapResidues);
            auto newOverlap = localOverlaps();
            auto diff       = newOverlap + (previousOverlap * -1);
            bool isWorse    = cds::compareOverlaps(newOverlap, previousOverlap) > 0;
            if (isWorse || (!acceptViaMetropolisCriterion(diff.count)))
            {
                cds::setShape(linkages, recordedShape);
            }
            else
            {
                newGlobalOverlap           += diff;
                glycositePreferences[site] = preferences;
                glycositeShape[site]       = currentShape;
                gmml::log(__LINE__, __FILE__, gmml::INF,
                          "RandomDescent accepted a change of " + std::to_string(diff.count));
            }
        }
        if (cds::compareOverlaps(globalOverlap, newGlobalOverlap) > 0)
        {
            cycle = 0;
        }
        globalOverlap = newGlobalOverlap;
        std::vector<size_t> newOverlapSites =
            determineSitesWithOverlap(overlapSites, glycosidicLinkages, overlapResidues, glycositeResidues);
        overlapSites = newOverlapSites;
    }
    return {globalOverlap, glycositePreferences, glycositeShape};
}
