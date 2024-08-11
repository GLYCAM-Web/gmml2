#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycoproteinOverlapResolution.hpp"

#include "includes/CentralDataStructure/Overlaps/overlaps.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralShape.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralAngleSearch.hpp"
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

    void wigglePermutationLinkage(double interval, cds::ResidueLinkage& linkage,
                                  const cds::PermutationShapePreference& shapePreference,
                                  std::array<std::vector<Residue*>, 2>& overlapInput)
    {
        auto& dihedrals = linkage.rotatableDihedrals;
        auto& metadata  = linkage.dihedralMetadata;
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
            auto best = cds::wiggleUsingRotamers(coordinates, codeUtils::indexVector(metadataVector), metadataVector,
                                                 preference, interval, input);
            cds::setDihedralAngle(dihedral, best.angle);
        }
    }

    void wiggleConformerLinkage(double interval, cds::ResidueLinkage& linkage,
                                const cds::ConformerShapePreference& shapePreference,
                                std::array<std::vector<Residue*>, 2>& overlapInput)
    {
        auto& dihedrals         = linkage.rotatableDihedrals;
        auto& metadata          = linkage.dihedralMetadata;
        size_t numberOfMetadata = metadata[0].size();
        std::vector<std::vector<cds::AngleWithMetadata>> results;
        results.resize(numberOfMetadata);
        std::vector<cds::AngleOverlap> bestOverlaps;
        bestOverlaps.resize(numberOfMetadata);
        auto initialShape = cds::currentShape(dihedrals, metadata);
        auto index        = codeUtils::indexVector(metadata[0]);
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
                auto best = cds::wiggleUsingRotamers(coordinates, index, metadata[n], preference, interval, input);
                results[k].push_back(best.angle);
                bestOverlaps[k] = best;
                cds::setDihedralAngle(dihedral, best.angle);
            }
        }
        size_t bestIndex = cds::bestOverlapResultIndex(bestOverlaps);
        for (size_t n = 0; n < dihedrals.size(); n++)
        {
            cds::setDihedralAngle(dihedrals[n], results[bestIndex][n]);
        }
    }

    void wiggleLinkage(double interval, cds::ResidueLinkage& linkage,
                       const cds::ResidueLinkageShapePreference& shapePreference,
                       const std::vector<Residue*>& overlapResidues)
    {
        std::array<std::vector<Residue*>, 2> overlapInput = {
            codeUtils::vectorAppend(linkage.nonReducingOverlapResidues, overlapResidues),
            linkage.reducingOverlapResidues};
        switch (linkage.rotamerType)
        {
            case RotamerType::permutation:
                {
                    auto preference = std::get<cds::PermutationShapePreference>(shapePreference);
                    wigglePermutationLinkage(interval, linkage, preference, overlapInput);
                    return;
                }
            case RotamerType::conformer:
                {
                    auto preference = std::get<cds::ConformerShapePreference>(shapePreference);
                    wiggleConformerLinkage(interval, linkage, preference, overlapInput);
                    return;
                }
        }
    }
} // namespace

cds::Overlap countGlycositeOverlaps(const std::vector<Residue*>& overlapResidues,
                                    const std::vector<Residue*>& glycositeResidues)
{
    return cds::CountOverlappingAtoms(overlapResidues, glycositeResidues);
}

cds::Overlap countTotalOverlaps(const std::vector<std::vector<Residue*>>& overlapResidues,
                                const std::vector<std::vector<Residue*>>& glycositeResidues)
{
    cds::Overlap overlap {0, 0.0};
    for (size_t n = 0; n < overlapResidues.size(); n++)
    {
        overlap += countGlycositeOverlaps(overlapResidues[n], glycositeResidues[n]);
    }
    return overlap;
}

std::vector<size_t> determineSitesWithOverlap(uint overlapTolerance,
                                              const std::vector<std::vector<Residue*>>& overlapResidues,
                                              const std::vector<std::vector<Residue*>>& glycositeResidues)
{
    std::vector<size_t> indices;
    cds::Overlap overlap {0, 0.0};
    for (size_t n = 0; n < overlapResidues.size(); n++)
    {
        overlap = countGlycositeOverlaps(overlapResidues[n], glycositeResidues[n]);
        if (overlap.count > overlapTolerance)
        {
            indices.push_back(n);
        }
    }
    return indices;
}

cds::Overlap
wiggleSitesWithOverlaps(pcg32& rng, uint overlapTolerance, int persistCycles, bool firstLinkageOnly, double interval,
                        std::vector<std::vector<cds::ResidueLinkage>>& glycosidicLinkages,
                        const std::vector<std::vector<cds::ResidueLinkageShapePreference>>& glycositePreferences,
                        const std::vector<std::vector<Residue*>>& overlapResidues,
                        const std::vector<std::vector<Residue*>>& glycositeResidues)
{
    int cycle            = 0;
    cds::Overlap overlap = countTotalOverlaps(overlapResidues, glycositeResidues);
    while (cycle < persistCycles)
    {
        ++cycle;
        gmml::log(__LINE__, __FILE__, gmml::INF,
                  "Cycle " + std::to_string(cycle) + "/" + std::to_string(persistCycles) +
                      ". Overlap: " + std::to_string(overlap.count));
        std::vector<size_t> sites_with_overlaps =
            determineSitesWithOverlap(overlapTolerance, overlapResidues, glycositeResidues);
        if (sites_with_overlaps.size() == 0)
        {
            return overlap;
        }
        for (auto& glycosite : codeUtils::shuffleVector(rng, sites_with_overlaps))
        {
            auto& linkages = glycosidicLinkages[glycosite];
            for (size_t n = 0; (n == 0 || !firstLinkageOnly) && (n < linkages.size()); n++)
            {
                wiggleLinkage(interval, linkages[n], glycositePreferences[glycosite][n], overlapResidues[glycosite]);
            }
        }
        cds::Overlap newOverlap = countTotalOverlaps(overlapResidues, glycositeResidues);
        if (cds::compareOverlaps(overlap, newOverlap) > 0)
        {
            overlap = newOverlap;
            cycle   = 1;
        }
    }
    return overlap;
}

cds::Overlap randomDescent(pcg32& rng, LinkageShapeRandomizer randomizeShape, bool monte_carlo, int persistCycles,
                           uint overlapTolerance, std::vector<std::vector<cds::ResidueLinkage>>& glycosidicLinkages,
                           std::vector<std::vector<cds::ResidueLinkageShapePreference>>& glycositePreferences,
                           const ResiduesByType& overlapResidues,
                           const std::vector<std::vector<Residue*>>& glycositeResidues)
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

    int cycle                          = 1;
    cds::Overlap lowest_global_overlap = countTotalOverlaps(overlapResidues.all, glycositeResidues);
    std::vector<size_t> sites_with_overlaps =
        determineSitesWithOverlap(overlapTolerance, overlapResidues.all, glycositeResidues);
    if (sites_with_overlaps.size() == 0)
    {
        gmml::log(__LINE__, __FILE__, gmml::INF, "Stopping RandomDesent with all overlaps resolved.");
        return lowest_global_overlap;
    }
    while (cycle < persistCycles)
    {
        gmml::log(__LINE__, __FILE__, gmml::INF,
                  "Cycle " + std::to_string(cycle) + "/" + std::to_string(persistCycles));
        ++cycle;
        for (auto& currentGlycosite : codeUtils::shuffleVector(rng, sites_with_overlaps))
        {
            auto siteResidues                     = glycositeResidues[currentGlycosite];
            auto& otherGlycanResidues             = overlapResidues.glycan[currentGlycosite];
            auto& otherProteinResidues            = overlapResidues.protein[currentGlycosite];
            auto& linkages                        = glycosidicLinkages[currentGlycosite];
            auto recordedShape                    = cds::currentShape(linkages);
            cds::Overlap previous_glycan_overlap  = countGlycositeOverlaps(otherGlycanResidues, siteResidues);
            cds::Overlap previous_protein_overlap = countGlycositeOverlaps(otherProteinResidues, siteResidues);
            auto preferences                      = randomizeShape(linkages);
            cds::setShapeToPreference(linkages, preferences);
            cds::Overlap new_glycan_overlap  = countGlycositeOverlaps(otherGlycanResidues, siteResidues);
            cds::Overlap new_protein_overlap = countGlycositeOverlaps(otherProteinResidues, siteResidues);
            int overlap_difference           = (new_glycan_overlap + (new_protein_overlap * 5)).count -
                                     (previous_glycan_overlap + (previous_protein_overlap * 5)).count;
            if (overlap_difference >= 0) // if the change made it worse
            {
                cds::setShape(linkages, recordedShape);
            }
            else if (monte_carlo && !acceptViaMetropolisCriterion(overlap_difference))
            {
                cds::setShape(linkages, recordedShape);
            }
            else
            {
                glycositePreferences[currentGlycosite] = preferences;
                gmml::log(__LINE__, __FILE__, gmml::INF,
                          "RandomDescent accepted a change of " + std::to_string(overlap_difference));
            }
        }
        cds::Overlap new_global_overlap = countTotalOverlaps(overlapResidues.all, glycositeResidues);
        sites_with_overlaps = determineSitesWithOverlap(overlapTolerance, overlapResidues.all, glycositeResidues);
        if (sites_with_overlaps.size() == 0)
        {
            gmml::log(__LINE__, __FILE__, gmml::INF, "Stopping RandomDesent with all overlaps resolved.");
            return new_global_overlap;
        }
        if (cds::compareOverlaps(lowest_global_overlap, new_global_overlap) > 0)
        {
            lowest_global_overlap = new_global_overlap;
            cycle                 = 1;
        }
    }
    return lowest_global_overlap;
}

bool dumbRandomWalk(LinkageShapeRandomizer randomizeShape, uint overlapTolerance, int maxCycles,
                    std::vector<std::vector<cds::ResidueLinkage>>& glycosidicLinkages,
                    std::vector<std::vector<cds::ResidueLinkageShapePreference>>& glycositePreferences,
                    const std::vector<std::vector<Residue*>>& overlapResidues,
                    const std::vector<std::vector<Residue*>>& glycositeResidues)
{
    gmml::log(__LINE__, __FILE__, gmml::INF, "Starting DumbRandomWalk.");
    int cycle = 1;
    std::vector<size_t> sites_with_overlaps =
        determineSitesWithOverlap(overlapTolerance, overlapResidues, glycositeResidues);
    while (cycle < maxCycles)
    {
        ++cycle;
        for (auto& currentGlycosite : sites_with_overlaps)
        {
            auto& linkages   = glycosidicLinkages[currentGlycosite];
            auto preferences = randomizeShape(linkages);
            cds::setShapeToPreference(linkages, preferences);
            glycositePreferences[currentGlycosite] = preferences;
        }
        gmml::log(__LINE__, __FILE__, gmml::INF, "Updating list of sites with overlaps.");
        sites_with_overlaps = determineSitesWithOverlap(overlapTolerance, overlapResidues, glycositeResidues);
        if (sites_with_overlaps.empty())
        {
            gmml::log(__LINE__, __FILE__, gmml::INF, "DumbRandomWalk resolved the overlaps. Stopping.");
            return true;
        }
    }
    gmml::log(__LINE__, __FILE__, gmml::INF, "DumbRandomWalk did not resolve the overlaps.");
    return false;
}
