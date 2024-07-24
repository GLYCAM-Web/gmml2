#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycoproteinOverlapResolution.hpp"

#include "includes/CentralDataStructure/Overlaps/overlaps.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralShape.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralAngleSearch.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CodeUtils/metropolisCriterion.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/External_Libraries/PCG/pcg_random.h"
#include "includes/External_Libraries/PCG/pcg_extras.h"

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

std::vector<size_t> determineSitesWithOverlap(int overlapTolerance,
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

cds::Overlap wiggleSitesWithOverlaps(int overlapTolerance, int persistCycles, bool firstLinkageOnly, int interval,
                                     std::vector<std::vector<cds::ResidueLinkage>>& glycosidicLinkages,
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
        std::random_shuffle(sites_with_overlaps.begin(), sites_with_overlaps.end());
        for (auto& glycosite : sites_with_overlaps)
        {
            auto& linkages = glycosidicLinkages[glycosite];
            for (size_t n = 0; (n == 0 || !firstLinkageOnly) && (n < linkages.size()); n++)
            {
                auto& linkage                                     = linkages[n];
                std::array<std::vector<Residue*>, 2> overlapInput = {
                    codeUtils::vectorAppend(linkage.nonReducingOverlapResidues, overlapResidues[glycosite]),
                    linkage.reducingOverlapResidues};
                //  Reverse as convention is Glc1-4Gal and I want to wiggle in opposite direction i.e. from first
                //  rotatable bond in Asn outwards
                for (auto& dihedral : codeUtils::reverse(linkage.rotatableDihedrals))
                {
                    auto coordinates = cds::dihedralCoordinates(dihedral);
                    auto input       = cds::dihedralRotationInputData(dihedral, overlapInput);
                    auto best        = cds::wiggleUsingRotamers(coordinates, dihedral.metadataVector, interval, input);
                    setDihedralAngle(dihedral, best.angle);
                }
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

cds::Overlap randomDescent(bool monte_carlo, int persistCycles, int overlapTolerance,
                           std::vector<std::vector<cds::ResidueLinkage>>& glycosidicLinkages,
                           const ResiduesByType& overlapResidues,
                           const std::vector<std::vector<Residue*>>& glycositeResidues)
{
    std::stringstream logss;
    logss << "Random Decent, persisting for " << persistCycles << " cycles and monte carlo is set as " << std::boolalpha
          << monte_carlo << ".\n";
    gmml::log(__LINE__, __FILE__, gmml::INF, logss.str());
    int cycle                          = 1;
    cds::Overlap lowest_global_overlap = countTotalOverlaps(overlapResidues.all, glycositeResidues);
    std::vector<size_t> sites_with_overlaps =
        determineSitesWithOverlap(overlapTolerance, overlapResidues.all, glycositeResidues);
    if (sites_with_overlaps.size() == 0)
    {
        gmml::log(__LINE__, __FILE__, gmml::INF, "Stopping RandomDesent with all overlaps resolved.");
        return lowest_global_overlap;
    }
    // Make a random number engine for std::shuffle;
    pcg_extras::seed_seq_from<std::random_device> metropolis_seed_source;
    pcg32 rng_engine(metropolis_seed_source);
    while (cycle < persistCycles)
    {
        gmml::log(__LINE__, __FILE__, gmml::INF,
                  "Cycle " + std::to_string(cycle) + "/" + std::to_string(persistCycles));
        ++cycle;
        std::shuffle(sites_with_overlaps.begin(), sites_with_overlaps.end(), rng_engine);
        for (auto& current_glycosite : sites_with_overlaps)
        {
            auto siteResidues                     = glycositeResidues[current_glycosite];
            auto& otherGlycanResidues             = overlapResidues.glycan[current_glycosite];
            auto& otherProteinResidues            = overlapResidues.protein[current_glycosite];
            auto& linkages                        = glycosidicLinkages[current_glycosite];
            auto recordedShape                    = cds::currentShape(linkages);
            cds::Overlap previous_glycan_overlap  = countGlycositeOverlaps(otherGlycanResidues, siteResidues);
            cds::Overlap previous_protein_overlap = countGlycositeOverlaps(otherProteinResidues, siteResidues);
            cds::setRandomShapeUsingMetadata(linkages);
            // logss << "Site: " << current_glycosite->GetResidueNumber() << "\n";
            cds::Overlap new_glycan_overlap  = countGlycositeOverlaps(otherGlycanResidues, siteResidues);
            cds::Overlap new_protein_overlap = countGlycositeOverlaps(otherProteinResidues, siteResidues);
            int overlap_difference           = (new_glycan_overlap + (new_protein_overlap * 5)).count -
                                     (previous_glycan_overlap + (previous_protein_overlap * 5)).count;
            if (overlap_difference >= 0) // if the change made it worse
            {
                cds::setShape(linkages, recordedShape);
            }
            else if ((monte_carlo) && (!monte_carlo::accept_via_metropolis_criterion(overlap_difference)))
            {
                cds::setShape(linkages, recordedShape);
            }
            else
            {
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

bool dumbRandomWalk(int overlapTolerance, int maxCycles,
                    std::vector<std::vector<cds::ResidueLinkage>>& glycosidicLinkages,
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
            cds::setRandomShapeUsingMetadata(glycosidicLinkages[currentGlycosite]);
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

void resolveOverlapsWithWiggler(bool randomize, int persistCycles, int overlapTolerance,
                                std::vector<std::vector<cds::ResidueLinkage>>& glycosidicLinkages,
                                const ResiduesByType& overlapResidues,
                                const std::vector<std::vector<Residue*>>& glycositeResidues)
{
    if (randomize)
    { // First try a very fast/cheap approach
        if (dumbRandomWalk(overlapTolerance, 10, glycosidicLinkages, overlapResidues.all,
                           glycositeResidues)) // returns true if it fully resolves overlaps.
        {
            return;
        }
    }
    bool useMonteCarlo          = true;
    bool wiggleFirstLinkageOnly = true; // Only happens when passed to Wiggle function.
    int angleIncrement          = 5;
    cds::Overlap currentOverlap =
        wiggleSitesWithOverlaps(overlapTolerance, persistCycles, wiggleFirstLinkageOnly, angleIncrement,
                                glycosidicLinkages, overlapResidues.all, glycositeResidues);
    gmml::log(__LINE__, __FILE__, gmml::INF, "1. Overlap: " + std::to_string(currentOverlap.count));
    if (randomize)
    {
        currentOverlap = randomDescent(useMonteCarlo, persistCycles, overlapTolerance, glycosidicLinkages,
                                       overlapResidues, glycositeResidues);
        gmml::log(__LINE__, __FILE__, gmml::INF, "1r. Overlap: " + std::to_string(currentOverlap.count));
    }
    currentOverlap = wiggleSitesWithOverlaps(overlapTolerance, persistCycles, false, angleIncrement, glycosidicLinkages,
                                             overlapResidues.all, glycositeResidues);
    gmml::log(__LINE__, __FILE__, gmml::INF, "2. Overlap: " + std::to_string(currentOverlap.count));
    currentOverlap = wiggleSitesWithOverlaps(overlapTolerance, persistCycles, wiggleFirstLinkageOnly, angleIncrement,
                                             glycosidicLinkages, overlapResidues.all, glycositeResidues);
    gmml::log(__LINE__, __FILE__, gmml::INF, "3. Overlap: " + std::to_string(currentOverlap.count));
}
