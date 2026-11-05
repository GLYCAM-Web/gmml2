#include "include/programs/GlycoproteinBuilder/overlapResolution.hpp"

#include "include/CentralDataStructure/graphInterface.hpp"
#include "include/assembly/assemblyBounds.hpp"
#include "include/assembly/assemblyGraph.hpp"
#include "include/assembly/assemblyIndices.hpp"
#include "include/assembly/assemblyOverlap.hpp"
#include "include/assembly/assemblySelection.hpp"
#include "include/carbohydrate/dihedralAngleSearch.hpp"
#include "include/carbohydrate/dihedralShape.hpp"
#include "include/carbohydrate/pdbWriter.hpp"
#include "include/external/pcg/pcg_random.h"
#include "include/fileType/off/offFileData.hpp"
#include "include/fileType/off/offFileWriter.hpp"
#include "include/fileType/pdb/pdbFile.hpp"
#include "include/fileType/pdb/pdbFileData.hpp"
#include "include/fileType/pdb/pdbFileWriter.hpp"
#include "include/fileType/pdb/pdbResidue.hpp"
#include "include/geometry/boundingSphere.hpp"
#include "include/geometry/geometryTypes.hpp"
#include "include/geometry/orientation.hpp"
#include "include/geometry/overlap.hpp"
#include "include/glycoprotein/glycanShape.hpp"
#include "include/glycoprotein/glycanWiggle.hpp"
#include "include/glycoprotein/glycoproteinStructs.hpp"
#include "include/glycoprotein/glycoproteinUtil.hpp"
#include "include/glycoprotein/overlapCount.hpp"
#include "include/glycoprotein/shapeRandomization.hpp"
#include "include/glycoprotein/sidechains.hpp"
#include "include/glycoprotein/writerInterface.hpp"
#include "include/graph/graphManipulation.hpp"
#include "include/graph/graphTypes.hpp"
#include "include/metadata/dihedralangledata.hpp"
#include "include/programs/GlycoproteinBuilder/gpBuilderStats.hpp"
#include "include/programs/GlycoproteinBuilder/randomDescent.hpp"
#include "include/util/constants.hpp"
#include "include/util/containers.hpp"
#include "include/util/files.hpp"
#include "include/util/filesystem.hpp"
#include "include/util/logging.hpp"
#include "include/util/random.hpp"
#include "include/util/strings.hpp"
#include "include/util/structuredFiles.hpp"
#include "include/util/threads.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <optional>
#include <ostream>
#include <sstream>
#include <stddef.h> // size_t
#include <stdexcept>
#include <string>
#include <utility>
#include <variant>
#include <vector>

namespace gmml
{
    namespace gpbuilder
    {
        namespace
        {
            struct OmegaConformerSite
            {
                size_t glycanId;
                size_t linkageIndex;
                size_t bondIndex;
                size_t ggMetadataIndex;
                size_t gtMetadataIndex;
                std::string label;
            };

            bool isAlphaC1Phi(const DihedralAngleData& entry)
            {
                return entry.dihedral_angle_name_ == "Phi" && entry.linking_atom1_ == "C1" &&
                       util::contains(entry.residue1_conditions_, std::string("alpha"));
            }

            // Each alpha1-6 linkage whose omega metadata contains both gg and gt.
            std::vector<OmegaConformerSite> alpha16OmegaSites(
                const DihedralAngleDataTable& table, const AssemblyData& data)
            {
                std::vector<OmegaConformerSite> sites;
                const size_t glycanCount = data.glycans.linkages.size();
                for (size_t glycanId = 0; glycanId < glycanCount; glycanId++)
                {
                    const std::vector<size_t>& linkages = data.glycans.linkages[glycanId];
                    for (size_t linkageIndex = 0; linkageIndex < linkages.size(); linkageIndex++)
                    {
                        size_t linkageId = linkages[linkageIndex];
                        if (data.residueLinkages.rotamerTypes[linkageId] != RotamerType::permutation)
                        {
                            continue;
                        }
                        const std::vector<size_t>& bonds = data.residueLinkages.rotatableBonds[linkageId];
                        bool alphaC1 = false;
                        for (size_t bondId : bonds)
                        {
                            for (size_t metadata : data.rotatableBonds.metadata[bondId])
                            {
                                if (isAlphaC1Phi(table.entries[metadata]))
                                {
                                    alphaC1 = true;
                                }
                            }
                        }
                        if (!alphaC1)
                        {
                            continue;
                        }
                        for (size_t bondIndex = 0; bondIndex < bonds.size(); bondIndex++)
                        {
                            const std::vector<size_t>& metadata = data.rotatableBonds.metadata[bonds[bondIndex]];
                            std::optional<size_t> gg;
                            std::optional<size_t> gt;
                            for (size_t n = 0; n < metadata.size(); n++)
                            {
                                const DihedralAngleData& entry = table.entries[metadata[n]];
                                if (entry.dihedral_angle_name_ != "Omg" || entry.linking_atom2_ != "O6")
                                {
                                    continue;
                                }
                                if (entry.rotamer_name_ == "gg")
                                {
                                    gg = n;
                                }
                                else if (entry.rotamer_name_ == "gt")
                                {
                                    gt = n;
                                }
                            }
                            if (gg.has_value() && gt.has_value())
                            {
                                sites.push_back(
                                    OmegaConformerSite {glycanId, linkageIndex, bondIndex, gg.value(), gt.value(), ""});
                            }
                        }
                    }
                }
                std::vector<size_t> perGlycan(glycanCount, 0);
                for (const OmegaConformerSite& site : sites)
                {
                    perGlycan[site.glycanId] += 1;
                }
                for (OmegaConformerSite& site : sites)
                {
                    size_t residue = data.glycans.attachmentResidue[site.glycanId];
                    const std::string& chain = data.residues.chainIds[residue];
                    site.label =
                        (chain.empty() ? std::string() : chain + "_") + std::to_string(data.residues.numbers[residue]);
                    if (perGlycan[site.glycanId] > 1)
                    {
                        site.label += "-L" + std::to_string(site.linkageIndex);
                    }
                }
                return sites;
            }

            size_t omegaBondId(const AssemblyData& data, const OmegaConformerSite& site)
            {
                size_t linkageId = data.glycans.linkages[site.glycanId][site.linkageIndex];
                return data.residueLinkages.rotatableBonds[linkageId][site.bondIndex];
            }

            double circularDistanceDegrees(double a, double b)
            {
                double distance = std::fmod(std::abs(a - b), 360.0);
                return distance > 180.0 ? 360.0 - distance : distance;
            }

            // Nearest of the three omega wells. Sampling can leave the canonical angle, including tg.
            std::string omegaRotamerName(double degrees)
            {
                const std::array<std::pair<const char*, double>, 3> rotamers {
                    {
                     {"gg", -60.0},
                     {"gt", 60.0},
                     {"tg", 180.0},
                     }
                };
                const char* bestName = rotamers[0].first;
                double bestDistance = circularDistanceDegrees(degrees, rotamers[0].second);
                for (size_t n = 1; n < rotamers.size(); n++)
                {
                    double distance = circularDistanceDegrees(degrees, rotamers[n].second);
                    if (distance < bestDistance)
                    {
                        bestName = rotamers[n].first;
                        bestDistance = distance;
                    }
                }
                return bestName;
            }

            double currentOmegaDegrees(
                const AssemblyData& data, const assembly::Bounds& bounds, const OmegaConformerSite& site)
            {
                return constants::toDegrees(angle(dihedralCoordinates(data, bounds, omegaBondId(data, site))));
            }

            // Move only this alpha1-6 omega inside the chosen rotamer. The search settings select the
            // range around that rotamer. Every other dihedral stays put.
            void wiggleSingleOmega(
                const DihedralAngleDataTable& dihedralAngleTable,
                const OverlapSettings& overlapSettings,
                const assembly::Graph& graph,
                const AssemblyData& data,
                const assembly::Selection& selection,
                MutableData& mutableData,
                const OmegaConformerSite& site,
                bool useGt,
                const AngleSearchSettings& searchSettings)
            {
                size_t bondId = omegaBondId(data, site);
                const std::vector<size_t>& metadata = data.rotatableBonds.metadata[bondId];
                size_t chosen = useGt ? site.gtMetadataIndex : site.ggMetadataIndex;
                std::vector<double> angles(metadata.size(), 0.0);
                angles[chosen] = dihedralAngleTable.entries[metadata[chosen]].default_angle;
                AngleSearchPreference preference {searchSettings.deviation, angles, {chosen}};
                const std::array<Coordinate, 4> coordinates = dihedralCoordinates(data, mutableData.bounds, bondId);
                const std::vector<size_t>& movingAtoms = data.rotatableBonds.movingAtoms[bondId];
                std::vector<bool> atomMoving = util::indicesToBools(atomCount(graph.source), movingAtoms);
                assembly::Selection moving =
                    assembly::intersection(graph, selection, assembly::selectByAtoms(graph, atomMoving));
                assembly::Selection nonMoving = assembly::intersection(
                    graph, selection, assembly::selectByAtoms(graph, util::vectorNot(atomMoving)));
                auto searchOverlap = [&](const assembly::Bounds& bounds)
                {
                    return overlapAboveThresholdSum(
                        overlapSettings.rejectionThreshold,
                        overlapsBetweenSelections(
                            overlapSettings.potentialTable,
                            overlapSettings.tolerance,
                            graph,
                            bounds,
                            moving,
                            nonMoving,
                            data.atoms.elements,
                            data.residueEdges.atomsCloseToEdge));
                };
                OverlapState best = wiggleUsingRotamers(
                    searchOverlap,
                    searchSettings.angles,
                    searchSettings.halfIntervalSearches,
                    dihedralAngleTable,
                    graph,
                    mutableData.bounds,
                    movingAtoms,
                    coordinates,
                    util::indexVector(metadata),
                    metadata,
                    preference);
                mutableData.bounds = best.bounds;
                mutableData.rotatableBondCurrentMetadata[bondId] = best.angle.metadataIndex;
            }

            std::string omegaCombinationPrefix(
                const std::vector<OmegaConformerSite>& sites, size_t mask, const std::string& base)
            {
                std::vector<std::string> parts;
                parts.reserve(sites.size());
                for (size_t n = 0; n < sites.size(); n++)
                {
                    bool useGt = (mask & (size_t(1) << n)) != 0;
                    parts.push_back(sites[n].label + "-" + (useGt ? "gt" : "gg"));
                }
                return base + "_" + util::join("_", parts);
            }

            // Linkage indices, within the glycan, strictly beyond this alpha1-6 linkage.
            std::vector<size_t> linkagesAfterOmega(
                const assembly::Graph& graph, const AssemblyData& data, const OmegaConformerSite& site)
            {
                size_t glycanId = site.glycanId;
                size_t linkageId = data.glycans.linkages[glycanId][site.linkageIndex];
                size_t edgeId = data.residueLinkages.edgeId[linkageId];
                const std::array<size_t, 2>& ends = graph.residues.edges.nodeAdjacencies[edgeId];
                size_t attachment = data.glycans.attachmentResidue[glycanId];
                const std::vector<size_t>& glycanResidues =
                    assembly::moleculeResidues(graph, data.glycans.moleculeId[glycanId]);
                const size_t residueCount = graph::nodeCount(graph.residues);
                std::vector<bool> inGlycan(residueCount, false);
                for (size_t residue : glycanResidues)
                {
                    inGlycan[residue] = true;
                }
                inGlycan[attachment] = true;

                std::vector<size_t> distance(residueCount, std::numeric_limits<size_t>::max());
                std::vector<size_t> queue {attachment};
                distance[attachment] = 0;
                for (size_t q = 0; q < queue.size(); q++)
                {
                    size_t residue = queue[q];
                    for (size_t neighbor : graph.residues.nodes.nodeAdjacencies[residue])
                    {
                        if (inGlycan[neighbor] && distance[neighbor] == std::numeric_limits<size_t>::max())
                        {
                            distance[neighbor] = distance[residue] + 1;
                            queue.push_back(neighbor);
                        }
                    }
                }
                size_t closer = distance[ends[0]] <= distance[ends[1]] ? ends[0] : ends[1];
                size_t farther = closer == ends[0] ? ends[1] : ends[0];

                std::vector<bool> inSubtree(residueCount, false);
                queue.clear();
                queue.push_back(farther);
                inSubtree[farther] = true;
                for (size_t q = 0; q < queue.size(); q++)
                {
                    size_t residue = queue[q];
                    const std::vector<size_t>& neighbors = graph.residues.nodes.nodeAdjacencies[residue];
                    const std::vector<size_t>& edges = graph.residues.nodes.edgeAdjacencies[residue];
                    for (size_t n = 0; n < neighbors.size(); n++)
                    {
                        if (edges[n] == edgeId || !inGlycan[neighbors[n]] || inSubtree[neighbors[n]])
                        {
                            continue;
                        }
                        inSubtree[neighbors[n]] = true;
                        queue.push_back(neighbors[n]);
                    }
                }

                std::vector<size_t> result;
                const std::vector<size_t>& linkages = data.glycans.linkages[glycanId];
                for (size_t n = 0; n < linkages.size(); n++)
                {
                    if (n == site.linkageIndex)
                    {
                        continue;
                    }
                    const std::array<size_t, 2>& linkageEnds =
                        graph.residues.edges.nodeAdjacencies[data.residueLinkages.edgeId[linkages[n]]];
                    if (inSubtree[linkageEnds[0]] && inSubtree[linkageEnds[1]])
                    {
                        result.push_back(n);
                    }
                }
                return result;
            }

            bool linkageContainsChi(const DihedralAngleDataTable& table, const AssemblyData& data, size_t linkageId)
            {
                for (size_t bondId : data.residueLinkages.rotatableBonds[linkageId])
                {
                    const std::vector<size_t>& metadata = data.rotatableBonds.metadata[bondId];
                    if (metadata.empty())
                    {
                        continue;
                    }
                    const std::string& name = table.entries[metadata[0]].dihedral_angle_name_;
                    if (name == "Chi1" || name == "Chi2")
                    {
                        return true;
                    }
                }
                return false;
            }

            // Linkage index, within the glycan, of the linkage immediately closer to the attachment than
            // this alpha1-6. Empty when there is no such linkage, or when it carries Chi1 or Chi2.
            std::optional<size_t> linkageBeforeOmega(
                const DihedralAngleDataTable& table,
                const assembly::Graph& graph,
                const AssemblyData& data,
                const OmegaConformerSite& site)
            {
                size_t glycanId = site.glycanId;
                size_t linkageId = data.glycans.linkages[glycanId][site.linkageIndex];
                size_t edgeId = data.residueLinkages.edgeId[linkageId];
                const std::array<size_t, 2>& ends = graph.residues.edges.nodeAdjacencies[edgeId];
                size_t attachment = data.glycans.attachmentResidue[glycanId];
                const std::vector<size_t>& glycanResidues =
                    assembly::moleculeResidues(graph, data.glycans.moleculeId[glycanId]);
                const size_t residueCount = graph::nodeCount(graph.residues);
                std::vector<bool> inGlycan(residueCount, false);
                for (size_t residue : glycanResidues)
                {
                    inGlycan[residue] = true;
                }
                inGlycan[attachment] = true;

                std::vector<size_t> distance(residueCount, std::numeric_limits<size_t>::max());
                std::vector<size_t> queue {attachment};
                distance[attachment] = 0;
                for (size_t q = 0; q < queue.size(); q++)
                {
                    size_t residue = queue[q];
                    for (size_t neighbor : graph.residues.nodes.nodeAdjacencies[residue])
                    {
                        if (inGlycan[neighbor] && distance[neighbor] == std::numeric_limits<size_t>::max())
                        {
                            distance[neighbor] = distance[residue] + 1;
                            queue.push_back(neighbor);
                        }
                    }
                }
                size_t closer = distance[ends[0]] <= distance[ends[1]] ? ends[0] : ends[1];
                if (closer == attachment || distance[closer] == std::numeric_limits<size_t>::max())
                {
                    return std::nullopt;
                }
                size_t parent = std::numeric_limits<size_t>::max();
                for (size_t neighbor : graph.residues.nodes.nodeAdjacencies[closer])
                {
                    if (inGlycan[neighbor] && distance[neighbor] + 1 == distance[closer])
                    {
                        parent = neighbor;
                        break;
                    }
                }
                if (parent == std::numeric_limits<size_t>::max())
                {
                    return std::nullopt;
                }
                const std::vector<size_t>& linkages = data.glycans.linkages[glycanId];
                for (size_t n = 0; n < linkages.size(); n++)
                {
                    if (n == site.linkageIndex)
                    {
                        continue;
                    }
                    const std::array<size_t, 2>& linkageEnds =
                        graph.residues.edges.nodeAdjacencies[data.residueLinkages.edgeId[linkages[n]]];
                    bool connectsParent = (linkageEnds[0] == closer && linkageEnds[1] == parent) ||
                                          (linkageEnds[1] == closer && linkageEnds[0] == parent);
                    if (!connectsParent)
                    {
                        continue;
                    }
                    if (linkageContainsChi(table, data, linkages[n]))
                    {
                        return std::nullopt;
                    }
                    return n;
                }
                return std::nullopt;
            }
        } // namespace

        void resolveOverlaps(
            std::vector<StructureStats>& stats,
            const SidechainRotamerData& sidechainRotamers,
            const DihedralAngleDataTable& dihedralAngleDataTable,
            const GlycoproteinAssembly& assembly,
            const OverlapSettings& overlapSettings,
            const ResolutionSettings& resolutionSettings,
            const std::string& outputDir,
            const std::vector<std::string>& headerLines,
            int numThreads)
        {
            const std::string structureDir = outputDir + "/samples";
            const std::string rejectDir = structureDir + "/rejected";
            const std::string deletionDir = structureDir + "/deletions";

            MetadataOrder initialMetadataOrder =
                [](pcg32&, const DihedralAngleDataTable&, const std::vector<size_t>& metadataVector, bool allowFallback)
            { return allowFallback ? util::indexVector(metadataVector) : std::vector<size_t> {0}; };
            // only picks out the first of randomly ordered metadata
            MetadataOrder randomMetadataSelection = [&resolutionSettings](
                                                        pcg32& rng,
                                                        const DihedralAngleDataTable& table,
                                                        const std::vector<size_t>& metadataIndices,
                                                        bool allowFallback)
            {
                std::vector<size_t> reordered =
                    util::weightedRandomOrder(rng, util::indicesToValues(table.weights, metadataIndices));
                return allowFallback ? reordered : std::vector<size_t> {reordered[0]};
            };
            bool freezeGlycositeResidueConformation = resolutionSettings.useInitialGlycositeResidueConformation;
            bool allowRotamerFallback = resolutionSettings.allowRotamerFallback;
            bool moveOverlappingSidechains = resolutionSettings.moveOverlappingSidechains;
            auto standardDeviation = [](const AngleSettings& settings, const DihedralAngleData& metadata)
            {
                std::function<std::pair<double, double>(const AngleLimit&)> onLimit = [&](const AngleLimit& dev)
                {
                    double max_std = settings.preferenceDeviation + settings.searchDeviation;
                    double lower_std = dev.lowerDeviationLimit / max_std;
                    double upper_std = dev.upperDeviationLimit / max_std;
                    return std::pair<double, double> {lower_std, upper_std};
                };
                std::function<std::pair<double, double>(const AngleStd&)> onStd = [&](const AngleStd& dev) {
                    return std::pair<double, double> {dev.lowerDeviationStd, dev.upperDeviationStd};
                };
                return onAngleDeviation(onLimit, onStd, metadata.angle_deviation);
            };
            auto randomAngle =
                [&standardDeviation](pcg32& rng, const AngleSettings& settings, const DihedralAngleData& metadata)
            {
                double stdCutoff = settings.preferenceDeviation;
                double num = util::normalDistributionRandomDoubleWithCutoff(rng, -stdCutoff, stdCutoff);
                auto std = standardDeviation(settings, metadata);
                return metadata.default_angle + num * (num < 0 ? std.first : std.second);
            };
            GlycanShapeRandomizer randomizeShape =
                [&dihedralAngleDataTable, &randomAngle, &freezeGlycositeResidueConformation, &allowRotamerFallback](
                    pcg32& rng,
                    const AngleSettings& settings,
                    const AssemblyData& data,
                    const assembly::Bounds& bounds,
                    const LinkageShapeSettings& shapeSettings,
                    const std::vector<GlycanShapePreference>& currentPreference,
                    size_t glycanId)
            {
                GlycanShapePreference result = randomLinkageShapePreference(
                    rng,
                    dihedralAngleDataTable,
                    settings,
                    data,
                    bounds,
                    shapeSettings,
                    glycanId,
                    randomAngle,
                    freezeGlycositeResidueConformation);
                const std::vector<bool>& retain = shapeSettings.retainInitialPreference[glycanId];
                for (size_t n = 0; n < retain.size(); n++)
                {
                    if (retain[n])
                    {
                        result[n] = currentPreference[glycanId][n];
                    }
                }
                return result;
            };

            // Default-structure sampling for a full glycan (all linkages): preferenceDev = 0,
            // initialMetadataOrder, no retain, allow rotamer fallback on every linkage.
            const AngleSettings defaultGlycositeAngleSettings {0.0, 3.0, 3.0, 2, initialMetadataOrder};
            FullGlycanDihedralRandomizer randomizeAllDihedrals =
                [&dihedralAngleDataTable,
                 &randomAngle,
                 &defaultGlycositeAngleSettings,
                 &freezeGlycositeResidueConformation](
                    pcg32& rng, const AssemblyData& data, const assembly::Bounds& bounds, size_t glycanId)
            {
                const size_t glycanCount = data.glycans.moleculeId.size();
                std::vector<std::vector<bool>> retain(glycanCount);
                std::vector<std::vector<bool>> allowFallback(glycanCount);
                for (size_t g = 0; g < glycanCount; g++)
                {
                    const size_t linkageCount = data.glycans.linkages[g].size();
                    retain[g] = std::vector<bool>(linkageCount, false);
                    allowFallback[g] = std::vector<bool>(linkageCount, true);
                }
                const LinkageShapeSettings fullRandomSettings {retain, allowFallback};
                return randomLinkageShapePreference(
                    rng,
                    dihedralAngleDataTable,
                    defaultGlycositeAngleSettings,
                    data,
                    bounds,
                    fullRandomSettings,
                    glycanId,
                    randomAngle,
                    freezeGlycositeResidueConformation);
            };

            auto getCoordinates = [](const std::vector<Sphere>& bounds)
            {
                std::vector<Coordinate> result;
                result.reserve(bounds.size());
                for (size_t n = 0; n < bounds.size(); n++)
                {
                    result.push_back(bounds[n].center);
                }
                return result;
            };

            auto sidechainRotationPreferences =
                [](pcg32& rng, const AssemblyData& data, const std::vector<size_t>& residues)
            {
                std::vector<std::vector<size_t>> result;
                result.reserve(residues.size());
                for (size_t residue : residues)
                {
                    result.push_back(util::weightedRandomOrder(rng, data.residues.sidechainRotationWeights[residue]));
                }
                return result;
            };

            SidechainAdjustment adjustSidechains = [&](pcg32& rng,
                                                       const assembly::Graph& graph,
                                                       const AssemblyData& data,
                                                       MutableData& mutableData,
                                                       const std::vector<size_t>& glycans)
            {
                std::vector<size_t> sidechainResiduesWithGlycanOverlap;
                sidechainResiduesWithGlycanOverlap.reserve(residueCount(graph.source));

                for (size_t residue : aliveIndices(graph.residues.nodes))
                {
                    if (data.residues.types[residue] == ResidueType::Protein &&
                        data.residues.sidechainDihedrals[residue].size() > 0 &&
                        sidechainHasGlycanOverlap(overlapSettings, graph, data, mutableData, glycans, residue))
                    {
                        sidechainResiduesWithGlycanOverlap.push_back(residue);
                    }
                }
                std::vector<size_t> residues = util::shuffleVector(rng, sidechainResiduesWithGlycanOverlap);
                std::vector<std::vector<size_t>> preferences = sidechainRotationPreferences(rng, data, residues);
                for (size_t n = 0; n < residues.size(); n++)
                {
                    setSidechainToLowestOverlapState(
                        sidechainRotamers, overlapSettings, graph, data, mutableData, preferences[n], residues[n]);
                }
            };

            SidechainAdjustment restoreSidechains = [&](pcg32& rng,
                                                        const assembly::Graph& graph,
                                                        const AssemblyData& data,
                                                        MutableData& mutableData,
                                                        const std::vector<size_t>& glycans)
            {
                std::vector<size_t> residues = util::boolsToIndices(mutableData.residueSidechainMoved);
                for (size_t residue : residues)
                {
                    restoreSidechainRotation(graph, data, mutableData, residue);
                }
                std::vector<std::vector<size_t>> preferences = sidechainRotationPreferences(rng, data, residues);
                for (size_t n = 0; n < residues.size(); n++)
                {
                    if (sidechainHasGlycanOverlap(overlapSettings, graph, data, mutableData, glycans, residues[n]))
                    {
                        setSidechainToLowestOverlapState(
                            sidechainRotamers, overlapSettings, graph, data, mutableData, preferences[n], residues[n]);
                    }
                }
            };

            SidechainAdjustment noSidechainAdjustment =
                [&](pcg32&, const assembly::Graph&, const AssemblyData&, MutableData&, const std::vector<size_t>&)
            { return; };

            auto writeOffFile = [](const assembly::Graph& graph,
                                   const AssemblyData& data,
                                   const std::vector<Coordinate>& coordinates,
                                   const std::vector<bool>& includedMolecules,
                                   const std::string& outputDir,
                                   const std::string& prefix)
            {
                assembly::Selection selection =
                    selectByAtomsAndMolecules(graph, graph.source.atomGraph.nodes.alive, includedMolecules);
                off::OffFileData offData = toOffFileData(graph, data, coordinates);
                std::string fileName = outputDir + "/" + prefix + ".off";
                util::writeToFile(
                    fileName,
                    [&](std::ostream& stream)
                    {
                        off::writeResiduesTogether(
                            stream, offData, graph, util::boolsToIndices(selection.residues), "GLYCOPROTEINBUILDER");
                    });
            };

            auto residueTER = [](const AssemblyData& data, const std::vector<std::vector<size_t>>& residueIndices)
            {
                std::vector<std::vector<bool>> result;
                result.reserve(residueIndices.size());
                for (auto& indices : residueIndices)
                {
                    result.push_back(pdb::residueTER(util::indicesToValues(data.residues.types, indices)));
                }
                return result;
            };

            auto writePdbFile = [&](const assembly::Graph& graph,
                                    const AssemblyData& data,
                                    const std::vector<Coordinate>& coordinates,
                                    const std::vector<uint>& atomNumbers,
                                    const std::vector<uint>& residueNumbers,
                                    const std::vector<bool>& includedMolecules,
                                    const std::vector<std::array<size_t, 2>>& connectionIndices,
                                    const std::string& outputDir,
                                    const std::string& prefix)
            {
                pdb::PdbFileData pdbData = toPdbFileData(
                    graph.source.indices,
                    data,
                    coordinates,
                    atomNumbers,
                    residueNumbers,
                    data.residues.chainIds,
                    headerLines);
                std::vector<std::vector<size_t>> residueIndices =
                    util::boolsToValues(graph.molecules.nodes.constituents, includedMolecules);
                std::vector<std::vector<bool>> TER = residueTER(data, residueIndices);
                util::createDirectories(outputDir);
                std::string filename = outputDir + "/" + prefix + ".pdb";
                util::writeToFile(
                    filename,
                    [&](std::ostream& stream)
                    {
                        pdb::writeAssemblyToPdb(stream, pdbData, graph, residueIndices, TER, connectionIndices);
                        pdb::theEnd(stream);
                    });
            };

            const assembly::Graph& graph = assembly.graph;
            const AssemblyData& data = assembly.data;
            const MutableData& initialState = assembly.mutableData;
            SidechainAdjustment sidechainAdjustment =
                moveOverlappingSidechains ? adjustSidechains : noSidechainAdjustment;
            SidechainAdjustment sidechainRestoration =
                moveOverlappingSidechains ? restoreSidechains : noSidechainAdjustment;

            auto isNonProteinResidue = [&graph, &data](size_t n)
            { return data.molecules.types[residueMolecule(graph.source.indices, n)] != MoleculeType::protein; };

            std::vector<std::array<size_t, 2>> atomPairsConnectingNonProteinResidues;
            for (size_t n = 0; n < edgeCount(graph.residues); n++)
            {
                const std::array<size_t, 2>& adj = graph.residues.edges.nodeAdjacencies[n];
                // Local lambda to deal make sure our substractions don't go wrong.
                auto index_diff = [](size_t a, size_t b) -> size_t { return (a > b) ? (a - b) : (b - a); };
                if (isNonProteinResidue(adj[0]) || isNonProteinResidue(adj[1]) || index_diff(adj[1], adj[0]) > 1)
                {
                    size_t atomEdgeId = sourceIndex(graph.residues.edges, n);
                    const std::array<size_t, 2> atomAdj = graph.atoms.edges.nodeAdjacencies[atomEdgeId];
                    atomPairsConnectingNonProteinResidues.push_back(atomAdj);
                }
            }

            auto nonViableSites = [&overlapSettings](
                                      const assembly::Graph& graph,
                                      const assembly::Selection& selection,
                                      const AssemblyData& data,
                                      const std::vector<double>& overlaps)
            {
                std::vector<bool> result(moleculeCount(graph.source), false);
                for (size_t molecule : includedGlycanMoleculeIds(data, selection.molecules))
                {
                    for (size_t n : assembly::moleculeSelectedAtoms(graph, selection, molecule))
                    {
                        if (overlaps[n] > overlapSettings.rejectionThreshold)
                        {
                            result[molecule] = true;
                        }
                    }
                }
                return result;
            };

            const std::vector<size_t> glycanIndices = util::indexVector(data.glycans.moleculeId);

            std::function<std::vector<bool>(const size_t&)> allLinkagesFalse = [&](size_t glycanId)
            { return std::vector<bool>(data.glycans.linkages[glycanId].size(), false); };

            std::function<std::vector<bool>(const size_t&)> allLinkagesTrue = [&](size_t glycanId)
            { return std::vector<bool>(data.glycans.linkages[glycanId].size(), true); };

            const std::vector<std::vector<bool>> dontRetainInitialPreference =
                util::vectorMap(allLinkagesFalse, glycanIndices);

            const std::vector<std::vector<bool>> defaultAllowFallback = util::vectorMap(allLinkagesTrue, glycanIndices);

            std::function<std::vector<bool>(const size_t&)> retainPreference = [&](size_t glycanId)
            {
                std::vector<bool> result(data.glycans.linkages[glycanId].size(), false);
                result[0] = resolutionSettings.forceGlycositeRotamerVariety;
                return result;
            };

            const std::vector<std::vector<bool>> retainInitialPreference =
                util::vectorMap(retainPreference, glycanIndices);

            std::function<std::vector<bool>(const size_t&)> glycanLinkageFallback = [&](size_t glycanId)
            {
                std::vector<bool> result(
                    data.glycans.linkages[glycanId].size(), resolutionSettings.allowRotamerFallback);
                result[0] = result[0] && !resolutionSettings.forceGlycositeRotamerVariety;
                return result;
            };

            const std::vector<std::vector<bool>> iterationAllowFallback =
                util::vectorMap(glycanLinkageFallback, glycanIndices);

            const std::vector<OmegaConformerSite> omegaSites = resolutionSettings.enumerateOmegaConformers
                                                                   ? alpha16OmegaSites(dihedralAngleDataTable, data)
                                                                   : std::vector<OmegaConformerSite> {};
            if (omegaSites.size() > 12)
            {
                throw std::runtime_error(
                    "enumerateOmegaConformers found " + std::to_string(omegaSites.size()) +
                    " alpha1-6 linkages. Each one doubles the number of conformations, and more than 12 is refused.");
            }

            auto runInitial = [&](pcg32 rng, StructureStats& stats)
            {
                uint persistCycles = resolutionSettings.persistCycles;
                std::function<AngleSettings(uint)> toAngleSettings = [&](uint cycle)
                {
                    double c = cycle;
                    double ratio = std::min(1.0, c / 4.0);
                    double preferenceDev = 0.5 * ratio;
                    double wiggleRoom = 3.0 - preferenceDev;
                    return AngleSettings {preferenceDev, wiggleRoom, 3.0, 2, initialMetadataOrder};
                };

                const LinkageShapeSettings shapeSettings {dontRetainInitialPreference, defaultAllowFallback};

                const std::vector<GlycanShapePreference> initialPreference = randomizeInitialShapePreference(
                    rng, toAngleSettings(0), randomizeShape, shapeSettings, data, initialState.bounds);

                GlycoproteinState state = resolveOverlapsWithWiggler(
                    rng,
                    dihedralAngleDataTable,
                    toAngleSettings,
                    sidechainAdjustment,
                    sidechainRestoration,
                    randomizeShape,
                    randomizeAllDihedrals,
                    shapeSettings,
                    initialPreference,
                    overlapSettings,
                    graph,
                    data,
                    initialState,
                    persistCycles,
                    resolutionSettings.randomizationCycles,
                    false);
                std::vector<Coordinate> resolvedCoords = getCoordinates(state.mutableData.bounds.atoms);
                assembly::Selection selection = assembly::selectAll(graph);
                std::vector<double> atomOverlaps =
                    totalOverlaps(overlapSettings, graph, data, selection, state.mutableData.bounds);
                bool serialized = resolutionSettings.prepareForMD;
                const std::string prefix = "default";
                writePdbFile(
                    graph,
                    data,
                    resolvedCoords,
                    atomNumbers(serialized, data),
                    residueNumbers(serialized, data),
                    state.mutableData.moleculeIncluded,
                    atomPairsConnectingNonProteinResidues,
                    outputDir,
                    prefix);
                if (resolutionSettings.writeOffFile)
                {
                    writeOffFile(graph, data, resolvedCoords, state.mutableData.moleculeIncluded, outputDir, prefix);
                }
                std::vector<bool> nonViable = nonViableSites(graph, selection, data, atomOverlaps);
                stats = StructureStats {
                    prefix,
                    false,
                    false,
                    state.mutableData.bounds,
                    selection,
                    nonViable,
                    state.mutableData.residueSidechainMoved,
                    atomOverlaps};
            };

            auto runIteration = [&](pcg32 rng, StructureStats& stats, const std::string& prefix)
            {
                uint persistCycles = resolutionSettings.persistCycles;
                std::function<AngleSettings(uint)> toAngleSettings = [&](uint cycle)
                {
                    double c = cycle;
                    double ratio = std::min(1.0, c / 4.0);
                    double preferenceDev = 2.0 - 1.5 * ratio;
                    double wiggleRoom = 3.0 - preferenceDev;
                    return AngleSettings {preferenceDev, wiggleRoom, 3.0, 2, randomMetadataSelection};
                };

                const LinkageShapeSettings initialShapeSettings {dontRetainInitialPreference, iterationAllowFallback};
                const std::vector<GlycanShapePreference> initialPreference = randomizeInitialShapePreference(
                    rng, toAngleSettings(0), randomizeShape, initialShapeSettings, data, initialState.bounds);

                const LinkageShapeSettings shapeSettings {retainInitialPreference, iterationAllowFallback};

                GlycoproteinState state = resolveOverlapsWithWiggler(
                    rng,
                    dihedralAngleDataTable,
                    toAngleSettings,
                    sidechainAdjustment,
                    sidechainRestoration,
                    randomizeShape,
                    randomizeAllDihedrals,
                    shapeSettings,
                    initialPreference,
                    overlapSettings,
                    graph,
                    data,
                    initialState,
                    persistCycles,
                    resolutionSettings.randomizationCycles,
                    resolutionSettings.deleteSitesUntilResolved);
                std::vector<Coordinate> coordinates = getCoordinates(state.mutableData.bounds.atoms);
                bool hasDeleted = util::contains(state.mutableData.moleculeIncluded, false);
                std::string directory = outputDir;
                assembly::Selection selection = assembly::selectByAtomsAndMolecules(
                    graph, data.atoms.includeInEachOverlapCheck, state.mutableData.moleculeIncluded);
                std::vector<double> atomOverlaps =
                    totalOverlaps(overlapSettings, graph, data, selection, state.mutableData.bounds);
                std::vector<bool> nonViable = nonViableSites(graph, selection, data, atomOverlaps);
                bool reject = util::contains(nonViable, true);
                directory = hasDeleted ? deletionDir : (reject ? rejectDir : structureDir);
                bool serialized = resolutionSettings.prepareForMD;
                writePdbFile(
                    graph,
                    data,
                    coordinates,
                    atomNumbers(serialized, data),
                    residueNumbers(serialized, data),
                    state.mutableData.moleculeIncluded,
                    atomPairsConnectingNonProteinResidues,
                    directory,
                    prefix);
                if (resolutionSettings.writeOffFile)
                {
                    writeOffFile(graph, data, coordinates, state.mutableData.moleculeIncluded, directory, prefix);
                }
                stats = StructureStats {
                    prefix,
                    reject,
                    hasDeleted,
                    state.mutableData.bounds,
                    selection,
                    nonViable,
                    state.mutableData.residueSidechainMoved,
                    atomOverlaps};
            };

            writePdbFile(
                graph,
                data,
                getCoordinates(data.atoms.initialState),
                data.atoms.numbers,
                data.residues.numbers,
                initialState.moleculeIncluded,
                atomPairsConnectingNonProteinResidues,
                outputDir,
                "unresolved");

            // order of parallel for loop is undefined, so we generate all seeds up front for determinism
            size_t sampleCount = resolutionSettings.numberOfSamples;
            size_t totalStructures = sampleCount + 1;
            if (!omegaSites.empty())
            {
                std::cout << "Each resolved structure gets gg/gt alternatives that set every alpha1-6 omega to the "
                             "named rotamer, then run the randomization and persist cycles on that omega within its "
                             "rotamer range, on the linkage immediately before it, and on linkages beyond it. A "
                             "chi1/chi2 linkage before the omega stays fixed"
                          << std::endl;
            }
            pcg32 seedingRng(resolutionSettings.rngSeed);
            std::vector<uint64_t> rngSeeds = util::randomIntegers<uint64_t>(totalStructures, seedingRng);
            std::vector<std::string> prefixes;
            prefixes.reserve(sampleCount);
            if (sampleCount == 1)
            {
                prefixes.push_back("glycoprotein");
            }
            else if (sampleCount > 1)
            {
                size_t digitCount = std::to_string(sampleCount - 1).size();
                for (size_t n = 0; n < sampleCount; n++)
                {
                    std::string digits = std::to_string(n);
                    std::string extraZeroes = std::string(digitCount - digits.size(), '0');
                    prefixes.push_back(extraZeroes + digits + "_glycoprotein");
                }
            }
            // will be initialized in loop
            stats.resize(totalStructures);

            util::setOpenMpNumberOfThreads(numThreads);
//        clang format doesn't align pragmas
/*        */#pragma omp parallel for
            for (size_t count = 0; count < totalStructures; count++)
            {
                pcg32 rng(rngSeeds[count]);
                if (count == 0)
                {
                    runInitial(rng, stats[0]);
                }
                else
                {
                    runIteration(rng, stats[count], prefixes[count - 1]);
                }
            }

            if (!omegaSites.empty())
            {
                const size_t combinationCount = size_t(1) << omegaSites.size();
                const std::string defaultRejectedDir = outputDir + "/rejected_defaults";
                auto randomizeLinkagesAfterOmega = [&](pcg32& rng,
                                                       MutableData& state,
                                                       const assembly::Selection& selection,
                                                       const std::vector<std::pair<size_t, size_t>>& linkages,
                                                       const std::vector<size_t>& changedOmegaSites,
                                                       size_t mask,
                                                       bool isSample)
                {
                    if (linkages.empty() && changedOmegaSites.empty())
                    {
                        return;
                    }
                    auto angleSettingsForCycle = [&](uint cycle)
                    {
                        double ratio = std::min(1.0, double(cycle) / 4.0);
                        if (isSample)
                        {
                            double preferenceDev = 2.0 - 1.5 * ratio;
                            double wiggleRoom = 3.0 - preferenceDev;
                            return AngleSettings {preferenceDev, wiggleRoom, 3.0, 2, randomMetadataSelection};
                        }
                        double preferenceDev = 0.5 * ratio;
                        double wiggleRoom = 3.0 - preferenceDev;
                        return AngleSettings {preferenceDev, wiggleRoom, 3.0, 2, initialMetadataOrder};
                    };
                    auto linkageSearchFor = [&](const AngleSettings& settings)
                    {
                        return AngleSearchSettings {
                            settings.searchDeviation,
                            settings.halfIntervalSearches,
                            [&, settings](
                                const DihedralAngleData& metadata,
                                double preference,
                                double deviation,
                                uint halfIntervalSearches)
                            {
                                std::pair<double, double> deviationStd = standardDeviation(settings, metadata);
                                return AngleSpacing {
                                    preference,
                                    deviation * deviationStd.first,
                                    deviation * deviationStd.second,
                                    settings.searchIncrement,
                                    halfIntervalSearches};
                            }};
                    };
                    auto stillOverlapping = [&](const MutableData& current)
                    {
                        OverlapSites sites = determineOverlapState(
                            overlapSettings.rejectionThreshold,
                            overlapSettings,
                            graph,
                            data,
                            selection,
                            current.bounds);
                        return !sites.indices.empty();
                    };
                    auto overlapSum = [&](const MutableData& current) {
                        return overlapVectorSum(totalOverlaps(overlapSettings, graph, data, selection, current.bounds));
                    };
                    auto lockNamedOmega =
                        [&](ResidueLinkageShapePreference& preference, size_t glycanId, size_t linkageIndex)
                    {
                        PermutationShapePreference* permutation = std::get_if<PermutationShapePreference>(&preference);
                        if (permutation == nullptr)
                        {
                            return;
                        }
                        for (size_t n = 0; n < omegaSites.size(); n++)
                        {
                            if (omegaSites[n].glycanId != glycanId || omegaSites[n].linkageIndex != linkageIndex)
                            {
                                continue;
                            }
                            bool useGt = (mask & (size_t(1) << n)) != 0;
                            size_t chosen = useGt ? omegaSites[n].gtMetadataIndex : omegaSites[n].ggMetadataIndex;
                            size_t bondId = omegaBondId(data, omegaSites[n]);
                            permutation->metadataOrder[omegaSites[n].bondIndex] = {chosen};
                            permutation->angles[omegaSites[n].bondIndex][chosen] =
                                dihedralAngleDataTable.entries[data.rotatableBonds.metadata[bondId][chosen]]
                                    .default_angle;
                        }
                    };
                    auto pinDownstreamOmegas = [&](MutableData& trial)
                    {
                        for (size_t n = 0; n < omegaSites.size(); n++)
                        {
                            bool downstream = false;
                            for (const std::pair<size_t, size_t>& linkage : linkages)
                            {
                                if (linkage.first == omegaSites[n].glycanId &&
                                    linkage.second == omegaSites[n].linkageIndex)
                                {
                                    downstream = true;
                                }
                            }
                            if (!downstream)
                            {
                                continue;
                            }
                            bool useGt = (mask & (size_t(1) << n)) != 0;
                            wiggleSingleOmega(
                                dihedralAngleDataTable,
                                overlapSettings,
                                graph,
                                data,
                                selection,
                                trial,
                                omegaSites[n],
                                useGt,
                                defaultSearchSettings);
                        }
                    };
                    auto wiggleAlternativeOmegas = [&](MutableData& trial, const AngleSearchSettings& linkageSearch)
                    {
                        for (size_t siteIndex : changedOmegaSites)
                        {
                            bool useGt = (mask & (size_t(1) << siteIndex)) != 0;
                            wiggleSingleOmega(
                                dihedralAngleDataTable,
                                overlapSettings,
                                graph,
                                data,
                                selection,
                                trial,
                                omegaSites[siteIndex],
                                useGt,
                                linkageSearch);
                        }
                    };
                    auto drawAndSet = [&](MutableData& trial, const AngleSettings& settings, bool allowFallback)
                    {
                        std::vector<ResidueLinkageShapePreference> preferences;
                        preferences.reserve(linkages.size());
                        for (const std::pair<size_t, size_t>& linkage : linkages)
                        {
                            size_t glycanId = linkage.first;
                            size_t linkageIndex = linkage.second;
                            size_t linkageId = data.glycans.linkages[glycanId][linkageIndex];
                            ResidueLinkageShapePreference preference = randomSingleLinkageShapePreference(
                                rng,
                                dihedralAngleDataTable,
                                settings,
                                data,
                                trial.bounds,
                                allowFallback,
                                glycanId,
                                linkageIndex,
                                randomAngle,
                                freezeGlycositeResidueConformation);
                            lockNamedOmega(preference, glycanId, linkageIndex);
                            setLinkageShapeToPreference(graph, data, trial, linkageId, preference);
                            preferences.push_back(std::move(preference));
                        }
                        return preferences;
                    };
                    auto wiggleAll = [&](MutableData& trial,
                                         const AngleSettings& settings,
                                         const assembly::Selection& wiggleSelection,
                                         size_t times,
                                         const std::vector<ResidueLinkageShapePreference>& preferences)
                    {
                        AngleSearchSettings linkageSearch = linkageSearchFor(settings);
                        for (size_t pass = 0; pass < times; pass++)
                        {
                            wiggleAlternativeOmegas(trial, linkageSearch);
                            for (size_t n = 0; n < linkages.size(); n++)
                            {
                                size_t linkageId = data.glycans.linkages[linkages[n].first][linkages[n].second];
                                wiggleLinkage(
                                    dihedralAngleDataTable,
                                    overlapSettings,
                                    graph,
                                    data,
                                    wiggleSelection,
                                    trial,
                                    linkageId,
                                    linkageSearch,
                                    preferences[n]);
                            }
                        }
                    };
                    auto persistDownstream = [&](MutableData& current)
                    {
                        uint cycle = 0;
                        uint maxCycle = 0;
                        MutableData bestState = current;
                        double bestOverlap = overlapSum(bestState);
                        bool allowFallback = isSample ? resolutionSettings.allowRotamerFallback : true;
                        const uint alternativePersistCycles = 10;
                        while (cycle < alternativePersistCycles && stillOverlapping(bestState))
                        {
                            AngleSettings settings = angleSettingsForCycle(maxCycle);
                            cycle++;
                            maxCycle = std::max(cycle, maxCycle);
                            MutableData trial = bestState;
                            std::vector<ResidueLinkageShapePreference> preferences =
                                drawAndSet(trial, settings, allowFallback);
                            assembly::Selection mainSelection = assembly::selectByAtomsAndMolecules(
                                graph, data.atoms.includeInMainOverlapCheck, trial.moleculeIncluded);
                            wiggleAll(trial, settings, mainSelection, 1, preferences);
                            wiggleAll(trial, settings, selection, 2, preferences);
                            pinDownstreamOmegas(trial);
                            double overlap = overlapSum(trial);
                            double improvement = bestOverlap - overlap;
                            if (improvement > 4.5)
                            {
                                bestState = trial;
                                bestOverlap = overlap;
                            }
                        }
                        current = bestState;
                    };
                    const size_t maxRandomizationCycles =
                        resolutionSettings.omegaCycleMultiplier * resolutionSettings.randomizationCycles;
                    AngleSettings initialAngleSettings = angleSettingsForCycle(0);
                    for (size_t randCycle = 0; randCycle < maxRandomizationCycles; randCycle++)
                    {
                        if (randCycle > 0)
                        {
                            std::vector<ResidueLinkageShapePreference> preferences =
                                drawAndSet(state, defaultGlycositeAngleSettings, true);
                            if (!stillOverlapping(state))
                            {
                                break;
                            }
                            wiggleAll(state, initialAngleSettings, selection, 2, preferences);
                            pinDownstreamOmegas(state);
                        }
                        persistDownstream(state);
                        if (!stillOverlapping(state))
                        {
                            break;
                        }
                    }
                };
                auto alternativesOf = [&](const StructureStats& parent, bool isSample, size_t parentIndex)
                {
                    std::vector<StructureStats> alternatives;
                    alternatives.reserve(combinationCount);
                    assembly::Selection selection = assembly::selectByAtomsAndMolecules(
                        graph, data.atoms.includeInEachOverlapCheck, parent.selection.molecules);
                    std::vector<std::string> parentRotamers;
                    parentRotamers.reserve(omegaSites.size());
                    for (const OmegaConformerSite& site : omegaSites)
                    {
                        parentRotamers.push_back(omegaRotamerName(currentOmegaDegrees(data, parent.bounds, site)));
                    }
                    for (size_t mask = 0; mask < combinationCount; mask++)
                    {
                        bool sameAsParent = true;
                        for (size_t n = 0; n < omegaSites.size(); n++)
                        {
                            bool useGt = (mask & (size_t(1) << n)) != 0;
                            if (parentRotamers[n] != (useGt ? "gt" : "gg"))
                            {
                                sameAsParent = false;
                            }
                        }
                        // The parent is already this gg/gt combination. A sample may instead be tg, or a mix.
                        if (sameAsParent)
                        {
                            continue;
                        }
                        MutableData state = initialState;
                        state.bounds = parent.bounds;
                        state.residueSidechainMoved = parent.residueSidechainMoved;
                        state.moleculeIncluded = parent.selection.molecules;
                        for (size_t n = 0; n < omegaSites.size(); n++)
                        {
                            bool useGt = (mask & (size_t(1) << n)) != 0;
                            if (parentRotamers[n] == (useGt ? "gt" : "gg"))
                            {
                                continue;
                            }
                            wiggleSingleOmega(
                                dihedralAngleDataTable,
                                overlapSettings,
                                graph,
                                data,
                                selection,
                                state,
                                omegaSites[n],
                                useGt,
                                defaultSearchSettings);
                        }
                        std::vector<size_t> changedOmegaSites;
                        std::vector<std::pair<size_t, size_t>> linkagesToWiggle;
                        std::vector<std::vector<bool>> alreadyIncluded(data.glycans.moleculeId.size());
                        auto includeLinkage = [&](size_t glycanId, size_t linkageIndex)
                        {
                            if (alreadyIncluded[glycanId].empty())
                            {
                                alreadyIncluded[glycanId] =
                                    std::vector<bool>(data.glycans.linkages[glycanId].size(), false);
                            }
                            if (!alreadyIncluded[glycanId][linkageIndex])
                            {
                                alreadyIncluded[glycanId][linkageIndex] = true;
                                linkagesToWiggle.push_back({glycanId, linkageIndex});
                            }
                        };
                        for (size_t n = 0; n < omegaSites.size(); n++)
                        {
                            bool useGt = (mask & (size_t(1) << n)) != 0;
                            if (parentRotamers[n] == (useGt ? "gt" : "gg"))
                            {
                                continue;
                            }
                            changedOmegaSites.push_back(n);
                            size_t glycanId = omegaSites[n].glycanId;
                            std::optional<size_t> before =
                                linkageBeforeOmega(dihedralAngleDataTable, graph, data, omegaSites[n]);
                            if (before.has_value())
                            {
                                includeLinkage(glycanId, before.value());
                            }
                            for (size_t linkageIndex : linkagesAfterOmega(graph, data, omegaSites[n]))
                            {
                                includeLinkage(glycanId, linkageIndex);
                            }
                        }
                        std::sort(linkagesToWiggle.begin(), linkagesToWiggle.end());
                        pcg32 omegaRng(resolutionSettings.rngSeed + parentIndex * 1000003ULL + mask);
                        randomizeLinkagesAfterOmega(
                            omegaRng, state, selection, linkagesToWiggle, changedOmegaSites, mask, isSample);
                        std::vector<Coordinate> coordinates = getCoordinates(state.bounds.atoms);
                        std::vector<double> atomOverlaps =
                            totalOverlaps(overlapSettings, graph, data, selection, state.bounds);
                        std::vector<bool> nonViable = nonViableSites(graph, selection, data, atomOverlaps);
                        bool reject = util::contains(nonViable, true);
                        bool hasDeleted = parent.deletions;
                        std::string prefix = omegaCombinationPrefix(omegaSites, mask, parent.filename);
                        std::string directory = isSample
                                                    ? (hasDeleted ? deletionDir : (reject ? rejectDir : structureDir))
                                                    : (reject ? defaultRejectedDir : outputDir);
                        bool serialized = resolutionSettings.prepareForMD;
                        writePdbFile(
                            graph,
                            data,
                            coordinates,
                            atomNumbers(serialized, data),
                            residueNumbers(serialized, data),
                            state.moleculeIncluded,
                            atomPairsConnectingNonProteinResidues,
                            directory,
                            prefix);
                        if (resolutionSettings.writeOffFile)
                        {
                            writeOffFile(graph, data, coordinates, state.moleculeIncluded, directory, prefix);
                        }
                        alternatives.push_back(StructureStats {
                            prefix,
                            reject,
                            hasDeleted,
                            state.bounds,
                            selection,
                            nonViable,
                            state.residueSidechainMoved,
                            atomOverlaps});
                    }
                    return alternatives;
                };

                std::vector<std::vector<StructureStats>> alternativesByParent(stats.size());
                util::setOpenMpNumberOfThreads(numThreads);
                util::createDirectories(outputDir);
                util::createDirectories(structureDir);
                util::createDirectories(rejectDir);
                util::createDirectories(deletionDir);
                util::createDirectories(defaultRejectedDir);
//        clang format doesn't align pragmas
/*        */#pragma omp parallel for
                for (size_t n = 0; n < stats.size(); n++)
                {
                    alternativesByParent[n] = alternativesOf(stats[n], n > 0, n);
                }
                std::vector<StructureStats> ordered;
                ordered.reserve(stats.size() * combinationCount);
                size_t rejectedCount = 0;
                size_t alternativeCount = 0;
                for (size_t n = 0; n < stats.size(); n++)
                {
                    ordered.push_back(stats[n]);
                    const std::vector<StructureStats>& alternatives = alternativesByParent[n];
                    for (const StructureStats& alternative : alternatives)
                    {
                        rejectedCount += alternative.rejected ? 1 : 0;
                    }
                    alternativeCount += alternatives.size();
                    ordered.insert(ordered.end(), alternatives.begin(), alternatives.end());
                }
                stats = std::move(ordered);
                std::cout << "Rejected " << rejectedCount << " of " << alternativeCount
                          << " alternative omega conformations" << std::endl;
            }
        }
    } // namespace gpbuilder
} // namespace gmml
