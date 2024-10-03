#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycoproteinOverlapResolution.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycoproteinStructs.hpp"
#include "includes/CentralDataStructure/Geometry/types.hpp"
#include "includes/CentralDataStructure/Geometry/boundingSphere.hpp"
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

#include <sstream>

namespace
{
    using GlycamMetadata::DihedralAngleData;
    using GlycamMetadata::DihedralAngleDataVector;
    using GlycamMetadata::RotamerType;

    cds::DihedralRotationData toRotationInputData(const AssemblyGraphs& graphs, AssemblyData& data,
                                                  OverlapWeight weight, size_t glycanId, size_t linkageId,
                                                  size_t dihedralId)
    {
        auto& linkage                          = graphs.residueLinkages[linkageId];
        auto& dihedral                         = graphs.rotatableDihedralIndices[dihedralId];
        const std::vector<size_t>& movingAtoms = graphs.rotatableDihedralIndices[dihedralId].movingAtoms;
        std::vector<bool> atomMoving(graphs.indices.atoms.size(), false);
        for (size_t atom : movingAtoms)
        {
            atomMoving[atom] = true;
        }
        std::vector<cds::Sphere>& atomBounds    = data.atoms.bounds;
        std::vector<cds::Sphere>& residueBounds = data.residues.bounds;
        std::vector<double> residueWeights      = data.residues.overlapWeights;
        cds::Sphere movingAtomBounds            = cds::boundingSphere(codeUtils::indexValues(atomBounds, movingAtoms));
        Coordinate origin                       = data.atoms.bounds[dihedral.atoms[1]].center;
        Coordinate axis                         = data.atoms.bounds[dihedral.atoms[2]].center - origin;
        Coordinate closestPointOnAxis           = origin + projection(movingAtomBounds.center - origin, axis);
        double distanceToAxis                   = length(closestPointOnAxis - movingAtomBounds.center);
        cds::Sphere movementBounds = cds::Sphere {movingAtomBounds.radius + distanceToAxis, closestPointOnAxis};

        std::vector<size_t> intersectingResidues;
        intersectingResidues.reserve(graphs.indices.residues.size());
        codeUtils::insertInto(intersectingResidues, linkage.nonReducingResidues);
        size_t glycanMolecule    = graphs.glycans[glycanId].glycanMolecule;
        size_t attachmentResidue = graphs.glycans[glycanId].attachmentResidue;
        for (size_t n = 0; n < graphs.indices.molecules.size(); n++)
        {
            if ((n != glycanMolecule) &&
                cds::spheresOverlap(constants::overlapTolerance, movementBounds, data.molecules.bounds[n]))
            {
                cds::insertIndicesOfIntersection(intersectingResidues, movementBounds, residueBounds,
                                                 graphs.molecules.nodes.elements[n]);
            }
        }

        for (size_t n : linkage.reducingResidues)
        {
            residueWeights[n] = 1.0;
        }
        for (size_t n : linkage.nonReducingResidues)
        {
            residueWeights[n] = weight.self;
        }
        return {
            atomMoving,
            atomBounds,
            residueBounds,
            residueWeights,
            graphs.residues.nodes.elements,
            {intersectingResidues, linkage.reducingResidues},
            data.residueLinkageData.overlapBonds[linkageId]
        };
    }

    std::vector<bool> groupContains(const std::vector<size_t>& group, const std::vector<size_t>& indices)
    {
        std::vector<bool> result(group.size(), false);
        for (size_t n : indices)
        {
            result[group[n]] = true;
        }
        return result;
    }

    std::array<cds::Coordinate, 4> dihedralCoordinates(const AssemblyGraphs& graphs, AssemblyData& data,
                                                       size_t dihedralId)
    {
        auto& atoms = graphs.rotatableDihedralIndices[dihedralId].atoms;
        auto coord  = [&](size_t n)
        {
            return data.atoms.bounds[atoms[n]].center;
        };
        return {coord(3), coord(2), coord(1), coord(0)};
    }

    void setDihedralAngle(const AssemblyGraphs& graphs, AssemblyData& data, size_t linkageId, size_t dihedralId,
                          const cds::AngleWithMetadata& target)
    {
        const RotatableDihedralIndices& dihedral      = graphs.rotatableDihedralIndices[dihedralId];
        std::array<cds::Coordinate, 4> dihedralCoords = dihedralCoordinates(graphs, data, dihedralId);
        auto matrix = cds::rotationTo(dihedralCoords, constants::toRadians(target.value));
        for (size_t atom : dihedral.movingAtoms)
        {
            cds::Coordinate& coord = data.atoms.bounds[atom].center;
            coord                  = matrix * coord;
        }
        size_t edgeId                           = graphs.residueLinkages[linkageId].residueEdge;
        const std::array<size_t, 2>& residueIds = graphs.residues.edges.nodeAdjacencies[edgeId];
        std::vector<bool> residueMoving         = groupContains(graphs.indices.atomResidue, dihedral.movingAtoms);
        for (size_t n = 0; n < graphs.indices.residues.size(); n++)
        {
            if (residueMoving[n] && (n != residueIds[0]) && (n != residueIds[1]))
            {
                cds::Coordinate& coord = data.residues.bounds[n].center;
                coord                  = matrix * coord;
            }
        }
        for (size_t n : residueIds)
        {
            data.residues.bounds[n] =
                cds::boundingSphere(codeUtils::indexValues(data.atoms.bounds, graphs.residues.nodes.elements[n]));
        }

        data.rotatableDihedralData.currentMetadataIndex[dihedralId] = target.metadataIndex;
    }

    std::vector<cds::AngleWithMetadata> wigglePermutationLinkage(const AssemblyGraphs& graphs, AssemblyData& data,
                                                                 size_t glycanId, size_t linkageId,
                                                                 const cds::AngleSearchSettings& settings,
                                                                 OverlapWeight weight,
                                                                 const cds::PermutationShapePreference& shapePreference)
    {
        auto& dihedrals = graphs.residueLinkages[linkageId].rotatableDihedrals;
        auto& metadata  = data.residueLinkageData.metadata[linkageId];
        std::vector<cds::AngleWithMetadata> shape;
        shape.resize(dihedrals.size());
        //  Reverse as convention is Glc1-4Gal and I want to wiggle in opposite direction i.e. from first
        //  rotatable bond in Asn outwards
        for (size_t rn = 0; rn < dihedrals.size(); rn++)
        {
            size_t n             = dihedrals.size() - 1 - rn;
            size_t dihedralId    = graphs.residueLinkages[linkageId].rotatableDihedrals[n];
            auto preference      = cds::AngleSearchPreference {settings.deviation, shapePreference.angles[n],
                                                          shapePreference.metadataOrder[n]};
            auto coordinates     = dihedralCoordinates(graphs, data, dihedralId);
            auto input           = toRotationInputData(graphs, data, weight, glycanId, linkageId, dihedralId);
            auto& metadataVector = metadata[n];
            auto best = cds::wiggleUsingRotamers(settings.angles, coordinates, codeUtils::indexVector(metadataVector),
                                                 metadataVector, preference, input);
            auto& residues = graphs.indices.residues;
            setDihedralAngle(graphs, data, linkageId, dihedralId, best.angle);
            shape[n] = best.angle;
        }
        return shape;
    }

    std::vector<cds::AngleWithMetadata> wiggleConformerLinkage(const AssemblyGraphs& graphs, AssemblyData& data,
                                                               size_t glycanId, size_t linkageId,
                                                               const cds::AngleSearchSettings& settings,
                                                               OverlapWeight weight,
                                                               const cds::ConformerShapePreference& shapePreference)
    {
        auto& dihedrals         = graphs.residueLinkages[linkageId].rotatableDihedrals;
        auto& metadata          = data.residueLinkageData.metadata[linkageId];
        size_t numberOfMetadata = shapePreference.metadataOrder.size();
        auto& preferenceAngles  = shapePreference.angles;
        auto& metadataOrder     = shapePreference.metadataOrder;
        auto& isFrozen          = shapePreference.isFrozen;
        std::vector<std::vector<cds::AngleWithMetadata>> results;
        results.resize(numberOfMetadata);
        std::vector<cds::AngleOverlap> bestOverlaps;
        bestOverlaps.resize(numberOfMetadata);
        auto index = codeUtils::indexVector(metadata[0]);
        for (size_t k = 0; k < numberOfMetadata; k++)
        {
            results[k].resize(dihedrals.size());
            cds::ConformerShapePreference pref {isFrozen, preferenceAngles, {metadataOrder[k]}};
            setLinkageShapeToPreference(graphs, data, linkageId, pref);
            //  Reverse as convention is Glc1-4Gal and I want to wiggle in opposite direction i.e. from first
            //  rotatable bond in Asn outwards
            for (size_t rn = 0; rn < dihedrals.size(); rn++)
            {
                size_t n          = dihedrals.size() - 1 - rn;
                auto& dihedral    = dihedrals[n];
                size_t dihedralId = graphs.residueLinkages[linkageId].rotatableDihedrals[n];
                auto preference   = cds::AngleSearchPreference {
                    isFrozen[n] ? 0.0 : settings.deviation, preferenceAngles[n], {metadataOrder[k]}};
                auto coordinates = dihedralCoordinates(graphs, data, dihedralId);

                auto input = toRotationInputData(graphs, data, weight, glycanId, linkageId, dihedralId);
                auto best =
                    cds::wiggleUsingRotamers(settings.angles, coordinates, index, metadata[n], preference, input);
                results[k][n]   = best.angle;
                bestOverlaps[k] = best;
                setDihedralAngle(graphs, data, linkageId, dihedralId, best.angle);
            }
        }
        size_t bestIndex = cds::bestOverlapResultIndex(bestOverlaps);
        std::vector<cds::AngleWithMetadata> shape;
        shape.reserve(dihedrals.size());
        for (size_t n = 0; n < dihedrals.size(); n++)
        {
            size_t dihedralId = graphs.residueLinkages[linkageId].rotatableDihedrals[n];
            auto& bestShape   = results[bestIndex][n];
            setDihedralAngle(graphs, data, linkageId, dihedralId, bestShape);
            shape.push_back(bestShape);
        }
        return shape;
    }

    std::vector<bool> ignoredAtomsOf(const AssemblyGraphs& graphs, size_t residueIndex, size_t atomIndex)
    {
        const std::vector<size_t>& atomAdj      = graphs.atoms.nodes.nodeAdjacencies[atomIndex];
        const std::vector<size_t>& residueAtoms = graphs.residues.nodes.elements[residueIndex];
        std::vector<bool> ignored(residueAtoms.size(), false);
        ignored[codeUtils::indexOf(residueAtoms, atomIndex)] = true;
        for (size_t adj : atomAdj)
        {
            if (graphs.indices.atomResidue[adj] == residueIndex)
            {
                ignored[codeUtils::indexOf(residueAtoms, adj)] = true;
            }
        }
        return ignored;
    };

    cds::BondedResidueOverlapInput bondedResidueOverlapInput(const AssemblyGraphs& graphs, size_t bondIndex)
    {
        auto& atomBond             = graphs.atoms.edges.nodeAdjacencies[bondIndex];
        size_t atomA               = atomBond[0];
        size_t atomB               = atomBond[1];
        size_t residueA            = graphs.indices.atomResidue[atomA];
        size_t residueB            = graphs.indices.atomResidue[atomB];
        std::vector<bool> ignoredA = ignoredAtomsOf(graphs, residueA, atomA);
        std::vector<bool> ignoredB = ignoredAtomsOf(graphs, residueB, atomB);
        return {
            {residueA, residueB},
            {ignoredA, ignoredB}
        };
    }
} // namespace

void setLinkageShape(const AssemblyGraphs& graphs, AssemblyData& data, size_t glycanId,
                     const std::vector<std::vector<cds::AngleWithMetadata>>& shape)
{
    auto& linkages = graphs.glycans[glycanId].linkages;
    for (size_t n = 0; n < linkages.size(); n++)
    {
        size_t linkageId = linkages[n];
        auto& dihedrals  = graphs.residueLinkages[linkageId].rotatableDihedrals;
        for (size_t k = 0; k < dihedrals.size(); k++)
        {
            setDihedralAngle(graphs, data, linkageId, dihedrals[k], shape[n][k]);
        }
    }
}

void setLinkageShapeToPreference(const AssemblyGraphs& graphs, AssemblyData& data, size_t linkageId,
                                 const cds::ResidueLinkageShapePreference& preference)
{
    const std::vector<size_t>& dihedrals = graphs.residueLinkages[linkageId].rotatableDihedrals;
    if (data.residueLinkageData.rotamerTypes[linkageId] == GlycamMetadata::RotamerType::conformer)
    {
        if (!std::holds_alternative<cds::ConformerShapePreference>(preference))
        {
            throw std::runtime_error("expected but did not receive ConformerShapePreference for conformer linkage");
        }
        cds::ConformerShapePreference pref = std::get<cds::ConformerShapePreference>(preference);
        size_t metadataIndex               = pref.metadataOrder[0];
        for (size_t n = 0; n < dihedrals.size(); n++)
        {
            double angle = pref.angles[n][metadataIndex];
            cds::AngleWithMetadata am {angle, angle, metadataIndex};
            setDihedralAngle(graphs, data, linkageId, dihedrals[n], am);
        }
    }
    else
    {
        if (!std::holds_alternative<cds::PermutationShapePreference>(preference))
        {
            throw std::runtime_error("expected but did not receive PermutationShapePreference for permutation linkage");
        }
        cds::PermutationShapePreference pref = std::get<cds::PermutationShapePreference>(preference);
        for (size_t n = 0; n < dihedrals.size(); n++)
        {
            size_t metadataIndex = pref.metadataOrder[n][0];
            double angle         = pref.angles[n][metadataIndex];
            cds::AngleWithMetadata am {angle, angle, metadataIndex};
            setDihedralAngle(graphs, data, linkageId, dihedrals[n], am);
        }
    }
}

cds::Overlap intraGlycanOverlaps(const AssemblyGraphs& graphs, const AssemblyData& data, size_t glycanId)
{
    auto& glycanLinkages = graphs.glycans[glycanId].linkages;
    cds::Overlap overlap = {0.0, 0.0};
    // skip first linkage as it connects to protein. We're only counting glycan atoms here
    for (size_t n = 1; n < glycanLinkages.size(); n++)
    {
        auto& linkage                        = graphs.residueLinkages[glycanLinkages[n]];
        // only take first non-reducing residue to avoid any double-counting
        size_t residueA                      = linkage.nonReducingResidues[0];
        const std::vector<size_t>& residuesB = linkage.reducingResidues;

        std::vector<cds::BondedResidueOverlapInput> bonds;
        auto& adjacencies     = graphs.residues.nodes.nodeAdjacencies[residueA];
        size_t adjacencyIndex = codeUtils::indexOf(adjacencies, residuesB[0]);
        if (adjacencyIndex < adjacencies.size())
        {
            size_t edgeIndex = graphs.residues.nodes.edgeAdjacencies[residueA][adjacencyIndex];
            size_t bondIndex = graphs.residues.edges.indices[edgeIndex];
            bonds.push_back(bondedResidueOverlapInput(graphs, bondIndex));
        }
        overlap += cds::CountOverlappingAtoms(
            cds::ResidueAtomOverlapInputReference {data.atoms.bounds, data.residues.bounds,
                                                   graphs.residues.nodes.elements, data.residues.overlapWeights},
            bonds, {residueA}, residuesB);
    }
    return overlap;
}

cds::Overlap moleculeOverlaps(const AssemblyGraphs& graphs, const AssemblyData& data, size_t moleculeA,
                              size_t moleculeB)
{
    const cds::Sphere& boundsA = data.molecules.bounds[moleculeA];
    const cds::Sphere& boundsB = data.molecules.bounds[moleculeB];
    if (!cds::spheresOverlap(constants::overlapTolerance, boundsA, boundsB))
    {
        return cds::Overlap {0.0, 0.0};
    }
    else
    {
        std::vector<cds::BondedResidueOverlapInput> bonds;
        auto& adjacencies     = graphs.molecules.nodes.nodeAdjacencies[moleculeA];
        size_t adjacencyIndex = codeUtils::indexOf(adjacencies, moleculeB);
        if (adjacencyIndex < adjacencies.size())
        {
            size_t edgeIndex        = graphs.molecules.nodes.edgeAdjacencies[moleculeA][adjacencyIndex];
            size_t residueBondIndex = graphs.molecules.edges.indices[edgeIndex];
            size_t atomBondIndex    = graphs.residues.edges.indices[residueBondIndex];
            bonds.push_back(bondedResidueOverlapInput(graphs, atomBondIndex));
        }
        std::vector<size_t> residuesA;
        std::vector<size_t> residuesB;
        cds::insertIndicesOfIntersection(residuesA, boundsB, data.residues.bounds,
                                         graphs.molecules.nodes.elements[moleculeA]);
        cds::insertIndicesOfIntersection(residuesB, boundsA, data.residues.bounds,
                                         graphs.molecules.nodes.elements[moleculeB]);
        return cds::CountOverlappingAtoms(
            cds::ResidueAtomOverlapInputReference {data.atoms.bounds, data.residues.bounds,
                                                   graphs.residues.nodes.elements, data.residues.overlapWeights},
            bonds, residuesA, residuesB);
    }
}

cds::Overlap moleculeResidueOverlaps(const AssemblyGraphs& graphs, const AssemblyData& data, size_t molecule,
                                     size_t residue)
{
    const cds::Sphere& moleculeBounds = data.molecules.bounds[molecule];
    const cds::Sphere& residueBounds  = data.residues.bounds[residue];
    size_t residueMolecule            = graphs.indices.residueMolecule[residue];
    if (!cds::spheresOverlap(constants::overlapTolerance, moleculeBounds, residueBounds))
    {
        return cds::Overlap {0.0, 0.0};
    }
    else
    {
        std::vector<cds::BondedResidueOverlapInput> bonds;
        auto& adjacencies     = graphs.molecules.nodes.nodeAdjacencies[molecule];
        size_t adjacencyIndex = codeUtils::indexOf(adjacencies, residueMolecule);
        if (adjacencyIndex < adjacencies.size())
        {
            size_t edgeIndex        = graphs.molecules.nodes.edgeAdjacencies[molecule][adjacencyIndex];
            size_t residueBondIndex = graphs.molecules.edges.indices[edgeIndex];
            size_t atomBondIndex    = graphs.residues.edges.indices[residueBondIndex];
            bonds.push_back(bondedResidueOverlapInput(graphs, atomBondIndex));
        }
        return cds::CountOverlappingAtoms(
            cds::ResidueAtomOverlapInputReference {data.atoms.bounds, data.residues.bounds,
                                                   graphs.residues.nodes.elements, data.residues.overlapWeights},
            bonds, {residue}, graphs.molecules.nodes.elements[molecule]);
    }
}

cds::Overlap totalOverlaps(OverlapWeight weight, const AssemblyGraphs& graphs, const AssemblyData& data)
{
    cds::Overlap overlap {0.0, 0.0};
    const std::vector<GlycanIndices>& glycosites = graphs.glycans;
    for (size_t n : graphs.proteinMolecules)
    {
        for (auto& glycan : glycosites)
        {
            overlap += moleculeOverlaps(graphs, data, n, glycan.glycanMolecule);
        }
    }

    for (size_t n = 0; n < glycosites.size(); n++)
    {
        overlap += intraGlycanOverlaps(graphs, data, n) * weight.self;
        for (size_t k = n + 1; k < glycosites.size(); k++)
        {
            overlap += moleculeOverlaps(graphs, data, glycosites[n].glycanMolecule, glycosites[k].glycanMolecule);
        }
    }

    return overlap;
}

std::vector<size_t> determineSitesWithOverlap(const std::vector<size_t>& movedSites, const AssemblyGraphs& graphs,
                                              const AssemblyData& data)
{
    const std::vector<GlycanIndices> glycans = graphs.glycans;
    auto hasProteinOverlap                   = [&](size_t n)
    {
        cds::Overlap overlap {0.0, 0.0};
        for (size_t k : graphs.proteinMolecules)
        {
            overlap += moleculeOverlaps(graphs, data, k, glycans[n].glycanMolecule);
        }
        return overlap.count > 0;
    };
    auto hasSelfOverlap = [&](size_t n)
    {
        auto overlap = intraGlycanOverlaps(graphs, data, n);
        return overlap.count > 0;
    };
    std::vector<bool> justMoved(glycans.size(), false);
    for (size_t n : movedSites)
    {
        justMoved[n] = true;
    }
    std::vector<bool> glycanOverlap(glycans.size(), false);
    for (size_t n : movedSites)
    {
        for (size_t k = n + 1; k < glycans.size(); k++)
        {
            if (!(glycanOverlap[n] && glycanOverlap[k]))
            {
                if (moleculeOverlaps(graphs, data, glycans[n].glycanMolecule, glycans[k].glycanMolecule).count > 0 ||
                    moleculeResidueOverlaps(graphs, data, glycans[k].glycanMolecule, glycans[n].attachmentResidue)
                            .count > 0)
                {
                    glycanOverlap[n] = true;
                    glycanOverlap[k] = true;
                }
            }
        }
    }
    std::vector<size_t> indices;
    for (size_t n = 0; n < glycans.size(); n++)
    {
        // glycans which haven't moved won't overlap with protein or themselves (at least not more than before)
        if (glycanOverlap[n] || (justMoved[n] && (hasProteinOverlap(n) || hasSelfOverlap(n))))
        {
            indices.push_back(n);
        }
    }
    return indices;
}

std::vector<cds::AngleWithMetadata> wiggleLinkage(const AssemblyGraphs& graphs, AssemblyData& data, size_t glycanId,
                                                  size_t linkageId, const cds::AngleSearchSettings& searchSettings,
                                                  OverlapWeight weight,
                                                  const cds::ResidueLinkageShapePreference& shapePreference)
{
    switch (data.residueLinkageData.rotamerTypes[linkageId])
    {
        case RotamerType::permutation:
            {
                auto preference = std::get<cds::PermutationShapePreference>(shapePreference);
                return wigglePermutationLinkage(graphs, data, glycanId, linkageId, searchSettings, weight, preference);
            }
        case RotamerType::conformer:
            {
                auto preference = std::get<cds::ConformerShapePreference>(shapePreference);
                return wiggleConformerLinkage(graphs, data, glycanId, linkageId, searchSettings, weight, preference);
            }
    }
    throw std::runtime_error("unhandled linkage shape preference in glycoproteinOverlapResolution wiggleLinkage");
}

std::vector<std::vector<cds::AngleWithMetadata>>
wiggleGlycan(const AssemblyGraphs& graphs, AssemblyData& data, size_t glycanId,
             const cds::AngleSearchSettings& searchSettings, OverlapWeight weight,
             const std::vector<cds::ResidueLinkageShapePreference>& preferences)
{
    auto toWeight = [](const std::vector<cds::Residue*>& residues, double a)
    {
        return std::vector<double>(residues.size(), a);
    };
    const std::vector<size_t>& linkages = graphs.glycans[glycanId].linkages;
    std::vector<std::vector<cds::AngleWithMetadata>> shape;
    shape.resize(linkages.size());
    // wiggling twice gives the first linkages a second chance to resolve in a better structure
    for (size_t k = 0; k < 2; k++)
    {
        for (size_t n = 0; n < linkages.size(); n++)
        {
            size_t linkageId = graphs.glycans[glycanId].linkages[n];
            shape[n]         = wiggleLinkage(graphs, data, glycanId, linkageId, searchSettings, weight, preferences[n]);
        }
    }
    updateGlycanBounds(graphs, data, glycanId);
    return shape;
}

void updateResidueBounds(const AssemblyGraphs& graphs, AssemblyData& data, size_t index)
{
    data.residues.bounds[index] =
        cds::boundingSphere(codeUtils::indexValues(data.atoms.bounds, graphs.residues.nodes.elements[index]));
}

void updateResidueMoleculeBounds(const AssemblyGraphs& graphs, AssemblyData& data, size_t index)
{
    size_t moleculeIndex = graphs.indices.residueMolecule[index];
    data.molecules.bounds[moleculeIndex] =
        cds::boundingSphereIncluding(data.molecules.bounds[moleculeIndex], data.residues.bounds[index]);
}

void updateMoleculeBounds(const AssemblyGraphs& graphs, AssemblyData& data, size_t index)
{
    data.molecules.bounds[index] =
        cds::boundingSphere(codeUtils::indexValues(data.residues.bounds, graphs.molecules.nodes.elements[index]));
}

void updateGlycanBounds(const AssemblyGraphs& graphs, AssemblyData& data, size_t glycanId)
{
    const GlycanIndices& glycan = graphs.glycans[glycanId];
    size_t siteResidue          = glycan.attachmentResidue;
    updateResidueBounds(graphs, data, siteResidue);
    updateResidueMoleculeBounds(graphs, data, siteResidue);
    size_t moleculeId = glycan.glycanMolecule;
    for (size_t residue : graphs.molecules.nodes.elements[moleculeId])
    {
        updateResidueBounds(graphs, data, residue);
    }
    updateMoleculeBounds(graphs, data, moleculeId);
}

GlycoproteinState randomDescent(pcg32 rng, LinkageShapeRandomizer randomizeShape, WiggleGlycan wiggleGlycan,
                                int persistCycles, OverlapWeight overlapWeight, const AssemblyGraphs& graphs,
                                AssemblyData& data, std::vector<std::vector<cds::ResidueLinkage>>& glycosidicLinkages,
                                const GlycoproteinState& initialState)
{
    std::stringstream logss;
    logss << "Random Decent, persisting for " << persistCycles << " cycles.\n";
    gmml::log(__LINE__, __FILE__, gmml::INF, logss.str());

    auto globalOverlap        = initialState.overlap;
    auto overlapSites         = initialState.overlapSites;
    auto glycositePreferences = initialState.preferences;
    auto glycositeShape       = initialState.shape;
    int cycle                 = 0;
    while ((cycle < persistCycles) && (!overlapSites.empty()))
    {
        gmml::log(__LINE__, __FILE__, gmml::INF,
                  "Cycle " + std::to_string(cycle) + "/" + std::to_string(persistCycles));
        ++cycle;
        cds::Overlap newGlobalOverlap = globalOverlap;
        for (auto& glycanId : codeUtils::shuffleVector(rng, overlapSites))
        {
            auto& linkages                        = glycosidicLinkages[glycanId];
            const std::vector<size_t>& linkageIds = graphs.glycans[glycanId].linkages;
            auto localOverlaps                    = [&]()
            {
                const std::vector<GlycanIndices>& glycans = graphs.glycans;
                auto& thisGlycan                          = glycans[glycanId];
                cds::Overlap overlap = intraGlycanOverlaps(graphs, data, glycanId) * overlapWeight.self;
                for (size_t n : graphs.proteinMolecules)
                {
                    overlap += moleculeOverlaps(graphs, data, n, thisGlycan.glycanMolecule);
                }

                for (size_t n = 0; n < glycans.size(); n++)
                {
                    if (n != glycanId)
                    {
                        size_t other = glycans[n].glycanMolecule;
                        overlap      += moleculeOverlaps(graphs, data, thisGlycan.glycanMolecule, other);
                        overlap      += moleculeResidueOverlaps(graphs, data, other, thisGlycan.attachmentResidue);
                    }
                }
                return overlap;
            };
            auto previousOverlap = localOverlaps();
            auto preferences     = randomizeShape(linkages);
            auto recordedShape   = glycositeShape[glycanId];
            for (size_t n = 0; n < linkageIds.size(); n++)
            {
                setLinkageShapeToPreference(graphs, data, linkageIds[n], preferences[n]);
            }
            auto currentShape = wiggleGlycan(graphs, data, glycanId, overlapWeight, preferences);
            updateGlycanBounds(graphs, data, glycanId);
            auto newOverlap = localOverlaps();
            auto diff       = newOverlap + (previousOverlap * -1);
            bool isWorse    = cds::compareOverlaps(newOverlap, previousOverlap) > 0;
            if (isWorse)
            {
                setLinkageShape(graphs, data, glycanId, recordedShape);
                updateGlycanBounds(graphs, data, glycanId);
            }
            else
            {
                newGlobalOverlap               += diff;
                glycositePreferences[glycanId] = preferences;
                glycositeShape[glycanId]       = currentShape;
                gmml::log(__LINE__, __FILE__, gmml::INF,
                          "RandomDescent accepted a change of " + std::to_string(diff.count));
            }
        }
        if (cds::compareOverlaps(globalOverlap, newGlobalOverlap) > 0)
        {
            cycle = 0;
        }
        globalOverlap = newGlobalOverlap;
        overlapSites  = determineSitesWithOverlap(overlapSites, graphs, data);
    }
    return {globalOverlap, overlapSites, glycositePreferences, glycositeShape};
}
