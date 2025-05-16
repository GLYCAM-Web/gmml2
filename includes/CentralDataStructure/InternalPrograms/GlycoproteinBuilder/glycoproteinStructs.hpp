#ifndef INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GLYCOPROTEINSTRUCTS_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GLYCOPROTEINSTRUCTS_HPP

#include "includes/Graph/graphTypes.hpp"
#include "includes/Assembly/assemblyGraph.hpp"
#include "includes/Assembly/assemblyBounds.hpp"
#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
#include "includes/CentralDataStructure/Geometry/overlap.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralAngleSearch.hpp"
#include "includes/CentralDataStructure/Overlaps/atomOverlaps.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/MolecularMetadata/GLYCAM/dihedralangledata.hpp"
#include "includes/External_Libraries/PCG/pcg_random.h"

#include <functional>
#include <vector>

namespace glycoproteinBuilder
{
    enum class MoleculeType
    {
        protein,
        glycan
    };

    struct AtomData
    {
        std::vector<std::string> names;
        std::vector<std::string> types;
        std::vector<uint> numbers;
        std::vector<uint> serializedNumbers;
        std::vector<uint> atomicNumbers;
        std::vector<std::string> elementStrings;
        std::vector<MolecularMetadata::Element> elements;
        std::vector<double> charges;
        std::vector<cds::Sphere> initialState;
        std::vector<bool> all;
        std::vector<bool> includeInEachOverlapCheck;
        std::vector<bool> includeInMainOverlapCheck;
        std::vector<bool> partOfMovableSidechain;
    };

    struct SidechainDihedral
    {
        std::array<size_t, 4> atoms;
        std::vector<size_t> movingAtoms;
    };

    struct ResidueData
    {
        std::vector<std::string> names;
        std::vector<cds::ResidueType> types;
        std::vector<bool> hasAllExpectedAtoms;
        std::vector<double> phiAngles;
        std::vector<double> psiAngles;
        std::vector<std::string> ids;
        std::vector<uint> numbers;
        std::vector<uint> serializedNumbers;
        std::vector<std::string> chainIds;
        std::vector<std::vector<SidechainDihedral>> sidechainDihedrals;
        std::vector<std::vector<size_t>> sidechainRotations;
        std::vector<std::vector<double>> sidechainRotationWeights;
        std::vector<cds::Sphere> sidechainPotentialBounds;
    };

    struct ResidueEdgeData
    {
        std::vector<std::array<std::vector<bool>, 2>> atomsCloseToEdge;
    };

    struct MoleculeData
    {
        std::vector<MoleculeType> types;
    };

    struct GlycanData
    {
        std::vector<size_t> attachmentResidue;
        std::vector<size_t> moleculeId;
        std::vector<std::vector<size_t>> linkages;
    };

    struct RotatableDihedralData
    {
        std::vector<std::vector<size_t>> metadata;
        std::vector<cds::AngleWithMetadata> initialShape;
    };

    struct ResidueLinkageData
    {
        std::vector<GlycamMetadata::RotamerType> rotamerTypes;
        std::vector<bool> isGlycositeLinkage;
    };

    struct RotatableDihedralIndices
    {
        std::array<size_t, 4> atoms;
        std::vector<size_t> movingAtoms;
    };

    struct ResidueLinkageIndices
    {
        size_t residueEdge;
        std::vector<size_t> rotatableDihedrals;
    };

    struct AssemblyIndices
    {
        std::vector<size_t> proteinMolecules;
        std::vector<RotatableDihedralIndices> rotatableDihedrals;
        std::vector<ResidueLinkageIndices> residueLinkages;
    };

    struct AssemblyData
    {
        AtomData atoms;
        ResidueData residues;
        MoleculeData molecules;
        ResidueEdgeData residueEdges;
        GlycanData glycans;
        RotatableDihedralData rotatableDihedralData;
        ResidueLinkageData residueLinkageData;
        AssemblyIndices indices;
        GlycamMetadata::DihedralAngleDataTable dihedralAngleTable;
        MolecularMetadata::PotentialTable potentialTable;
        double overlapTolerance;
        double overlapRejectionThreshold;
    };

    struct MutableData
    {
        assembly::Bounds bounds;
        std::vector<size_t> dihedralCurrentMetadata;
        std::vector<bool> moleculeIncluded;
        std::vector<bool> residueSidechainMoved;
    };

    struct GlycoproteinAssembly
    {
        assembly::Graph graph;
        AssemblyData data;
        MutableData mutableData;
    };

    inline std::vector<bool> glycanIncluded(const AssemblyData& data, const std::vector<bool>& includedMolecules)
    {
        return codeUtils::indicesToValues(includedMolecules, data.glycans.moleculeId);
    }

    inline std::vector<size_t> includedGlycanMoleculeIds(const AssemblyData& data,
                                                         const std::vector<bool>& includedMolecules)
    {
        return codeUtils::indicesToValues(data.glycans.moleculeId,
                                          codeUtils::boolsToIndices(glycanIncluded(data, includedMolecules)));
    }

    inline std::vector<size_t> includedGlycanIndices(const AssemblyData& data,
                                                     const std::vector<bool>& includedMolecules)
    {
        return codeUtils::boolsToIndices(glycanIncluded(data, includedMolecules));
    }

    inline const std::vector<uint>& atomNumbers(bool serialized, const AssemblyData& data)
    {
        return serialized ? data.atoms.serializedNumbers : data.atoms.numbers;
    }

    inline const std::vector<uint>& residueNumbers(bool serialized, const AssemblyData& data)
    {
        return serialized ? data.residues.serializedNumbers : data.residues.numbers;
    }

    typedef std::function<std::vector<size_t>(pcg32&, const GlycamMetadata::DihedralAngleDataTable&,
                                              const std::vector<size_t>&)>
        MetadataOrder;

    struct AngleSettings
    {
        double preferenceDeviation;
        double searchDeviation;
        double searchIncrement;
        MetadataOrder randomMetadata;
    };

    typedef std::function<void(const assembly::Graph&, const AssemblyData&, const assembly::Selection&,
                               const AngleSettings&, const cds::GlycanShapePreference&, MutableData&, size_t)>
        WiggleGlycan;

    typedef std::function<void(pcg32&, const AngleSettings& settings, WiggleGlycan, const assembly::Graph&,
                               const AssemblyData&, MutableData&, const std::vector<cds::GlycanShapePreference>&,
                               const std::vector<size_t>&)>
        SidechainAdjustment;

} // namespace glycoproteinBuilder
#endif
