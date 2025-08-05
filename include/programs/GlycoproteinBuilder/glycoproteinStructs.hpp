#ifndef INCLUDE_PROGRAMS_GLYCOPROTEINBUILDER_GLYCOPROTEINSTRUCTS_HPP
#define INCLUDE_PROGRAMS_GLYCOPROTEINBUILDER_GLYCOPROTEINSTRUCTS_HPP

#include "include/assembly/assemblyTypes.hpp"
#include "include/carbohydrate/dihedralAngleSearchTypes.hpp"
#include "include/external/pcg/pcg_random.h"
#include "include/geometry/geometryTypes.hpp"
#include "include/geometry/overlap.hpp"
#include "include/graph/graphTypes.hpp"
#include "include/metadata/dihedralangledata.hpp"
#include "include/structure/atomOverlaps.hpp"

#include <functional>
#include <vector>

namespace gmml
{
    namespace gpbuilder
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
            std::vector<Element> elements;
            std::vector<double> charges;
            std::vector<Sphere> initialState;
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
            std::vector<ResidueType> types;
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
            std::vector<Sphere> sidechainPotentialBounds;
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
            std::vector<AngleWithMetadata> initialShape;
        };

        struct ResidueLinkageData
        {
            std::vector<RotamerType> rotamerTypes;
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
            DihedralAngleDataTable dihedralAngleTable;
            PotentialTable potentialTable;
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

        typedef std::function<std::vector<size_t>(pcg32&, const DihedralAngleDataTable&, const std::vector<size_t>&)>
            MetadataOrder;

        struct AngleSettings
        {
            double preferenceDeviation;
            double searchDeviation;
            double searchIncrement;
            MetadataOrder randomMetadata;
        };

        typedef std::function<void(
            const assembly::Graph&,
            const AssemblyData&,
            const assembly::Selection&,
            const AngleSettings&,
            const GlycanShapePreference&,
            MutableData&,
            size_t)>
            WiggleGlycan;

        typedef std::function<void(
            pcg32&,
            const AngleSettings& settings,
            WiggleGlycan,
            const assembly::Graph&,
            const AssemblyData&,
            MutableData&,
            const std::vector<GlycanShapePreference>&,
            const std::vector<size_t>&)>
            SidechainAdjustment;
    } // namespace gpbuilder
} // namespace gmml

#endif
