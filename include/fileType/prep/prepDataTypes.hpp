#ifndef INCLUDE_FILETYPE_PREP_PREPDATATYPES_HPP
#define INCLUDE_FILETYPE_PREP_PREPDATATYPES_HPP

#include "include/geometry/geometryTypes.hpp"
#include "include/graph/graphTypes.hpp"

#include <string>
#include <vector>

namespace gmml
{
    namespace prep
    {
        enum TopologicalType
        {
            kTopTypeE,
            kTopTypeS,
            kTopTypeB,
            kTopType3,
            kTopType4,
            kTopTypeM
        };

        static const std::vector<std::string> topologicalTypeNames {"E", "S", "B", "3", "4", "M"};

        enum CoordinateType
        {
            kINT,
            kXYZ
        };

        static const std::vector<std::string> coordinateTypeNames {"INT", "XYZ"};

        enum OutputFormat
        {
            kFormatted = 0,
            kBinary = 1
        };

        static const std::vector<std::string> outputFormatNames {"NBinary", "Binary"};

        enum GeometryType
        {
            kGeometryCorrect,
            kGeometryChange
        };

        static const std::vector<std::string> geometryTypesNames {"CORRECT", "CHANGE"};

        enum DummyAtomPosition
        {
            kPositionAll,
            kPositionBeg
        };

        static const std::vector<std::string> dummyAtomPositionNames {"ALL", "BEG"};

        enum DummyAtomOmission
        {
            kOmit,
            kNomit
        };

        static const std::vector<std::string> dummyAtomOmissionNames {"OMIT", "NOMIT"};
        static const std::vector<std::string> dummyAtomOmissionYesNo {"YES", "NO"};

        enum SectionType
        {
            kSectionLoop,
            kSectionImproper,
            kSectionDone,
            kSectionBlank,
            kSectionOther
        };

        struct AtomData
        {
            /*!< Sample line of the atom section of a prep file: 4 H1   H1  M  3  2  1  1.000    90.0     180.0     0.0
             */
            std::vector<uint> number;
            std::vector<std::string> name;
            std::vector<std::string> type;
            std::vector<TopologicalType> topologicalType;
            std::vector<uint> bondIndex;
            std::vector<uint> angleIndex;
            std::vector<uint> dihedralIndex;
            std::vector<double> bondLength;
            std::vector<double> angle;
            std::vector<double> dihedral;
            std::vector<double> charge;
            std::vector<Coordinate> coordinate;
        };

        typedef std::vector<std::string> Dihedral; // This looks poorly named
        typedef std::vector<Dihedral> DihedralVector;

        struct ResidueData
        {
            std::vector<std::string> name;
            std::vector<std::string> title;
            std::vector<CoordinateType> coordinateType;
            std::vector<OutputFormat> outputFormat;
            std::vector<GeometryType> geometryType;
            std::vector<DummyAtomOmission> dummyAtomOmission;
            std::vector<std::string> dummyAtomType;
            std::vector<DummyAtomPosition> dummyAtomPosition;
            std::vector<double> charge;
            std::vector<DihedralVector> improperDihedrals;
            std::vector<std::vector<std::pair<std::string, std::string>>> loops;
            //!< Loops; fill by all lines between LOOP title in each residue section of the file and a blank
            //!< line in that section
            //	!< End of each residue section gets marked by DONE
            //	! \example
            //	 * A Sample of residue section in a prep file:
            //
            //            ROH for aglycon
            //
            //            ROH    INT 0
            //            CORRECT OMIT DU BEG
            //            -0.194
            //             1 DUMM DU  M  0 -1 -2  0.000     0.0       0.0     0.0
            //             2 DUMM DU  M  1  0 -1  1.000     0.0       0.0     0.0
            //             3 DUMM DU  M  2  1  0  1.000    90.0       0.0     0.0
            //             4 HO1  HO  M  3  2  1  1.000    90.0     180.0     0.445
            //             5 O1   OH  M  4  3  2  0.960   107.0     180.0    -0.639
            //
            //            DONE
        };

        struct PrepData
        {
            AtomData atoms;
            ResidueData residues;
            graph::Database atomGraph;
            std::vector<size_t> atomResidue;
            size_t atomCount = 0;
            size_t residueCount = 0;
        };
    } // namespace prep
} // namespace gmml

#endif
