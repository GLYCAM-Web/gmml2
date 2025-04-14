#ifndef INCLUDES_CENTRALDATASTRUCTURE_READERS_PREP_PREPRESIDUE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_READERS_PREP_PREPRESIDUE_HPP

#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/Readers/Prep/prepAtom.hpp"
#include <string>
#include <map>
#include <vector>
#include <istream>
#include <ostream>

namespace prep
{
    enum CoordinateType
    {
        kINT,
        kXYZ
    };

    enum OutputFormat
    {
        kFormatted = 0,
        kBinary    = 1
    };

    enum GeometryType
    {
        kGeometryCorrect,
        kGeometryChange
    };

    enum DummyAtomPosition
    {
        kPositionAll,
        kPositionBeg
    };

    enum DummyAtomOmission
    {
        kOmit,
        kNomit
    };

    enum SectionType
    {
        kSectionLoop,
        kSectionImproper,
        kSectionDone,
        kSectionBlank,
        kSectionOther
    };

    typedef std::vector<std::string> Dihedral; // This looks poorly named
    typedef std::vector<Dihedral> DihedralVector;

    struct PrepResidueProperties
    {

        std::string title = ""; //!< Residue title; fill by the first line of each residue section of the file
        CoordinateType coordinateType = prep::kINT;       //!< Coordinate type(INT, XYZ); fill by the 2nd column of the
                                                          //!< third line of each residue section of the file
        OutputFormat outputFormat     = prep::kFormatted; //!< Output format(Binary=1,Formatted=1); fill by the third
                                                      //!< column of the 3rd line of each residue section of the file
        GeometryType geometryType =
            prep::kGeometryCorrect; //!< Geometry type(CORRECT, CHANGE); fill by the first column of the 4th line of
                                    //!< each residue section of the file
        DummyAtomOmission dummyAtomOmission = prep::kOmit; //!< Dummy atom omission(OMIT, NOMIT); fill by the 3rd column
                                                           //!< of the 4th line of each residue section of the file
        std::string dummyAtomType =
            "DU"; //!< Dummy atom type; fill by the 4th column of the 4th line of each residue section of the file
        DummyAtomPosition dummyAtomPosition =
            prep::kPositionBeg; //!< Dummy atom position(ALL, BEG); fill by the 5th column of the 4th line of each
                                //!< residue section of the file
        double charge = 0.0; //!< Total charge of the residue; fill by the 5th line of each residue section of the file
        DihedralVector improperDihedrals; //!< Improper dihedrals; fill by all lines between IMPROPER title in each
                                          //!< residue section of the file and a blank line in that section
        std::vector<std::pair<std::string, std::string>>
            loops; //!< Loops; fill by all lines between LOOP title in each residue section of the file and a blank
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

    class PrepResidue : public cds::Residue
    {
      public:
        PrepResidue(std::istream& in_file, std::string& line);
        PrepResidue()  = default;
        ~PrepResidue() = default;
        PrepResidue(PrepResidue&& other) noexcept; // Move Ctor
        PrepResidue(const PrepResidue& other);     // Copy Ctor
        PrepResidue& operator=(PrepResidue other); // Move and Copy assignment operator

        friend void swap(PrepResidue& lhs, PrepResidue& rhs)
        {
            using std::swap;
            swap(static_cast<cds::Residue&>(lhs), static_cast<cds::Residue&>(rhs));
            swap(lhs.properties, rhs.properties);
        }

        void Generate3dStructure();
        void DeleteDummyAtoms();
        void SetConnectivities();
        std::vector<std::string> GetAtomNames() const;
        std::vector<std::string> GetHeavyAtomNames() const;
        double CalculatePrepResidueCharge();
        std::string toString() const;
        void Write(std::ostream& stream);

        void ExtractResidueName(std::istream& ss);
        void ExtractResidueCoordinateType(std::istream& ss);
        void ExtractResidueOutputFormat(std::istream& ss);
        void ExtractResidueGeometryType(std::istream& ss);
        void ExtractResidueDummyAtomOmission(std::istream& ss);
        void ExtractResidueDummyAtomPosition(std::istream& ss);
        prep::SectionType ExtractSectionType(std::string& line);
        void ExtractLoops(std::istream& in_file);
        void ExtractImproperDihedral(std::istream& in_file);

        std::string GetStringFormatOfCoordinateType(CoordinateType coordinate_type) const;
        std::string GetStringFormatOfGeometryType(GeometryType geometry_type) const;
        std::string GetStringFormatOfDummyAtomPosition(DummyAtomPosition dummy_atom_position) const;
        std::string GetStringFormatOfDummyAtomOmission(DummyAtomOmission dummy_atom_omission) const;

        PrepResidueProperties properties;
    };
} // namespace prep
#endif
