#include <iostream>
#include <regex>
#include <string>
#include "includes/MolecularMetadata/GLYCAM/dihedralangledata.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/containers.hpp"

namespace
{
    using GlycamMetadata::AngleLimit;
    using GlycamMetadata::AngleStd;
    using GlycamMetadata::DihedralAngleData;
    using GlycamMetadata::RotamerType;

    // Struct is copied here for reference.
    // struct DihedralAngleData
    //{
    //    std::string linking_atom1_ ;
    //    std::string linking_atom2_ ;
    //    std::string dihedral_angle_name_ ;
    //    double default_angle;
    //    std::variant<AngleLimit, AngleStd> angle_deviation;
    //    double weight_;
    //    std::string rotamer_type_ ; // permutation or conformer
    //    std::string rotamer_name_ ;
    //    int number_of_bonds_from_anomeric_carbon_;
    //    int index_ ; // Used to indicate whether multiple entries are meant to overwrite each other or generate an
    //    additional angle StringVector residue1_conditions_ ; StringVector residue2_conditions_ ; std::string atom1_ ;
    //    std::string atom2_ ;
    //    std::string atom3_ ;
    //    std::string atom4_ ;
    //} ;

    /*
     * A note on index_.
     * number_of_bonds_from_anomeric_carbon_. So Phi is 1, Psi is 2, Omg is 3.
     * The index refers the rotamer number. If there are two Phi rotamers, they will have an index of 1 and 2.
     * Chi angle index numbering varies depending on side chain length, so in ASN the Chi1 is 4 bonds away from the
     * sugar, so would be 4 If two entries have matching regex and have the same index_ number (e.g. 1), the first will
     * be overwritten. If two entries have matching regex and different index_ numbers (e.g. 1,2,3) they will all be
     * used to create multiple rotamers/conformers If two entries have different regex, but apply to the same dihedral
     * angle (e.g. Phi), give them the same index_ number (e.g. 1).
     */

    // clang-format off
    std::vector<DihedralAngleData> dihedralAngleDataVector_ {
      { // Regex1  , Regex2   , Name   , Angle  , Deviation               , Weight, Entry Type               , Name , B , I , Res1 Condition , Res2 Conditions           , Atom names                                                               // Atom names this applies to
          { "C1"   , "O[1-9]" , "Phi"  , 180.0  , AngleLimit{25.0, 25.0}  , 1.0   , RotamerType::permutation , "g"  , 1 , 1 , {"aldose"}     , {"monosaccharide"}                  , "C2" , "C1" , "O." , "C."  }, // Phi should be C2-C1(ano)-Ox-Cx, or C1-C2(ano)-Ox-Cx
          { "C2"   , "O[1-9]" , "Phi"  , -60.0  , AngleLimit{25.0, 25.0}  , 1.0   , RotamerType::permutation , "g"  , 1 , 1 , {"n-carbon=6", "ketose", "alpha"}     , {"monosaccharide"}            , "C1" , "C2" , "O." , "C."  }, // Phi is defined by C1-C2(ano)-Ox-Cx for ketoses like Fru
          { "C2"   , "O[1-9]" , "Phi"  ,  60.0  , AngleLimit{25.0, 25.0}  , 1.0   , RotamerType::permutation , "-g" , 1 , 1 , {"n-carbon=6", "ketose", "beta"}     , {"monosaccharide"}            , "C1" , "C2" , "O." , "C."  }, // Phi is defined by C1-C2(ano)-Ox-Cx for ketoses like Fru
          { "C2"   , "O[1-9]" , "Phi"  , 180.0  , AngleLimit{25.0, 25.0}  , 1.0   , RotamerType::permutation , "t"  , 1 , 1 , {"ketose", "ulosonate"}     , {"monosaccharide"}      , "C1" , "C2" , "O." , "C."  }, // Phi should be C2-C1(ano)-Ox-Cx, or C1-C2(ano)-Ox-Cx
          // Sialic acid type 2- linkages:
          { "C2"   , "O[3-6]" , "Phi"  , -60.0  , AngleLimit{25.0, 25.0}  , 1.0   , RotamerType::permutation , "-g" , 1 , 2 , {"ulosonate" }  , {"monosaccharide"}         , "C1" , "C2" , "O." , "C."  },
          // Generic psi linkages, why is this labelled "ap"?
          { "C."   , "O[1-5]" , "Psi"  ,   0.0  , AngleLimit{40.0, 40.0}  , 1.0   , RotamerType::permutation , "ap" , 2 , 1 , {"monosaccharide"}       , {"monosaccharide"}                  , "C." , "O." , "C." , "H."  }, // Psi should be C(ano)-Ox-Cx-Hx, if Cx is ring, otherwise, C(ano)-Ox-Cx-C(x-1)
          { "C."   , "O[6-9]" , "Psi"  , 180.0  , AngleLimit{40.0, 40.0}  , 1.0   , RotamerType::permutation , "t"  , 2 , 1 , {"monosaccharide"}       , {"monosaccharide"}                  , "C." , "O." , "C." , "C."  },
          // Omega angle in x-6 linkages.
          { "C.*"  , "O6"     , "Omg"  , -60.0  , AngleLimit{20.0, 20.0}  , 1.0   , RotamerType::permutation , "gg" , 3 , 1 , {}             , {"pyranose"}                  , "O6" , "C6" , "C5" , "O5"  }, // omg is O6-C5-C5-O5
          { "C.*"  , "O6"     , "Omg"  ,  60.0  , AngleLimit{20.0, 20.0}  , 1.0   , RotamerType::permutation , "gt" , 3 , 2 , {}             , {"pyranose"}                  , "O6" , "C6" , "C5" , "O5"  },
          { "C.*"  , "O6"     , "Omg"  , 180.0  , AngleLimit{20.0, 20.0}  , 1.0   , RotamerType::permutation , "tg" , 3 , 3 , {}             , {"pyranose"} 				 , "O6" , "C6" , "C5" , "O5"  },
          { "C.*"  , "O6"     , "Omg"  , 180.0  , AngleLimit{20.0, 20.0}  , 0.001 , RotamerType::permutation , "tg" , 3 , 3 , {}             , {"gauche-effect=gluco"}       , "O6" , "C6" , "C5" , "O5"  },
          // Omega angle in x-5 linkages.
          { "C.*"  , "O5"     , "Omg"  , -60.0  , AngleLimit{20.0, 20.0}  , 1.0   , RotamerType::permutation , "gg" , 3 , 1 , {}             , {"furanose"}                  , "O5" , "C5" , "C4" , "O4"  },
          { "C.*"  , "O5"     , "Omg"  ,  60.0  , AngleLimit{20.0, 20.0}  , 1.0   , RotamerType::permutation , "gt" , 3 , 2 , {}             , {"furanose"}                  , "O5" , "C5" , "C4" , "O4"  },
		  // Ketose 1-1 linkages. Copied from old GLYCAM-Web builder.
 		  { "C1"   , "O1"     , "Phi"  , 180.0  , AngleLimit{20.0, 20.0}  , 1.0   , RotamerType::permutation , ""   , 1 , 1 , {}             , {"ketose"}                , "C2" , "C1" , "O1" , "C1" },
  		  { "C1"   , "O1"     , "Psi"  ,-150.0  , AngleLimit{20.0, 20.0}  , 1.0   , RotamerType::permutation , ""   , 2 , 1 , {}             , {"ketose"}                , "C1" , "O1" , "C1" , "C2" },
  		  { "C1"   , "O1"     , "Omg"  ,  60.0  , AngleLimit{20.0, 20.0}  , 1.0   , RotamerType::permutation , "gt" , 3 , 1 , {}             , {"ketose"}                , "O1" , "C1" , "C2" , "O5" },
		  // Ketose 2-2 linkages to other ketoses
		  { "C2"   , "O2"     , "Phi"  ,    0.0  , AngleLimit{20.0, 20.0}  , 1.0   , RotamerType::permutation , ""   , 1 , 1 , {"ketose"}  , {"ketose"}                , "C1" , "C2" , "O2" , "C2" },
		  { "C2"   , "O2"     , "Psi"  ,   90.0  , AngleLimit{20.0, 20.0}  , 1.0   , RotamerType::permutation , ""   , 2 , 1 , {"ketose"}  , {"ketose"}                , "C2" , "O2" , "C2" , "C1" },
		  // Omega angle in x-1 linkages to ketofuranoses
          { "C.*"  , "O1"     , "Omg"  , -60.0  , AngleLimit{20.0, 20.0}  , 1.0   , RotamerType::permutation , "gg" , 3 , 1 , {}             , {"furanose", "ketose"}        , "O1" , "C1" , "C2" , "O6"  },
          { "C.*"  , "O1"     , "Omg"  ,  60.0  , AngleLimit{20.0, 20.0}  , 1.0   , RotamerType::permutation , "gt" , 3 , 2 , {}             , {"furanose", "ketose"}        , "O1" , "C1" , "C2" , "O6"  },
		  // Sucrose from: Bock K, Lemieux RU. Carbohydrate Research, 100 (1982) 63-74. Now allowing ordering of sequence either way, thus four entries instead of two.
		  { "C1"   , "O2"     , "Phi"  ,-135.0  , AngleLimit{ 5.0,  5.0}  , 1.0   , RotamerType::permutation , ""   , 1 , 1 , {"pyranose", "aldose"}    , {"furanose", "ketose"}         , "C2" , "C1" , "O2" , "C2" },
		  { "C1"   , "O2"     , "Psi"  ,  80.0  , AngleLimit{10.0, 10.0}  , 1.0   , RotamerType::permutation , ""   , 2 , 1 , {"pyranose", "aldose"}    , {"furanose", "ketose"}         , "C1" , "O2" , "C2" , "C1" },
		  { "C2"   , "O1"     , "Phi"  ,  80.0  , AngleLimit{10.0, 10.0}  , 1.0   , RotamerType::permutation , ""   , 1 , 1 , {"furanose", "ketose"}    , {"pyranose", "aldose"}         , "C1" , "C2" , "O1" , "C1" },
		  { "C2"   , "O1"     , "Psi"  ,-135.0  , AngleLimit{ 5.0,  5.0}  , 1.0   , RotamerType::permutation , ""   , 2 , 1 , {"furanose", "ketose"}    , {"pyranose", "aldose"}         , "C2" , "O1" , "C1" , "C2" },
		  // x-5 (and x-6) linkages to aldofuranoses, they are likely flexible but I have no data. Used values seen in 0LD prep file.
		  { "C.*"  , "O5"     , "Omg"  , 180.0  , AngleLimit{20.0, 20.0}  , 1.0   , RotamerType::permutation , "tg" , 3 , 1 , {}             , {"furanose"}                  , "O5" , "C5" , "C4" , "O4"  },
		  { "C.*"  , "O5"     , "Omg6" ,  60.0  , AngleLimit{20.0, 20.0}  , 1.0   , RotamerType::permutation , "gt" , 4 , 1 , {}             , {"furanose"}                  , "O6" , "C6" , "C5" , "C4"  },
		  // x-6
		  { "C.*"  , "O6"     , "Omg"  , -60.0  , AngleLimit{20.0, 20.0}  , 1.0   , RotamerType::permutation , "gg" , 3 , 1 , {}             , {"furanose"}                  , "O6" , "C6" , "C5" , "C4"  },
		  { "C.*"  , "O6"     , "Omg"  ,  60.0  , AngleLimit{20.0, 20.0}  , 1.0   , RotamerType::permutation , "gt" , 3 , 2 , {}             , {"furanose"}                  , "O6" , "C6" , "C5" , "C4"  },
		  { "C.*"  , "O6"     , "Omg"  , 180.0  , AngleLimit{20.0, 20.0}  , 1.0   , RotamerType::permutation , "tg" , 3 , 3 , {}             , {"furanose"}                  , "O6" , "C6" , "C5" , "C4"  },
		  { "C.*"  , "O6"     , "Omg4" ,  60.0  , AngleLimit{20.0, 20.0}  , 1.0   , RotamerType::permutation , "gt" , 4 , 1 , {}             , {"aldose", "furanose"}                  , "C6" , "C5" , "C4", "O4"  },
		  { "C.*"  , "O6"     , "Omg5" ,  60.0  , AngleLimit{20.0, 20.0}  , 1.0   , RotamerType::permutation , "gt" , 5 , 1 , {}             , {"aldose", "furanose"}                  , "O6" , "C6" , "C5" , "O5"  },
		  // 2-7 linkages copied from GlycamWeb Jan 2021. Branching in the linkage causes oddities with the bond number and which atoms are chosen for the torsion.
          { "C2"   , "O7"     , "Phi"  , -60.0  , AngleLimit{20.0, 20.0}  , 1.0   , RotamerType::permutation , "-g" , 1 , 2 , {"ulosonate"}  , {"ulosonate"}    , "C3" , "C2" , "O." , "C."  },
          { "C2"   , "O7"     , "Psi"  ,   0.0  , AngleLimit{20.0, 20.0}  , 1.0   , RotamerType::permutation , "c"  , 2 , 1 , {"ulosonate"}  , {"ulosonate"}    , "C." , "O." , "C." , "H."  },
          { "C.*"  , "O7"     , "Omg7" ,  60.0  , AngleLimit{20.0, 20.0}  , 1.0   , RotamerType::permutation , "g"  , 3 , 1 , {}             , {"monosaccharide"}    , "O7" , "C7" , "C6" , "O6"  },
          { "C.*"  , "O7"     , "Omg9" , 180.0  , AngleLimit{20.0, 20.0}  , 1.0   , RotamerType::permutation , "t"  , 4 , 1 , {}             , {"monosaccharide"}    , "O9" , "C9" , "C8" , "C7"  },
          { "C.*"  , "O7"     , "Omg8" ,  60.0  , AngleLimit{20.0, 20.0}  , 1.0   , RotamerType::permutation , "g"  , 5 , 1 , {}             , {"monosaccharide"}    , "C9" , "C8" , "C7" , "O7"  },
          // 2-9 linkages copied from GlycamWeb Jan 2021.
          { "C2"   , "O9"     , "Phi"  , -60.0  , AngleLimit{20.0, 20.0}  , 1.0   , RotamerType::permutation , "-g" , 1 , 2 , {"ulosonate"}  , {"ulosonate"}    , "C3" , "C2" , "O." , "C."  },
          { "C.*"  , "O9"     , "Omg9" , 180.0  , AngleLimit{20.0, 20.0}  , 1.0   , RotamerType::permutation , "t"  , 3 , 1 , {}                      , {"monosaccharide"}    , "O9" , "C9" , "C8" , "C7"  },
          { "C.*"  , "O9"     , "Omg8" , 180.0  , AngleLimit{20.0, 20.0}  , 1.0   , RotamerType::permutation , "t"  , 4 , 1 , {}                      , {"monosaccharide"}    , "C9" , "C8" , "C7" , "C6"  },
          { "C.*"  , "O9"     , "Omg7" , -60.0  , AngleLimit{20.0, 20.0}  , 1.0   , RotamerType::permutation , "-g" , 5 , 1 , {}                      , {"monosaccharide"}    , "O7" , "C7" , "C6" , "O6"  },
           // Generic x-8 linkages. Values copied from external conformer A below. Arbitrary, but I have no data for them.
          { "C.*"   , "O8"     , "Phi"  , -79.5  ,  AngleStd{20.0, 20.0}  , 0.42  , RotamerType::conformer   , "A"  , 1 , 1 , {"monosaccharide"}  , {"ulosonate"}    , "C1" , "C2" , "O8" , "C8"  },
          { "C.*"   , "O8"     , "Psi"  ,  88.1  ,  AngleStd{20.0, 20.0}  , 0.42  , RotamerType::conformer   , "A"  , 2 , 1 , {"monosaccharide"}  , {"ulosonate"}    , "C2" , "O8" , "C8" , "C7"  },
          { "C.*"   , "O8"     , "Omg8" ,-170.1  ,  AngleStd{20.0, 20.0}  , 0.42  , RotamerType::conformer   , "A"  , 3 , 1 , {"monosaccharide"}  , {"monosaccharide"}    , "O8" , "C8" , "C7" , "C6"  },
          { "C.*"   , "O8"     , "Omg7" , -61.8  ,  AngleStd{20.0, 20.0}  , 0.42  , RotamerType::conformer   , "A"  , 4 , 1 , {"monosaccharide"}  , {"monosaccharide"}    , "C8" , "C7" , "C6" , "O6"  },
          { "C.*"   , "O8"     , "Omg9" ,  65.8  ,  AngleStd{20.0, 20.0}  , 0.42  , RotamerType::conformer   , "A"  , 5 , 1 , {"monosaccharide"}  , {"monosaccharide"}    , "O9" , "C9" , "C8" , "O8"  },
// Oliver commenting out B-E as gems/website can't handle conformers properly yet. This is a fine stopgap to just create one shape.
		  // Internal 2-8 linkages
          { "C2"   , "O8"     , "Phi"  , -79.5  ,   AngleStd{20.0, 20.0}  , 0.42  , RotamerType::conformer   , "A"  , 1 , 1 , {"ulosonate", "internal"}  , {"ulosonate"}    , "C1" , "C2" , "O8" , "C8"  },
          { "C2"   , "O8"     , "Psi"  ,  88.1  ,   AngleStd{20.0, 20.0}  , 0.42  , RotamerType::conformer   , "A"  , 2 , 1 , {"ulosonate", "internal"}  , {"ulosonate"}    , "C2" , "O8" , "C8" , "C7"  },
          { "C2"   , "O8"     , "Omg8" ,-170.1  ,   AngleStd{20.0, 20.0}  , 0.42  , RotamerType::conformer   , "A"  , 3 , 1 , {"ulosonate", "internal"}  , {"monosaccharide"}    , "O8" , "C8" , "C7" , "C6"  },
          { "C2"   , "O8"     , "Omg7" , -61.8  ,   AngleStd{20.0, 20.0}  , 0.42  , RotamerType::conformer   , "A"  , 4 , 1 , {"ulosonate", "internal"}  , {"monosaccharide"}    , "C8" , "C7" , "C6" , "O6"  },
          { "C2"   , "O8"     , "Omg9" ,  65.8  ,   AngleStd{20.0, 20.0}  , 0.42  , RotamerType::conformer   , "A"  , 5 , 1 , {"ulosonate", "internal"}  , {"monosaccharide"}    , "O9" , "C9" , "C8" , "O8"  },
//          { "C2"   , "O8"     , "Phi"  , -73.6  ,   AngleStd{20.0, 20.0}  , 0.24  , RotamerType::conformer   , "B"  , 1 , 2 , {"ulosonate", "alpha", "internal"}  , {"ulosonate"}    , "C1" , "C2" , "O8" , "C8"  },
//          { "C2"   , "O8"     , "Psi"  , 109.2  ,   AngleStd{20.0, 20.0}  , 0.24  , RotamerType::conformer   , "B"  , 2 , 2 , {"ulosonate", "alpha", "internal"}  , {"ulosonate"}    , "C2" , "O8" , "C8" , "C7"  },
//          { "C2"   , "O8"     , "Omg8" ,-169.2  ,   AngleStd{20.0, 20.0}  , 0.24  , RotamerType::conformer   , "B"  , 3 , 2 , {"ulosonate", "alpha", "internal"}  , {"ulosonate"}    , "O8" , "C8" , "C7" , "C6"  },
//          { "C2"   , "O8"     , "Omg7" , -61.9  ,   AngleStd{20.0, 20.0}  , 0.24  , RotamerType::conformer   , "B"  , 4 , 2 , {"ulosonate", "alpha", "internal"}  , {"ulosonate"}    , "C8" , "C7" , "C6" , "O6"  },
//          { "C2"   , "O8"     , "Omg9" , -61.0  ,   AngleStd{20.0, 20.0}  , 0.24  , RotamerType::conformer   , "B"  , 5 , 2 , {"ulosonate", "alpha", "internal"}  , {"ulosonate"}    , "O9" , "C9" , "C8" , "O8"  },
//          { "C2"   , "O8"     , "Phi"  ,-172.0  ,   AngleStd{20.0, 20.0}  , 0.08  , RotamerType::conformer   , "C"  , 1 , 3 , {"ulosonate", "alpha", "internal"}  , {"ulosonate"}    , "C1" , "C2" , "O8" , "C8"  },
//          { "C2"   , "O8"     , "Psi"  , 108.7  ,   AngleStd{20.0, 20.0}  , 0.08  , RotamerType::conformer   , "C"  , 2 , 3 , {"ulosonate", "alpha", "internal"}  , {"ulosonate"}    , "C2" , "O8" , "C8" , "C7"  },
//          { "C2"   , "O8"     , "Omg8" ,-164.2  ,   AngleStd{20.0, 20.0}  , 0.08  , RotamerType::conformer   , "C"  , 3 , 3 , {"ulosonate", "alpha", "internal"}  , {"ulosonate"}    , "O8" , "C8" , "C7" , "C6"  },
//          { "C2"   , "O8"     , "Omg7" , -52.9  ,   AngleStd{20.0, 20.0}  , 0.08  , RotamerType::conformer   , "C"  , 4 , 3 , {"ulosonate", "alpha", "internal"}  , {"ulosonate"}    , "C8" , "C7" , "C6" , "O6"  },
//          { "C2"   , "O8"     , "Omg9" ,  68.5  ,   AngleStd{20.0, 20.0}  , 0.08  , RotamerType::conformer   , "C"  , 5 , 3 , {"ulosonate", "alpha", "internal"}  , {"ulosonate"}    , "O9" , "C9" , "C8" , "O8"  },
//          { "C2"   , "O8"     , "Phi"  , -72.1  ,   AngleStd{20.0, 20.0}  , 0.07  , RotamerType::conformer   , "D"  , 1 , 4 , {"ulosonate", "alpha", "internal"}  , {"ulosonate"}    , "C1" , "C2" , "O8" , "C8"  },
//          { "C2"   , "O8"     , "Psi"  , 127.1  ,   AngleStd{20.0, 20.0}  , 0.07  , RotamerType::conformer   , "D"  , 2 , 4 , {"ulosonate", "alpha", "internal"}  , {"ulosonate"}    , "C2" , "O8" , "C8" , "C7"  },
//          { "C2"   , "O8"     , "Omg8" ,  58.5  ,   AngleStd{20.0, 20.0}  , 0.07  , RotamerType::conformer   , "D"  , 3 , 4 , {"ulosonate", "alpha", "internal"}  , {"ulosonate"}    , "O8" , "C8" , "C7" , "C6"  },
//          { "C2"   , "O8"     , "Omg7" , -66.1  ,   AngleStd{20.0, 20.0}  , 0.07  , RotamerType::conformer   , "D"  , 4 , 4 , {"ulosonate", "alpha", "internal"}  , {"ulosonate"}    , "C8" , "C7" , "C6" , "O6"  },
//          { "C2"   , "O8"     , "Omg9" , 178.9  ,   AngleStd{20.0, 20.0}  , 0.07  , RotamerType::conformer   , "D"  , 5 , 4 , {"ulosonate", "alpha", "internal"}  , {"ulosonate"}    , "O9" , "C9" , "C8" , "O8"  },
//          { "C2"   , "O8"     , "Phi"  , -42.9  ,   AngleStd{20.0, 20.0}  , 0.07  , RotamerType::conformer   , "E"  , 1 , 5 , {"ulosonate", "alpha", "internal"}  , {"ulosonate"}    , "C1" , "C2" , "O8" , "C8"  },
//          { "C2"   , "O8"     , "Psi"  , 162.3  ,   AngleStd{20.0, 20.0}  , 0.07  , RotamerType::conformer   , "E"  , 2 , 5 , {"ulosonate", "alpha", "internal"}  , {"ulosonate"}    , "C2" , "O8" , "C8" , "C7"  },
//          { "C2"   , "O8"     , "Omg8" , -58.3  ,   AngleStd{20.0, 20.0}  , 0.07  , RotamerType::conformer   , "E"  , 3 , 5 , {"ulosonate", "alpha", "internal"}  , {"ulosonate"}    , "O8" , "C8" , "C7" , "C6"  },
//          { "C2"   , "O8"     , "Omg7" , -55.7  ,   AngleStd{20.0, 20.0}  , 0.07  , RotamerType::conformer   , "E"  , 4 , 5 , {"ulosonate", "alpha", "internal"}  , {"ulosonate"}    , "C8" , "C7" , "C6" , "O6"  },
//          { "C2"   , "O8"     , "Omg9" , -58.5  ,   AngleStd{20.0, 20.0}  , 0.07  , RotamerType::conformer   , "E"  , 5 , 5 , {"ulosonate", "alpha", "internal"}  , {"ulosonate"}    , "O9" , "C9" , "C8" , "O8"  },
//          // External 2-8 linkages
          { "C2"   , "O8"     , "Phi"  , -66.0  ,   AngleStd{20.0, 20.0}  , 0.42  , RotamerType::conformer   , "A"  , 1 , 1 , {"ulosonate", "external"}  , {"ulosonate"}    , "C1" , "C2" , "O8" , "C8"  },
          { "C2"   , "O8"     , "Psi"  , 118.2  ,   AngleStd{20.0, 20.0}  , 0.42  , RotamerType::conformer   , "A"  , 2 , 1 , {"ulosonate", "external"}  , {"ulosonate"}    , "C2" , "O8" , "C8" , "C7"  },
          { "C2"   , "O8"     , "Omg8" ,-160.4  ,   AngleStd{20.0, 20.0}  , 0.42  , RotamerType::conformer   , "A"  , 3 , 1 , {"ulosonate", "external"}  , {"ulosonate"}    , "O8" , "C8" , "C7" , "C6"  },
          { "C2"   , "O8"     , "Omg7" , -59.6  ,   AngleStd{20.0, 20.0}  , 0.42  , RotamerType::conformer   , "A"  , 4 , 1 , {"ulosonate", "external"}  , {"ulosonate"}    , "C8" , "C7" , "C6" , "O6"  },
          { "C2"   , "O8"     , "Omg9" ,  73.4  ,   AngleStd{20.0, 20.0}  , 0.42  , RotamerType::conformer   , "A"  , 5 , 1 , {"ulosonate", "external"}  , {"ulosonate"}    , "O9" , "C9" , "C8" , "O8"  },
//          { "C2"   , "O8"     , "Phi"  , -68.4  ,   AngleStd{20.0, 20.0}  , 0.24  , RotamerType::conformer   , "B"  , 1 , 2 , {"ulosonate", "alpha", "external"}  , {"ulosonate"}    , "C1" , "C2" , "O8" , "C8"  },
//          { "C2"   , "O8"     , "Psi"  , 132.7  ,   AngleStd{20.0, 20.0}  , 0.24  , RotamerType::conformer   , "B"  , 2 , 2 , {"ulosonate", "alpha", "external"}  , {"ulosonate"}    , "C2" , "O8" , "C8" , "C7"  },
//          { "C2"   , "O8"     , "Omg8" ,  64.4  ,   AngleStd{20.0, 20.0}  , 0.24  , RotamerType::conformer   , "B"  , 3 , 2 , {"ulosonate", "alpha", "external"}  , {"ulosonate"}    , "O8" , "C8" , "C7" , "C6"  },
//          { "C2"   , "O8"     , "Omg7" , -61.1  ,   AngleStd{20.0, 20.0}  , 0.24  , RotamerType::conformer   , "B"  , 4 , 2 , {"ulosonate", "alpha", "external"}  , {"ulosonate"}    , "C8" , "C7" , "C6" , "O6"  },
//          { "C2"   , "O8"     , "Omg9" ,-179.0  ,   AngleStd{20.0, 20.0}  , 0.24  , RotamerType::conformer   , "B"  , 5 , 2 , {"ulosonate", "alpha", "external"}  , {"ulosonate"}    , "O9" , "C9" , "C8" , "O8"  },
//          { "C2"   , "O8"     , "Phi"  , -53.2  ,   AngleStd{20.0, 20.0}  , 0.08  , RotamerType::conformer   , "C"  , 1 , 3 , {"ulosonate", "alpha", "external"}  , {"ulosonate"}    , "C1" , "C2" , "O8" , "C8"  },
//          { "C2"   , "O8"     , "Psi"  , 161.6  ,   AngleStd{20.0, 20.0}  , 0.08  , RotamerType::conformer   , "C"  , 2 , 3 , {"ulosonate", "alpha", "external"}  , {"ulosonate"}    , "C2" , "O8" , "C8" , "C7"  },
//          { "C2"   , "O8"     , "Omg8" , -67.8  ,   AngleStd{20.0, 20.0}  , 0.08  , RotamerType::conformer   , "C"  , 3 , 3 , {"ulosonate", "alpha", "external"}  , {"ulosonate"}    , "O8" , "C8" , "C7" , "C6"  },
//          { "C2"   , "O8"     , "Omg7" , -61.1  ,   AngleStd{20.0, 20.0}  , 0.08  , RotamerType::conformer   , "C"  , 4 , 3 , {"ulosonate", "alpha", "external"}  , {"ulosonate"}    , "C8" , "C7" , "C6" , "O6"  },
//          { "C2"   , "O8"     , "Omg9" , -57.9  ,   AngleStd{20.0, 20.0}  , 0.08  , RotamerType::conformer   , "C"  , 5 , 3 , {"ulosonate", "alpha", "external"}  , {"ulosonate"}    , "O9" , "C9" , "C8" , "O8"  },
//          { "C2"   , "O8"     , "Phi"  , -72.7  ,   AngleStd{20.0, 20.0}  , 0.07  , RotamerType::conformer   , "D"  , 1 , 4 , {"ulosonate", "alpha", "external"}  , {"ulosonate"}    , "C1" , "C2" , "O8" , "C8"  },
//          { "C2"   , "O8"     , "Psi"  , 124.9  ,   AngleStd{20.0, 20.0}  , 0.07  , RotamerType::conformer   , "D"  , 2 , 4 , {"ulosonate", "alpha", "external"}  , {"ulosonate"}    , "C2" , "O8" , "C8" , "C7"  },
//          { "C2"   , "O8"     , "Omg8" ,  62.4  ,   AngleStd{20.0, 20.0}  , 0.07  , RotamerType::conformer   , "D"  , 3 , 4 , {"ulosonate", "alpha", "external"}  , {"ulosonate"}    , "O8" , "C8" , "C7" , "C6"  },
//          { "C2"   , "O8"     , "Omg7" , -59.5  ,   AngleStd{20.0, 20.0}  , 0.07  , RotamerType::conformer   , "D"  , 4 , 4 , {"ulosonate", "alpha", "external"}  , {"ulosonate"}    , "C8" , "C7" , "C6" , "O6"  },
//          { "C2"   , "O8"     , "Omg9" ,  78.8  ,   AngleStd{20.0, 20.0}  , 0.07  , RotamerType::conformer   , "D"  , 5 , 4 , {"ulosonate", "alpha", "external"}  , {"ulosonate"}    , "O9" , "C9" , "C8" , "O8"  },
//          { "C2"   , "O8"     , "Phi"  ,-164.9  ,   AngleStd{20.0, 20.0}  , 0.07  , RotamerType::conformer   , "E"  , 1 , 5 , {"ulosonate", "alpha", "external"}  , {"ulosonate"}    , "C1" , "C2" , "O8" , "C8"  },
//          { "C2"   , "O8"     , "Psi"  , 103.0  ,   AngleStd{20.0, 20.0}  , 0.07  , RotamerType::conformer   , "E"  , 2 , 5 , {"ulosonate", "alpha", "external"}  , {"ulosonate"}    , "C2" , "O8" , "C8" , "C7"  },
//          { "C2"   , "O8"     , "Omg8" ,-161.3  ,   AngleStd{20.0, 20.0}  , 0.07  , RotamerType::conformer   , "E"  , 3 , 5 , {"ulosonate", "alpha", "external"}  , {"ulosonate"}    , "O8" , "C8" , "C7" , "C6"  },
//          { "C2"   , "O8"     , "Omg7" , -54.3  ,   AngleStd{20.0, 20.0}  , 0.07  , RotamerType::conformer   , "E"  , 4 , 5 , {"ulosonate", "alpha", "external"}  , {"ulosonate"}    , "C8" , "C7" , "C6" , "O6"  },
//          { "C2"   , "O8"     , "Omg9" ,  71.1  ,   AngleStd{20.0, 20.0}  , 0.07  , RotamerType::conformer   , "E"  , 5 , 5 , {"ulosonate", "alpha", "external"}  , {"ulosonate"}    , "O9" , "C9" , "C8" , "O8"  },

        // Common sugar derivatives
        // Phosphate/sulfate
          { "[SP]1", "N[2]"   , "Phi"  , 180.0  , AngleLimit{20.0, 20.0}  , 1.0   , RotamerType::permutation , ""   , 1 , 1 , {"derivative"}       , {"monosaccharide"}                  , "O." , ".1" , "N." , "C."  },
          { "[SP]1", "N[2]"   , "Psi"  , -40.0  , AngleLimit{20.0, 20.0}  , 1.0   , RotamerType::permutation , ""   , 2 , 1 , {"derivative"}       , {"monosaccharide"}                  , ".1" , "N." , "C." , "H."  },
          { "[SP]1", "O[1-9]" , "Phi"  , 180.0  , AngleLimit{20.0, 20.0}  , 1.0   , RotamerType::permutation , ""   , 1 , 1 , {"derivative"}       , {"monosaccharide"}                  , "O." , ".1" , "O." , "C."  },
          { "[SP]1", "O[1-5]" , "Psi"  ,   0.0  , AngleLimit{20.0, 20.0}  , 1.0   , RotamerType::permutation , ""   , 2 , 1 , {"derivative"}       , {"monosaccharide"}                  , ".1" , "O." , "C." , "H."  },
          { "[SP]1", "O[6-9]" , "Psi"  , 180.0  , AngleLimit{20.0, 20.0}  , 1.0   , RotamerType::permutation , ""   , 2 , 1 , {"derivative"}       , {"monosaccharide"}                  , ".1" , "O." , "C." , "C."  },
          { "[SP]1", "O6"     , "Omg"  , -60.0  , AngleLimit{20.0, 20.0}  , 1.0   , RotamerType::permutation , "gg" , 3 , 1 , {"derivative"}       , {"monosaccharide"}                  , "O6" , "C6" , "C5" , "O5"  },
  //        These two were removed to stop them showing up as options the carb builder, but for some applications like grafting you may want them enabled. No functionality for that context yet tho.
  //        { "[SP]1", "O6"     , "Omg"  ,  60.0  ,  20.0  ,  20.0  , 1.0   , RotamerType::permutation , "gt" , 3 , 2 , {"derivative"}       , {"monosaccharide"}                  , "O6" , "C6" , "C5" , "O5"  },
  //        { "[SP]1", "O6"     , "Omg"  , 180.0  ,  20.0  ,  20.0  , 1.0   , RotamerType::permutation , "tg" , 3 , 3 , {"derivative"}       , {"gauche-effect=galacto"} , "O6" , "C6" , "C5" , "O5"  },
        // Ac ACX
          { "C1A"  , "O[1-9]" , "Phi"  , 180.0  , AngleLimit{20.0, 20.0}  , 1.0   , RotamerType::permutation , "t"  , 1 , 1 , {"derivative"}       , {"monosaccharide"}                  , "C2A", "C1A", "O." , "C."  },
          { "C1A"  , "O[1-5]" , "Psi"  ,   0.0  , AngleLimit{20.0, 20.0}  , 1.0   , RotamerType::permutation , "c"  , 2 , 1 , {"derivative"}       , {"monosaccharide"}                  , "C1A", "O." , "C." , "H."  },
          { "C1A"  , "O[6-9]" , "Psi"  , 180.0  , AngleLimit{20.0, 20.0}  , 1.0   , RotamerType::permutation , "-g" , 2 , 1 , {"derivative"}       , {"monosaccharide"}                  , "C1A", "O." , "C." , "C."  },
        // Me MEX
          { "CH3"  , "O[1-9]" , "Phi"  , 180.0  , AngleLimit{20.0, 20.0}  , 1.0   , RotamerType::permutation , "t"  , 1 , 1 , {"derivative"}       , {"monosaccharide"}                  , "CH3", "O." , "C." , "H."  },
          { "CH3"  , "O[1-5]" , "Psi"  ,   0.0  , AngleLimit{20.0, 20.0}  , 1.0   , RotamerType::permutation , "c"  , 2 , 1 , {"derivative"}       , {"monosaccharide"}                  , "CH3", "O." , "C." , "H."  },
          { "CH3"  , "O[6-9]" , "Psi"  , 180.0  , AngleLimit{20.0, 20.0}  , 1.0   , RotamerType::permutation , "-g" , 2 , 1 , {"derivative"}       , {"monosaccharide"}                  , "CH3", "O." , "C." , "C."  },
        // Any derivative linked to O8 of Sia. Copied from 2-8 linkage values.
          { ".*"   , "O8"     , "Omg8" ,-160.4  , AngleLimit{20.0, 20.0}  , 1.0   , RotamerType::permutation , ""   , 3 , 1 , {"derivative"}                        , {"ulosonate"}    , "O8" , "C8" , "C7" , "C6"  },
          { ".*"   , "O8"     , "Omg7" , -59.6  , AngleLimit{20.0, 20.0}  , 1.0   , RotamerType::permutation , ""   , 4 , 1 , {"derivative"}                        , {"ulosonate"}    , "C8" , "C7" , "C6" , "O6"  },
          { ".*"   , "O8"     , "Omg9" ,  73.4  , AngleLimit{20.0, 20.0}  , 1.0   , RotamerType::permutation , ""   , 5 , 1 , {"derivative"}                        , {"ulosonate"}    , "O9" , "C9" , "C8" , "O8"  },
          // Protein linkages
          // ASN // Values are from Petrescu et al 2004.
          { "C."   , "ND2"    , "Chi1" , 191.6  ,   AngleStd{14.4, 14.4}  , 0.497 , RotamerType::conformer   , "A"  , 4 , 1 , {}             , {"amino-acid"}            , "CG" , "CB" , "CA" , "N"   },
          { "C."   , "ND2"    , "Chi2" , 177.6  ,   AngleStd{43.0, 43.0}  , 0.497 , RotamerType::conformer   , "A"  , 3 , 1 , {}             , {"amino-acid"}            , "ND2", "CG" , "CB" , "CA"  },
          { "C."   , "ND2"    , "Psi"  , 177.3  ,   AngleStd{12.3, 12.3}  , 0.497 , RotamerType::conformer   , "A"  , 2 , 1 , {}             , {"amino-acid"}            , "C." , "ND2", "CG" , "CB"  },
          { "C1"   , "ND2"    , "Phi"  , 141.0  ,   AngleStd{21.3, 21.3}  , 0.497 , RotamerType::conformer   , "A"  , 1 , 1 , {}             , {"amino-acid"}            , "C." , "C." , "ND2", "CG"  },
          { "C."   , "ND2"    , "Chi1" ,  63.6  ,   AngleStd{ 8.9,  8.9}  , 0.178 , RotamerType::conformer   , "B"  , 4 , 2 , {}             , {"amino-acid"}            , "CG" , "CB" , "CA" , "N"   },
          { "C."   , "ND2"    , "Chi2" , 191.1  ,   AngleStd{31.6, 31.6}  , 0.178 , RotamerType::conformer   , "B"  , 3 , 2 , {}             , {"amino-acid"}            , "ND2", "CG" , "CB" , "CA"  },
          { "C."   , "ND2"    , "Psi"  , 178.5  ,   AngleStd{13.9, 13.9}  , 0.178 , RotamerType::conformer   , "B"  , 2 , 2 , {}             , {"amino-acid"}            , "C." , "ND2", "CG" , "CB"  },
          { "C1"   , "ND2"    , "Phi"  , 133.7  ,   AngleStd{21.5, 21.5}  , 0.178 , RotamerType::conformer   , "B"  , 1 , 2 , {}             , {"amino-acid"}            , "C." , "C." , "ND2", "CG"  },
          { "C."   , "ND2"    , "Chi1" , 290.6  ,   AngleStd{12.7, 12.7}  , 0.235 , RotamerType::conformer   , "C"  , 4 , 3 , {}             , {"amino-acid"}            , "CG" , "CB" , "CA" , "N"   },
          { "C."   , "ND2"    , "Chi2" , 152.9  ,   AngleStd{23.9, 23.9}  , 0.235 , RotamerType::conformer   , "C"  , 3 , 3 , {}             , {"amino-acid"}            , "ND2", "CG" , "CB" , "CA"  },
          { "C."   , "ND2"    , "Psi"  , 173.1  ,   AngleStd{12.2, 12.2}  , 0.235 , RotamerType::conformer   , "C"  , 2 , 3 , {}             , {"amino-acid"}            , "C." , "ND2", "CG" , "CB"  },
          { "C1"   , "ND2"    , "Phi"  , 148.0  ,   AngleStd{20.3, 20.3}  , 0.235 , RotamerType::conformer   , "C"  , 1 , 3 , {}             , {"amino-acid"}            , "C." , "C." , "ND2", "CG"  },
          { "C."   , "ND2"    , "Chi1" , 302.3  ,   AngleStd{11.5, 11.5}  , 0.090 , RotamerType::conformer   , "D"  , 4 , 4 , {}             , {"amino-acid"}            , "CG" , "CB" , "CA" , "N"   },
          { "C."   , "ND2"    , "Chi2" , 255.0  ,   AngleStd{28.8, 28.8}  , 0.090 , RotamerType::conformer   , "D"  , 3 , 4 , {}             , {"amino-acid"}            , "ND2", "CG" , "CB" , "CA"  },
          { "C."   , "ND2"    , "Psi"  , 178.1  ,   AngleStd{11.5, 11.5}  , 0.090 , RotamerType::conformer   , "D"  , 2 , 4 , {}             , {"amino-acid"}            , "C." , "ND2", "CG" , "CB"  },
          { "C1"   , "ND2"    , "Phi"  , 147.5  ,   AngleStd{23.9, 23.9}  , 0.090 , RotamerType::conformer   , "D"  , 1 , 4 , {}             , {"amino-acid"}            , "C." , "C." , "ND2", "CG"  },
          // THR // Values are from Lovell et al "PENULTIMATE ROTAMER LIBRARY"
          { "C."   , "OG1"    , "Chi1" ,  59.0  ,   AngleStd{10.0, 10.0}  , 0.490 , RotamerType::conformer   , "p"  , 3 , 1 , {}             , {"amino-acid"}            , "OG1", "CB" , "CA" , "N"   },
          { "C."   , "OG1"    , "Psi"  , -60.0  ,   AngleStd{20.0, 20.0}  , 0.490 , RotamerType::conformer   , "p"  , 2 , 1 , {}             , {"amino-acid"}            , "C." , "OG1", "CB" , "CA"  },
          { "C."   , "OG1"    , "Phi"  , 180.0  ,   AngleStd{20.0, 20.0}  , 0.490 , RotamerType::conformer   , "p"  , 1 , 1 , {}             , {"amino-acid"}            , "C." , "C." , "OG1", "CB"  },
		  { "C."   , "OG1"    , "Chi1" ,-171.0  ,   AngleStd{ 6.0,  6.0}  , 0.070 , RotamerType::conformer   , "t"  , 3 , 2 , {}             , {"amino-acid"}            , "OG1", "CB" , "CA" , "N"   },
		  { "C."   , "OG1"    , "Psi"  , -60.0  ,   AngleStd{20.0, 20.0}  , 0.070 , RotamerType::conformer   , "t"  , 2 , 2 , {}             , {"amino-acid"}            , "C." , "OG1", "CB" , "CA"  },
		  { "C."   , "OG1"    , "Phi"  , 180.0  ,   AngleStd{20.0, 20.0}  , 0.070 , RotamerType::conformer   , "t"  , 1 , 2 , {}             , {"amino-acid"}            , "C." , "C." , "OG1", "CB"  },
	      { "C."   , "OG1"    , "Chi1" , -61.0  ,   AngleStd{ 7.0,  7.0}  , 0.430 , RotamerType::conformer   , "m"  , 3 , 3 , {}             , {"amino-acid"}            , "OG1", "CB" , "CA" , "N"   },
		  { "C."   , "OG1"    , "Psi"  , -60.0  ,   AngleStd{20.0, 20.0}  , 0.430 , RotamerType::conformer   , "m"  , 2 , 3 , {}             , {"amino-acid"}            , "C." , "OG1", "CB" , "CA"  },
		  { "C."   , "OG1"    , "Phi"  , 180.0  ,   AngleStd{20.0, 20.0}  , 0.430 , RotamerType::conformer   , "m"  , 1 , 3 , {}             , {"amino-acid"}            , "C." , "C." , "OG1", "CB"  },
           // SER // Values are from Lovell et al "PENULTIMATE ROTAMER LIBRARY"
          { "C."   , "OG"     , "Chi1" ,  64.0  ,   AngleStd{10.0, 10.0}  , 0.480 , RotamerType::conformer   , "p"  , 3 , 1 , {}             , {"amino-acid"}            , "OG" , "CB" , "CA" , "N"   },
          { "C."   , "OG"     , "Psi"  ,   0.0  ,   AngleStd{20.0, 20.0}  , 0.480 , RotamerType::conformer   , "p"  , 2 , 1 , {}             , {"amino-acid"}            , "C." , "OG" , "CB" , "CA"  },
          { "C."   , "OG"     , "Phi"  , 180.0  ,   AngleStd{20.0, 20.0}  , 0.480 , RotamerType::conformer   , "p"  , 1 , 1 , {}             , {"amino-acid"}            , "OG" , "CB" , "CA" , "N"   },
          { "C."   , "OG"     , "Chi1" , 178.0  ,   AngleStd{11.0, 11.0}  , 0.220 , RotamerType::conformer   , "t"  , 3 , 2 , {}             , {"amino-acid"}            , "OG" , "CB" , "CA" , "N"   },
          { "C."   , "OG"     , "Psi"  ,   0.0  ,   AngleStd{20.0, 20.0}  , 0.220 , RotamerType::conformer   , "t"  , 2 , 2 , {}             , {"amino-acid"}            , "C." , "OG" , "CB" , "CA"  },
          { "C."   , "OG"     , "Phi"  , 180.0  ,   AngleStd{20.0, 20.0}  , 0.220 , RotamerType::conformer   , "t"  , 1 , 2 , {}             , {"amino-acid"}            , "OG" , "CB" , "CA" , "N"   },
          { "C."   , "OG"     , "Chi1" , -65.0  ,   AngleStd{ 9.0,  9.0}  , 0.290 , RotamerType::conformer   , "m"  , 3 , 3 , {}             , {"amino-acid"}            , "OG" , "CB" , "CA" , "N"   },
          { "C."   , "OG"     , "Psi"  ,   0.0  ,   AngleStd{20.0, 20.0}  , 0.290 , RotamerType::conformer   , "m"  , 2 , 3 , {}             , {"amino-acid"}            , "C." , "OG" , "CB" , "CA"  },
          { "C."   , "OG"     , "Phi"  , 180.0  ,   AngleStd{20.0, 20.0}  , 0.290 , RotamerType::conformer   , "m"  , 1 , 3 , {}             , {"amino-acid"}            , "C." , "C." , "OG1", "CB"  },
           // TYR // Values are from Lovell et al "PENULTIMATE ROTAMER LIBRARY"
          { "C."   , "OH"     , "Chi1" ,  63.0  ,   AngleStd{13.0, 13.0}  , 0.130 , RotamerType::conformer   , "p90", 7 , 1 , {}             , {"amino-acid"}            , "CG" , "CB" , "CA" , "N"   },
          { "C."   , "OH"     , "Chi2" ,  89.0  ,   AngleStd{13.0, 13.0}  , 0.130 , RotamerType::conformer   , "p90", 6 , 1 , {}             , {"amino-acid"}            , "CG" , "CB" , "CA" , "N"   },
		  { "C."   , "OH"     , "Psi"  , -60.0  ,   AngleStd{20.0, 20.0}  , 0.130 , RotamerType::conformer   , "p90", 2 , 1 , {}             , {"amino-acid"}            , "C." , "OH" , "CZ" , "CE1" },
		  { "C."   , "OH"     , "Phi"  , 180.0  ,   AngleStd{20.0, 20.0}  , 0.130 , RotamerType::conformer   , "p90", 1 , 1 , {}             , {"amino-acid"}            , "C." , "C." , "OH ", "CZ"  },
		  { "C."   , "OH"     , "Chi1" ,  63.0  ,   AngleStd{11.0, 11.0}  , 0.340 , RotamerType::conformer   , "t80", 7 , 2 , {}             , {"amino-acid"}            , "CG" , "CB" , "CA" , "N"   },
          { "C."   , "OH"     , "Chi2" ,  89.0  ,   AngleStd{14.0, 14.0}  , 0.340 , RotamerType::conformer   , "t80", 6 , 2 , {}             , {"amino-acid"}            , "CG" , "CB" , "CA" , "N"   },
		  { "C."   , "OH"     , "Psi"  , -60.0  ,   AngleStd{20.0, 20.0}  , 0.340 , RotamerType::conformer   , "t80", 2 , 2 , {}             , {"amino-acid"}            , "C." , "OH" , "CZ" , "CE1" },
		  { "C."   , "OH"     , "Phi"  , 180.0  ,   AngleStd{20.0, 20.0}  , 0.340 , RotamerType::conformer   , "t80", 1 , 2 , {}             , {"amino-acid"}            , "C." , "C." , "OH ", "CZ"  },
		  { "C."   , "OH"     , "Chi1" , -65.0  ,   AngleStd{11.0, 11.0}  , 0.430 , RotamerType::conformer   ,"m-85", 7 , 3 , {}             , {"amino-acid"}            , "CG" , "CB" , "CA" , "N"   },
		  { "C."   , "OH"     , "Chi2" , -87.0  ,   AngleStd{21.0, 21.0}  , 0.430 , RotamerType::conformer   ,"m-85", 6 , 3 , {}             , {"amino-acid"}            , "CG" , "CB" , "CA" , "N"   },
		  { "C."   , "OH"     , "Psi"  , -60.0  ,   AngleStd{20.0, 20.0}  , 0.430 , RotamerType::conformer   ,"m-85", 2 , 3 , {}             , {"amino-acid"}            , "C." , "OH" , "CZ" , "CE1" },
		  { "C."   , "OH"     , "Phi"  , 180.0  ,   AngleStd{20.0, 20.0}  , 0.430 , RotamerType::conformer   ,"m-85", 1 , 3 , {}             , {"amino-acid"}            , "C." , "C." , "OH ", "CZ"  },
		  { "C."   , "OH"     , "Chi1" , -64.0  ,   AngleStd{11.0, 11.0}  , 0.090 , RotamerType::conformer   ,"m-30", 7 , 4 , {}             , {"amino-acid"}            , "CG" , "CB" , "CA" , "N"   },
		  { "C."   , "OH"     , "Chi2" , -30.0  ,   AngleStd{18.0, 18.0}  , 0.090 , RotamerType::conformer   ,"m-30", 6 , 4 , {}             , {"amino-acid"}            , "CG" , "CB" , "CA" , "N"   },
		  { "C."   , "OH"     , "Psi"  , -60.0  ,   AngleStd{20.0, 20.0}  , 0.090 , RotamerType::conformer   ,"m-30", 2 , 4 , {}             , {"amino-acid"}            , "C." , "OH" , "CZ" , "CE1" },
		  { "C."   , "OH"     , "Phi"  , 180.0  ,   AngleStd{20.0, 20.0}  , 0.090 , RotamerType::conformer   ,"m-30", 1 , 4 , {}             , {"amino-acid"}            , "C." , "C." , "OH ", "CZ"  },

          // ROH
          { "C1"   , "O1"     , "Phi"  , 180.0  , AngleLimit{20.0, 20.0}  , 1.0   , RotamerType::permutation , "g"  , 1 , 1 , {}           , {"aglycon"}               , "C2" , "C1" , "O1" , "H1"  },
          { "C2"   , "O1"     , "Phi"  , 180.0  , AngleLimit{20.0, 20.0}  , 1.0   , RotamerType::permutation , "g"  , 1 , 1 , {}           , {"aglycon"}               , "C2" , "C1" , "O1" , "H1"  },
          { "C1"   , "O"      , "Phi"  , 180.0  , AngleLimit{20.0, 20.0}  , 1.0   , RotamerType::permutation , "g"  , 1 , 1 , {}           , {"aglycon"}               , "C2" , "C1" , "O1" , "H1"  },
          { "C2"   , "O"      , "Phi"  , 180.0  , AngleLimit{20.0, 20.0}  , 1.0   , RotamerType::permutation , "g"  , 1 , 1 , {}           , {"aglycon"}               , "C2" , "C1" , "O1" , "H1"  }}};

    // clang-format on

    //    Statistical analysis of the protein environment of N-glycosylation sites: implications for occupancy,
    //    structure, and folding Andrei-J. Petrescu  Adina-L. Milac  Stefana M. Petrescu  Raymond A. Dwek Mark R.
    //    Wormald Glycobiology, Volume 14, Issue 2, 1 February 2004, Pages 103–114,

    // Some entries have conditions for the first or second residue to have a particular type (aka tag).
    // Most entries have "none" for condition. This checks first if condition is "none", and therefore satisfied.
    // Otherwise (else if) it checks if any of the residue_types match the condition for the entry, e.g.
    // gauche_effect=galacto.

    std::function<double(const DihedralAngleData&)> metadataWeight = [](const DihedralAngleData& entry)
    {
        return entry.weight_;
    };

    const GlycamMetadata::DihedralAngleDataTable dihedralAngleDataTable_ {
        dihedralAngleDataVector_, codeUtils::vectorMap(metadataWeight, dihedralAngleDataVector_)};

    bool checkIfResidueConditionsAreSatisfied(const std::vector<std::string>& residue_types,
                                              const std::vector<std::string>& entry_conditions)
    {
        for (auto& entry_condition : entry_conditions)
        {
            if (!codeUtils::contains(residue_types, entry_condition))
            {
                return false;
            }
        }
        return true;
    }
} // namespace

const GlycamMetadata::DihedralAngleDataTable& GlycamMetadata::dihedralAngleDataTable()
{
    return dihedralAngleDataTable_;
}

// All of this is rubbish but I want to match what's already in the metadata (strings).
// Much better to e.g. put the enum type into the metadata than convert to strings.
// Also how the metadata is structured and applied here is insane.
std::vector<std::string> getTagsForResidue(const cds::ResidueAttributes& residueAttributes)
{
    std::string message = "Searching for attributes for " + residueAttributes.name;
    gmml::log(__LINE__, __FILE__, gmml::INF, message);
    std::vector<std::string> foundAttributes;
    std::vector<std::string> nCarbonSix        = {"Tal", "All", "Alt", "Fuc", "Gal", "Glc", "Gul", "Man",
                                                  "Qui", "Rha", "Ido", "Fru", "Sor", "Tag", "Psi"};
    std::vector<std::string> glucoGauche       = {"Glc", "All", "Alt", "Man"};
    std::vector<std::string> aldoseResidues    = {"All", "Alt", "Ara", "Fuc", "Gal", "Glc", "Gul", "Ido",
                                                  "Lyx", "Man", "Qui", "Rha", "Rib", "Tal", "Xyl", "Tyv",
                                                  "dUA", "Bac", "Abe", "Oli", "AAT", "Mur", "man"}; // Note lowercase man
                                                                                                    // as in LDmanHep :(
    std::vector<std::string> ketoseResidues    = {"Fru", "Psi", "Sor", "Tag", "Neu", "KDN",
                                                  "KDO", "K3O", "Aci", "Fus", "Leg", "Pse"};
    std::vector<std::string> ulosonateResidues = {"Neu", "KDN", "KDO", "K3O", "Aci", "Fus", "Leg", "Pse"};
    if (codeUtils::contains(glucoGauche, residueAttributes.name))
    {
        foundAttributes.push_back("gauche-effect=gluco");
    }
    if (codeUtils::contains(nCarbonSix, residueAttributes.name))
    {
        foundAttributes.push_back("n-carbon=6");
    }
    if (codeUtils::contains(aldoseResidues, residueAttributes.name))
    {
        foundAttributes.push_back("aldose");
    }
    if (codeUtils::contains(ketoseResidues, residueAttributes.name))
    {
        foundAttributes.push_back("ketose");
    }
    if (codeUtils::contains(ulosonateResidues, residueAttributes.name))
    {
        foundAttributes.push_back("ulosonate");
    }
    if (residueAttributes.configuration == "a")
    {
        foundAttributes.push_back("alpha");
    };
    if (residueAttributes.configuration == "b")
    {
        foundAttributes.push_back("beta");
    };
    if (residueAttributes.ringType == "p")
    {
        foundAttributes.push_back("pyranose");
    };
    if (residueAttributes.ringType == "f")
    {
        foundAttributes.push_back("furanose");
    };
    foundAttributes.push_back(residueAttributes.isInternal ? "internal" : "external");
    foundAttributes.push_back(residueTypeToString(residueAttributes.type));
    return foundAttributes;
}

// Pass in the two atoms on either side the residue-residue linkage
std::vector<std::vector<size_t>> GlycamMetadata::getDihedralAngleDataEntriesForLinkage(
    const std::string& atom1Name, const cds::ResidueAttributes& residue1Attributes, const std::string& atom2Name,
    const cds::ResidueAttributes& residue2Attributes)
{
    const DihedralAngleDataTable& table = dihedralAngleDataTable();
    std::vector<size_t> matching_entries;
    std::vector<std::string> residue1_types = getTagsForResidue(residue1Attributes);
    std::vector<std::string> residue2_types = getTagsForResidue(residue2Attributes);
    // Go through each entry in the metadata
    for (size_t n = 0; n < table.entries.size(); n++)
    {
        const DihedralAngleData& entry = table.entries[n];
        // Create a regex of each entry's linking_atom1_ and 2_. These are regex queries.
        std::regex regex1(entry.linking_atom1_, std::regex_constants::ECMAScript);
        std::regex regex2(entry.linking_atom2_, std::regex_constants::ECMAScript);
        // If metadata entry matches (regex query) to the two linking atom names
        if ((std::regex_match(atom1Name, regex1)) && (std::regex_match(atom2Name, regex2)))
        {
            // Some entries have conditions for the residue, that they have certain tags. Make sure any conditions are
            // met:
            if (checkIfResidueConditionsAreSatisfied(residue1_types, entry.residue1_conditions_) &&
                checkIfResidueConditionsAreSatisfied(residue2_types, entry.residue2_conditions_))
            {

                // Always add a later entry, but remove earlier match if number_of_bonds_from_anomeric_carbon_ AND index
                // number are the same. I've overloaded the == and != operators in the DihedralAngleData struct to
                // evaluate those.
                auto areEquivalent = [&](size_t k)
                {
                    const DihedralAngleData& a = table.entries[n];
                    const DihedralAngleData& b = table.entries[k];
                    return a.index_ == b.index_ &&
                           a.number_of_bonds_from_anomeric_carbon_ == b.number_of_bonds_from_anomeric_carbon_;
                };
                matching_entries.erase(std::remove_if(matching_entries.begin(), matching_entries.end(), areEquivalent),
                                       matching_entries.end());
                matching_entries.push_back(n);
            }
        }
    }
    unsigned int maxMetadataDihedral = 0;
    for (auto& entry : matching_entries)
    {
        maxMetadataDihedral = std::max(maxMetadataDihedral, table.entries[entry].number_of_bonds_from_anomeric_carbon_);
    }
    std::vector<std::vector<size_t>> orderedEntries;
    orderedEntries.resize(maxMetadataDihedral);
    for (size_t n = 0; n < maxMetadataDihedral; n++)
    {
        std::copy_if(matching_entries.begin(), matching_entries.end(), std::back_inserter(orderedEntries[n]),
                     [&](size_t entry)
                     {
                         return table.entries[entry].number_of_bonds_from_anomeric_carbon_ - 1 == n;
                     });
    }
    return orderedEntries;
}

std::vector<size_t> GlycamMetadata::likelyMetadata(const DihedralAngleDataTable& table,
                                                   const std::vector<size_t>& entries)
{
    std::vector<size_t> returningMetadata;
    returningMetadata.reserve(entries.size());
    for (auto& entry : entries)
    {
        if (table.entries[entry].weight_ >= 0.01) // HARDCODE EVERYTHING.
        {
            returningMetadata.push_back(entry);
        }
    }
    return returningMetadata;
}
