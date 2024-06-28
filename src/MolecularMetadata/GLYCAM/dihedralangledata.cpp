#include "includes/MolecularMetadata/GLYCAM/dihedralangledata.hpp"
#include "includes/MolecularMetadata/GLYCAM/glycam06residueinfo.hpp"
#include "includes/CodeUtils/logging.hpp"
#include <regex>
//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
using gmml::MolecularMetadata::GLYCAM::DihedralAngleDataContainer;
using gmml::MolecularMetadata::GLYCAM::DihedralAngleDataVector;

//////////////////////////////////////////////////////////
//                      QUERY FUNCTIONS                 //
//////////////////////////////////////////////////////////
// Pass in the two atoms on either side the residue-residue linkage
std::vector<DihedralAngleDataVector>
DihedralAngleDataContainer::GetEntriesForLinkage(const std::string atom1Name, const std::string residue1Name,
                                                 const std::string atom2Name, const std::string residue2Name) const
{
    DihedralAngleDataVector matching_entries;
    Glycam06NamesToTypesLookupContainer metadata_residueNamesToTypes;
    std::vector<std::string> residue1_types = metadata_residueNamesToTypes.GetTypesForResidue(residue1Name);
    std::vector<std::string> residue2_types = metadata_residueNamesToTypes.GetTypesForResidue(residue2Name);
    // Go through each entry in the metadata
    for (const auto& entry : dihedralAngleDataVector_)
    {
        // Create a regex of each entry's linking_atom1_ and 2_. These are regex queries.
        // std::cout << "Compare entry " << entry.linking_atom1_ << "-" << entry.linking_atom2_ << " : " <<
        // linking_atom1->GetName() << "-" << linking_atom2->GetName() <<"\n";
        std::regex regex1(entry.linking_atom1_, std::regex_constants::ECMAScript);
        std::regex regex2(entry.linking_atom2_, std::regex_constants::ECMAScript);
        // If metadata entry matches (regex query) to the two linking atom names
        if ((std::regex_match(atom1Name, regex1)) && (std::regex_match(atom2Name, regex2)))
        {
            // Some entries have conditions for the residue, that they have certain tags. Make sure any conditions are
            // met:
            // std::cout << "Matched. Checking if conditions apply.\n";
            if ((checkIfResidueConditionsAreSatisfied(residue1_types, entry.residue1_conditions_)) &&
                (checkIfResidueConditionsAreSatisfied(residue2_types, entry.residue2_conditions_)))
            {
                // std::cout << "Found a match: " << entry.linking_atom1_ << "-" << entry.linking_atom2_ << ", " <<
                // entry.dihedral_angle_name_ << "\n"; Always add a later entry, but remove earlier match if
                // number_of_bonds_from_anomeric_carbon_ AND index number are the same.
                //  I've overloaded the == and != operators in the DihedralAngleData struct to evaluate those.
                //  This next line removes any elements of matching_entries that match "entry", then the line after adds
                //  entry.
                matching_entries.erase(std::remove(matching_entries.begin(), matching_entries.end(), entry),
                                       matching_entries.end());
                matching_entries.push_back(entry);
            }
        }
    }

    unsigned int maxMetadataDihedral = 0;
    for (auto& entry : matching_entries)
    {
        maxMetadataDihedral = std::max(maxMetadataDihedral, entry.number_of_bonds_from_anomeric_carbon_);
    }
    std::vector<DihedralAngleDataVector> orderedEntries;
    orderedEntries.resize(maxMetadataDihedral);
    for (size_t n = 0; n < maxMetadataDihedral; n++)
    {
        std::copy_if(matching_entries.begin(), matching_entries.end(), std::back_inserter(orderedEntries[n]),
                     [&](auto& entry)
                     {
                         return entry.number_of_bonds_from_anomeric_carbon_ - 1 == n;
                     });
    }
    return orderedEntries;
}

//////////////////////////////////////////////////////////
//                    PRIVATE FUNCTIONS                 //
//////////////////////////////////////////////////////////
// Some entries have conditions for the first or second residue to have a particular type (aka tag).
// Most entries have "none" for condition. This checks first if condition is "none", and therefore satisfied.
// Otherwise (else if) it checks if any of the residue_types match the condition for the entry, e.g.
// gauche_effect=galacto.
bool DihedralAngleDataContainer::checkIfResidueConditionsAreSatisfied(std::vector<std::string> residue_types,
                                                                      std::vector<std::string> entry_conditions) const
{
    for (const auto& entry_condition : entry_conditions)
    { // If no condition, return true. If can't find the condition in the list return false, otherwise, having found the
      // condition(s), return true.
        // gmml::log(__LINE__,__FILE__,gmml::INF, "Entry condition: " + entry_condition);
        if (entry_condition == "none")
        {
            // gmml::log(__LINE__,__FILE__,gmml::INF, "Returning true as conditions are none");
            return true;
        }
        if (!(std::find(residue_types.begin(), residue_types.end(), entry_condition) != residue_types.end()))
        {
            // gmml::log(__LINE__,__FILE__,gmml::INF, "Returning false as did not find the condition in residue tags");
            return false; // If any condition isn't satisified. return false.
        }
    }
    // gmml::log(__LINE__,__FILE__,gmml::INF, "All residue conditions are satisfied");
    return true;
}

//////////////////////////////////////////////////////////
//                    INITIALIZER                       //
//////////////////////////////////////////////////////////

// Struct is copied here for reference.
// struct DihedralAngleData
//{
//    std::string linking_atom1_ ;
//    std::string linking_atom2_ ;
//    std::string dihedral_angle_name_ ;
//    double default_angle_value_ ;
//    double lower_deviation_ ;
//    double upper_deviation_ ;
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
 * Chi angle index numbering varies depending on side chain length, so in ASN the Chi1 is 4 bonds away from the sugar,
 * so would be 4 If two entries have matching regex and have the same index_ number (e.g. 1), the first will be
 * overwritten. If two entries have matching regex and different index_ numbers (e.g. 1,2,3) they will all be used to
 * create multiple rotamers/conformers If two entries have different regex, but apply to the same dihedral angle (e.g.
 * Phi), give them the same index_ number (e.g. 1).
 */

// clang-format off
DihedralAngleDataContainer::DihedralAngleDataContainer()
{ // const AmberAtomTypeInfo Glycam06j1AtomTypes[] =
    dihedralAngleDataVector_ =
      { // Regex1  , Regex2   , Name   , Angle  , Upper  , Lower  , Weight, Entry Type    , Name , B , I , Res1 Condition , Res2 Conditions           , Atom names                                                               // Atom names this applies to
          { "C1"   , "O[1-9]" , "Phi"  , 180.0  ,  25.0  ,  25.0  , 1.0   , "permutation" , "g"  , 1 , 1 , {"aldose"}     , {"monosaccharide"}                  , "C2" , "C1" , "O." , "C."  }, // Phi should be C2-C1(ano)-Ox-Cx, or C1-C2(ano)-Ox-Cx
          { "C2"   , "O[1-9]" , "Phi"  , -60.0  ,  25.0  ,  25.0  , 1.0   , "permutation" , "g"  , 1 , 1 , {"n-carbon=6", "ketose", "alpha"}     , {"monosaccharide"}            , "C1" , "C2" , "O." , "C."  }, // Phi is defined by C1-C2(ano)-Ox-Cx for ketoses like Fru
          { "C2"   , "O[1-9]" , "Phi"  ,  60.0  ,  25.0  ,  25.0  , 1.0   , "permutation" , "-g" , 1 , 1 , {"n-carbon=6", "ketose", "beta"}     , {"monosaccharide"}            , "C1" , "C2" , "O." , "C."  }, // Phi is defined by C1-C2(ano)-Ox-Cx for ketoses like Fru
          { "C2"   , "O[1-9]" , "Phi"  , 180.0  ,  25.0  ,  25.0  , 1.0   , "permutation" , "t"  , 1 , 1 , {"ketose", "ulosonate"}     , {"monosaccharide"}      , "C1" , "C2" , "O." , "C."  }, // Phi should be C2-C1(ano)-Ox-Cx, or C1-C2(ano)-Ox-Cx
          // Sialic acid type 2- linkages:
          { "C2"   , "O[3-6]" , "Phi"  , -60.0  ,  25.0  ,  25.0  , 1.0   , "permutation" , "-g" , 1 , 2 , {"ulosonate", "alpha"}  , {"monosaccharide"}         , "C1" , "C2" , "O." , "C."  },
          // Generic psi linkages, why is this labelled "ap"?
          { "C."   , "O[1-5]" , "Psi"  ,   0.0  ,  40.0  ,  40.0  , 1.0   , "permutation" , "ap" , 2 , 1 , {"monosaccharide"}       , {"monosaccharide"}                  , "C." , "O." , "C." , "H."  }, // Psi should be C(ano)-Ox-Cx-Hx, if Cx is ring, otherwise, C(ano)-Ox-Cx-C(x-1)
          { "C."   , "O[6-9]" , "Psi"  , 180.0  ,  40.0  ,  40.0  , 1.0   , "permutation" , "t"  , 2 , 1 , {"monosaccharide"}       , {"monosaccharide"}                  , "C." , "O." , "C." , "C."  },
          // Omega angle in x-6 linkages.
          { "C.*"  , "O6"     , "Omg"  , -60.0  ,  20.0  ,  20.0  , 1.0   , "permutation" , "gg" , 3 , 1 , {"none"}       , {"monosaccharide"}                  , "O6" , "C6" , "C5" , "O5"  }, // omg is O6-C5-C5-O5
          { "C.*"  , "O6"     , "Omg"  ,  60.0  ,  20.0  ,  20.0  , 1.0   , "permutation" , "gt" , 3 , 2 , {"none"}       , {"monosaccharide"}                  , "O6" , "C6" , "C5" , "O5"  },
          { "C.*"  , "O6"     , "Omg"  , 180.0  ,  20.0  ,  20.0  , 1.0   , "permutation" , "tg" , 3 , 3 , {"none"}       , {"gauche-effect=galacto"} , "O6" , "C6" , "C5" , "O5"  },
          { "C.*"  , "O6"     , "Omg"  , 180.0  ,  20.0  ,  20.0  , 0.001 , "permutation" , "tg" , 3 , 3 , {"none"}       , {"gauche-effect=gluco"}   , "O6" , "C6" , "C5" , "O5"  },
          // Omega angle in x-5 linkages.
          { "C.*"  , "O5"     , "Omg"  , -60.0  ,  20.0  ,  20.0  , 1.0   , "permutation" , "gg" , 3 , 1 , {"none"}       , {"furanose"}                  , "O5" , "C5" , "C4" , "O4"  },
          { "C.*"  , "O5"     , "Omg"  ,  60.0  ,  20.0  ,  20.0  , 1.0   , "permutation" , "gt" , 3 , 2 , {"none"}       , {"furanose"}                  , "O5" , "C5" , "C4" , "O4"  },
          // Omega angle in x-1 linkages in ketose, furanoses
          { "C.*"  , "O1"     , "Omg"  , -60.0  ,  20.0  ,  20.0  , 1.0   , "permutation" , "gg" , 3 , 1 , {"none"}       , {"furanose", "ketose"}        , "O1" , "C1" , "C2" , "O6"  },
          { "C.*"  , "O1"     , "Omg"  ,  60.0  ,  20.0  ,  20.0  , 1.0   , "permutation" , "gt" , 3 , 2 , {"none"}       , {"furanose", "ketose"}        , "O1" , "C1" , "C2" , "O6"  },
          // 2-7 linkages copied from GlycamWeb Jan 2021. Branching in the linkage causes oddities with the bond number and which atoms are chosen for the torsion.
          { "C2"   , "O7"     , "Phi"  , -60.0  ,  20.0  ,  20.0  , 1.0   , "permutation" , "-g" , 1 , 2 , {"ulosonate", "alpha"}  , {"ulosonate"}    , "C3" , "C2" , "O." , "C."  },
          { "C2"   , "O7"     , "Psi"  ,   0.0  ,  20.0  ,  20.0  , 1.0   , "permutation" , "c"  , 2 , 1 , {"ulosonate", "alpha"}  , {"ulosonate"}    , "C." , "O." , "C." , "H."  },
          { "C.*"  , "O7"     , "Omg7" ,  60.0  ,  20.0  ,  20.0  , 1.0   , "permutation" , "g"  , 3 , 1 , {"none"}                , {"ulosonate"}    , "O7" , "C7" , "C6" , "O6"  },
          { "C.*"  , "O7"     , "Omg9" , 180.0  ,  20.0  ,  20.0  , 1.0   , "permutation" , "t"  , 4 , 1 , {"none"}                , {"ulosonate"}    , "O9" , "C9" , "C8" , "C7"  },
          { "C.*"  , "O7"     , "Omg8" ,  60.0  ,  20.0  ,  20.0  , 1.0   , "permutation" , "g"  , 5 , 1 , {"none"}                , {"ulosonate"}    , "C9" , "C8" , "C7" , "O7"  },
          // 2-9 linkages copied from GlycamWeb Jan 2021.
          { "C2"   , "O9"     , "Phi"  , -60.0  ,  20.0  ,  20.0  , 1.0   , "permutation" , "-g" , 1 , 2 , {"ulosonate", "alpha"}  , {"ulosonate"}    , "C3" , "C2" , "O." , "C."  },
          { "C.*"  , "O9"     , "Omg9" , 180.0  ,  20.0  ,  20.0  , 1.0   , "permutation" , "t"  , 3 , 1 , {"none"}                , {"ulosonate"}    , "O9" , "C9" , "C8" , "C7"  },
          { "C.*"  , "O9"     , "Omg8" , 180.0  ,  20.0  ,  20.0  , 1.0   , "permutation" , "t"  , 4 , 1 , {"none"}                , {"ulosonate"}    , "C9" , "C8" , "C7" , "C6"  },
          { "C.*"  , "O9"     , "Omg7" , -60.0  ,  20.0  ,  20.0  , 1.0   , "permutation" , "-g" , 5 , 1 , {"none"}                , {"ulosonate"}    , "O7" , "C7" , "C6" , "O6"  },
           // Generic x-8 linkages. Values copied from external conformer A below. Arbitrary, but I have no data for them.
          { "C.*"   , "O8"     , "Phi"  , -79.5  ,  20.0  ,  20.0  , 0.42  , "conformer"   , "A"  , 1 , 1 , {"monosaccharide"}  , {"ulosonate"}    , "C1" , "C2" , "O8" , "C8"  },
          { "C.*"   , "O8"     , "Psi"  ,  88.1  ,  20.0  ,  20.0  , 0.42  , "conformer"   , "A"  , 2 , 1 , {"monosaccharide"}  , {"ulosonate"}    , "C2" , "O8" , "C8" , "C7"  },
          { "C.*"   , "O8"     , "Omg8" ,-170.1  ,  20.0  ,  20.0  , 0.42  , "conformer"   , "A"  , 3 , 1 , {"monosaccharide"}  , {"ulosonate"}    , "O8" , "C8" , "C7" , "C6"  },
          { "C.*"   , "O8"     , "Omg7" , -61.8  ,  20.0  ,  20.0  , 0.42  , "conformer"   , "A"  , 4 , 1 , {"monosaccharide"}  , {"ulosonate"}    , "C8" , "C7" , "C6" , "O6"  },
          { "C.*"   , "O8"     , "Omg9" ,  65.8  ,  20.0  ,  20.0  , 0.42  , "conformer"   , "A"  , 5 , 1 , {"monosaccharide"}  , {"ulosonate"}    , "O9" , "C9" , "C8" , "O8"  },
           // Internal 2-8 linkages
          { "C2"   , "O8"     , "Phi"  , -79.5  ,  20.0  ,  20.0  , 0.42  , "conformer"   , "A"  , 1 , 1 , {"ulosonate", "alpha", "internal"}  , {"ulosonate"}    , "C1" , "C2" , "O8" , "C8"  },
          { "C2"   , "O8"     , "Psi"  ,  88.1  ,  20.0  ,  20.0  , 0.42  , "conformer"   , "A"  , 2 , 1 , {"ulosonate", "alpha", "internal"}  , {"ulosonate"}    , "C2" , "O8" , "C8" , "C7"  },
          { "C2"   , "O8"     , "Omg8" ,-170.1  ,  20.0  ,  20.0  , 0.42  , "conformer"   , "A"  , 3 , 1 , {"ulosonate", "alpha", "internal"}  , {"ulosonate"}    , "O8" , "C8" , "C7" , "C6"  },
          { "C2"   , "O8"     , "Omg7" , -61.8  ,  20.0  ,  20.0  , 0.42  , "conformer"   , "A"  , 4 , 1 , {"ulosonate", "alpha", "internal"}  , {"ulosonate"}    , "C8" , "C7" , "C6" , "O6"  },
          { "C2"   , "O8"     , "Omg9" ,  65.8  ,  20.0  ,  20.0  , 0.42  , "conformer"   , "A"  , 5 , 1 , {"ulosonate", "alpha", "internal"}  , {"ulosonate"}    , "O9" , "C9" , "C8" , "O8"  },
          { "C2"   , "O8"     , "Phi"  , -73.6  ,  20.0  ,  20.0  , 0.24  , "conformer"   , "B"  , 1 , 2 , {"ulosonate", "alpha", "internal"}  , {"ulosonate"}    , "C1" , "C2" , "O8" , "C8"  },
          { "C2"   , "O8"     , "Psi"  , 109.2  ,  20.0  ,  20.0  , 0.24  , "conformer"   , "B"  , 2 , 2 , {"ulosonate", "alpha", "internal"}  , {"ulosonate"}    , "C2" , "O8" , "C8" , "C7"  },
          { "C2"   , "O8"     , "Omg8" ,-169.2  ,  20.0  ,  20.0  , 0.24  , "conformer"   , "B"  , 3 , 2 , {"ulosonate", "alpha", "internal"}  , {"ulosonate"}    , "O8" , "C8" , "C7" , "C6"  },
          { "C2"   , "O8"     , "Omg7" , -61.9  ,  20.0  ,  20.0  , 0.24  , "conformer"   , "B"  , 4 , 2 , {"ulosonate", "alpha", "internal"}  , {"ulosonate"}    , "C8" , "C7" , "C6" , "O6"  },
          { "C2"   , "O8"     , "Omg9" , -61.0  ,  20.0  ,  20.0  , 0.24  , "conformer"   , "B"  , 5 , 2 , {"ulosonate", "alpha", "internal"}  , {"ulosonate"}    , "O9" , "C9" , "C8" , "O8"  },
          { "C2"   , "O8"     , "Phi"  ,-172.0  ,  20.0  ,  20.0  , 0.08  , "conformer"   , "C"  , 1 , 3 , {"ulosonate", "alpha", "internal"}  , {"ulosonate"}    , "C1" , "C2" , "O8" , "C8"  },
          { "C2"   , "O8"     , "Psi"  , 108.7  ,  20.0  ,  20.0  , 0.08  , "conformer"   , "C"  , 2 , 3 , {"ulosonate", "alpha", "internal"}  , {"ulosonate"}    , "C2" , "O8" , "C8" , "C7"  },
          { "C2"   , "O8"     , "Omg8" ,-164.2  ,  20.0  ,  20.0  , 0.08  , "conformer"   , "C"  , 3 , 3 , {"ulosonate", "alpha", "internal"}  , {"ulosonate"}    , "O8" , "C8" , "C7" , "C6"  },
          { "C2"   , "O8"     , "Omg7" , -52.9  ,  20.0  ,  20.0  , 0.08  , "conformer"   , "C"  , 4 , 3 , {"ulosonate", "alpha", "internal"}  , {"ulosonate"}    , "C8" , "C7" , "C6" , "O6"  },
          { "C2"   , "O8"     , "Omg9" ,  68.5  ,  20.0  ,  20.0  , 0.08  , "conformer"   , "C"  , 5 , 3 , {"ulosonate", "alpha", "internal"}  , {"ulosonate"}    , "O9" , "C9" , "C8" , "O8"  },
          { "C2"   , "O8"     , "Phi"  , -72.1  ,  20.0  ,  20.0  , 0.07  , "conformer"   , "D"  , 1 , 4 , {"ulosonate", "alpha", "internal"}  , {"ulosonate"}    , "C1" , "C2" , "O8" , "C8"  },
          { "C2"   , "O8"     , "Psi"  , 127.1  ,  20.0  ,  20.0  , 0.07  , "conformer"   , "D"  , 2 , 4 , {"ulosonate", "alpha", "internal"}  , {"ulosonate"}    , "C2" , "O8" , "C8" , "C7"  },
          { "C2"   , "O8"     , "Omg8" ,  58.5  ,  20.0  ,  20.0  , 0.07  , "conformer"   , "D"  , 3 , 4 , {"ulosonate", "alpha", "internal"}  , {"ulosonate"}    , "O8" , "C8" , "C7" , "C6"  },
          { "C2"   , "O8"     , "Omg7" , -66.1  ,  20.0  ,  20.0  , 0.07  , "conformer"   , "D"  , 4 , 4 , {"ulosonate", "alpha", "internal"}  , {"ulosonate"}    , "C8" , "C7" , "C6" , "O6"  },
          { "C2"   , "O8"     , "Omg9" , 178.9  ,  20.0  ,  20.0  , 0.07  , "conformer"   , "D"  , 5 , 4 , {"ulosonate", "alpha", "internal"}  , {"ulosonate"}    , "O9" , "C9" , "C8" , "O8"  },
          { "C2"   , "O8"     , "Phi"  , -42.9  ,  20.0  ,  20.0  , 0.07  , "conformer"   , "E"  , 1 , 5 , {"ulosonate", "alpha", "internal"}  , {"ulosonate"}    , "C1" , "C2" , "O8" , "C8"  },
          { "C2"   , "O8"     , "Psi"  , 162.3  ,  20.0  ,  20.0  , 0.07  , "conformer"   , "E"  , 2 , 5 , {"ulosonate", "alpha", "internal"}  , {"ulosonate"}    , "C2" , "O8" , "C8" , "C7"  },
          { "C2"   , "O8"     , "Omg8" , -58.3  ,  20.0  ,  20.0  , 0.07  , "conformer"   , "E"  , 3 , 5 , {"ulosonate", "alpha", "internal"}  , {"ulosonate"}    , "O8" , "C8" , "C7" , "C6"  },
          { "C2"   , "O8"     , "Omg7" , -55.7  ,  20.0  ,  20.0  , 0.07  , "conformer"   , "E"  , 4 , 5 , {"ulosonate", "alpha", "internal"}  , {"ulosonate"}    , "C8" , "C7" , "C6" , "O6"  },
          { "C2"   , "O8"     , "Omg9" , -58.5  ,  20.0  ,  20.0  , 0.07  , "conformer"   , "E"  , 5 , 5 , {"ulosonate", "alpha", "internal"}  , {"ulosonate"}    , "O9" , "C9" , "C8" , "O8"  },
          // External 2-8 linkages
          { "C2"   , "O8"     , "Phi"  , -66.0  ,  20.0  ,  20.0  , 0.42  , "conformer"   , "A"  , 1 , 1 , {"ulosonate", "alpha", "external"}  , {"ulosonate"}    , "C1" , "C2" , "O8" , "C8"  },
          { "C2"   , "O8"     , "Psi"  , 118.2  ,  20.0  ,  20.0  , 0.42  , "conformer"   , "A"  , 2 , 1 , {"ulosonate", "alpha", "external"}  , {"ulosonate"}    , "C2" , "O8" , "C8" , "C7"  },
          { "C2"   , "O8"     , "Omg8" ,-160.4  ,  20.0  ,  20.0  , 0.42  , "conformer"   , "A"  , 3 , 1 , {"ulosonate", "alpha", "external"}  , {"ulosonate"}    , "O8" , "C8" , "C7" , "C6"  },
          { "C2"   , "O8"     , "Omg7" , -59.6  ,  20.0  ,  20.0  , 0.42  , "conformer"   , "A"  , 4 , 1 , {"ulosonate", "alpha", "external"}  , {"ulosonate"}    , "C8" , "C7" , "C6" , "O6"  },
          { "C2"   , "O8"     , "Omg9" ,  73.4  ,  20.0  ,  20.0  , 0.42  , "conformer"   , "A"  , 5 , 1 , {"ulosonate", "alpha", "external"}  , {"ulosonate"}    , "O9" , "C9" , "C8" , "O8"  },
          { "C2"   , "O8"     , "Phi"  , -68.4  ,  20.0  ,  20.0  , 0.24  , "conformer"   , "B"  , 1 , 2 , {"ulosonate", "alpha", "external"}  , {"ulosonate"}    , "C1" , "C2" , "O8" , "C8"  },
          { "C2"   , "O8"     , "Psi"  , 132.7  ,  20.0  ,  20.0  , 0.24  , "conformer"   , "B"  , 2 , 2 , {"ulosonate", "alpha", "external"}  , {"ulosonate"}    , "C2" , "O8" , "C8" , "C7"  },
          { "C2"   , "O8"     , "Omg8" ,  64.4  ,  20.0  ,  20.0  , 0.24  , "conformer"   , "B"  , 3 , 2 , {"ulosonate", "alpha", "external"}  , {"ulosonate"}    , "O8" , "C8" , "C7" , "C6"  },
          { "C2"   , "O8"     , "Omg7" , -61.1  ,  20.0  ,  20.0  , 0.24  , "conformer"   , "B"  , 4 , 2 , {"ulosonate", "alpha", "external"}  , {"ulosonate"}    , "C8" , "C7" , "C6" , "O6"  },
          { "C2"   , "O8"     , "Omg9" ,-179.0  ,  20.0  ,  20.0  , 0.24  , "conformer"   , "B"  , 5 , 2 , {"ulosonate", "alpha", "external"}  , {"ulosonate"}    , "O9" , "C9" , "C8" , "O8"  },
          { "C2"   , "O8"     , "Phi"  , -53.2  ,  20.0  ,  20.0  , 0.08  , "conformer"   , "C"  , 1 , 3 , {"ulosonate", "alpha", "external"}  , {"ulosonate"}    , "C1" , "C2" , "O8" , "C8"  },
          { "C2"   , "O8"     , "Psi"  , 161.6  ,  20.0  ,  20.0  , 0.08  , "conformer"   , "C"  , 2 , 3 , {"ulosonate", "alpha", "external"}  , {"ulosonate"}    , "C2" , "O8" , "C8" , "C7"  },
          { "C2"   , "O8"     , "Omg8" , -67.8  ,  20.0  ,  20.0  , 0.08  , "conformer"   , "C"  , 3 , 3 , {"ulosonate", "alpha", "external"}  , {"ulosonate"}    , "O8" , "C8" , "C7" , "C6"  },
          { "C2"   , "O8"     , "Omg7" , -61.1  ,  20.0  ,  20.0  , 0.08  , "conformer"   , "C"  , 4 , 3 , {"ulosonate", "alpha", "external"}  , {"ulosonate"}    , "C8" , "C7" , "C6" , "O6"  },
          { "C2"   , "O8"     , "Omg9" , -57.9  ,  20.0  ,  20.0  , 0.08  , "conformer"   , "C"  , 5 , 3 , {"ulosonate", "alpha", "external"}  , {"ulosonate"}    , "O9" , "C9" , "C8" , "O8"  },
          { "C2"   , "O8"     , "Phi"  , -72.7  ,  20.0  ,  20.0  , 0.07  , "conformer"   , "D"  , 1 , 4 , {"ulosonate", "alpha", "external"}  , {"ulosonate"}    , "C1" , "C2" , "O8" , "C8"  },
          { "C2"   , "O8"     , "Psi"  , 124.9  ,  20.0  ,  20.0  , 0.07  , "conformer"   , "D"  , 2 , 4 , {"ulosonate", "alpha", "external"}  , {"ulosonate"}    , "C2" , "O8" , "C8" , "C7"  },
          { "C2"   , "O8"     , "Omg8" ,  62.4  ,  20.0  ,  20.0  , 0.07  , "conformer"   , "D"  , 3 , 4 , {"ulosonate", "alpha", "external"}  , {"ulosonate"}    , "O8" , "C8" , "C7" , "C6"  },
          { "C2"   , "O8"     , "Omg7" , -59.5  ,  20.0  ,  20.0  , 0.07  , "conformer"   , "D"  , 4 , 4 , {"ulosonate", "alpha", "external"}  , {"ulosonate"}    , "C8" , "C7" , "C6" , "O6"  },
          { "C2"   , "O8"     , "Omg9" ,  78.8  ,  20.0  ,  20.0  , 0.07  , "conformer"   , "D"  , 5 , 4 , {"ulosonate", "alpha", "external"}  , {"ulosonate"}    , "O9" , "C9" , "C8" , "O8"  },
          { "C2"   , "O8"     , "Phi"  ,-164.9  ,  20.0  ,  20.0  , 0.07  , "conformer"   , "E"  , 1 , 5 , {"ulosonate", "alpha", "external"}  , {"ulosonate"}    , "C1" , "C2" , "O8" , "C8"  },
          { "C2"   , "O8"     , "Psi"  , 103.0  ,  20.0  ,  20.0  , 0.07  , "conformer"   , "E"  , 2 , 5 , {"ulosonate", "alpha", "external"}  , {"ulosonate"}    , "C2" , "O8" , "C8" , "C7"  },
          { "C2"   , "O8"     , "Omg8" ,-161.3  ,  20.0  ,  20.0  , 0.07  , "conformer"   , "E"  , 3 , 5 , {"ulosonate", "alpha", "external"}  , {"ulosonate"}    , "O8" , "C8" , "C7" , "C6"  },
          { "C2"   , "O8"     , "Omg7" , -54.3  ,  20.0  ,  20.0  , 0.07  , "conformer"   , "E"  , 4 , 5 , {"ulosonate", "alpha", "external"}  , {"ulosonate"}    , "C8" , "C7" , "C6" , "O6"  },
          { "C2"   , "O8"     , "Omg9" ,  71.1  ,  20.0  ,  20.0  , 0.07  , "conformer"   , "E"  , 5 , 5 , {"ulosonate", "alpha", "external"}  , {"ulosonate"}    , "O9" , "C9" , "C8" , "O8"  },
          // Sucrose from: Bock K, Lemieux RU. Carbohydrate Research, 100 (1982) 63-74. Now allowing ordering of sequence either way, thus four entries instead of two.
          { "C1"   , "O2"     , "Phi"  ,-135.0  ,   5.0  ,   5.0  , 1.0   , "permutation" , ""   , 1 , 1 , {"pyranose", "aldose"}    , {"furanose", "ketose"}         , "C2" , "C1" , "O2" , "C2" },
          { "C1"   , "O2"     , "Psi"  ,  80.0  ,  10.0  ,  10.0  , 1.0   , "permutation" , ""   , 2 , 1 , {"pyranose", "aldose"}    , {"furanose", "ketose"}         , "C1" , "O2" , "C2" , "C1" },
          { "C2"   , "O1"     , "Phi"  ,  80.0  ,  10.0  ,  10.0  , 1.0   , "permutation" , ""   , 1 , 1 , {"furanose", "ketose"}    , {"pyranose", "aldose"}         , "C1" , "C2" , "O1" , "C1" },
          { "C2"   , "O1"     , "Psi"  ,-135.0  ,   5.0  ,   5.0  , 1.0   , "permutation" , ""   , 2 , 1 , {"furanose", "ketose"}    , {"pyranose", "aldose"}         , "C2" , "O1" , "C1" , "C2" },
          // Ketose 1-1 linkages. Copied from old GLYCAM-Web builder.
          { "C1"   , "O1"     , "Phi"  , 180.0  ,  20.0  ,  20.0  , 1.0   , "permutation" , ""   , 1 , 1 , {"none"}       , {"ketose"}                , "C2" , "C1" , "O1" , "C1" },
          { "C1"   , "O1"     , "Psi"  ,-150.0  ,  20.0  ,  20.0  , 1.0   , "permutation" , ""   , 2 , 1 , {"none"}       , {"ketose"}                , "C1" , "O1" , "C1" , "C2" },
          { "C1"   , "O1"     , "Omg"  ,  60.0  ,  20.0  ,  20.0  , 1.0   , "permutation" , ""   , 3 , 1 , {"none"}       , {"ketose"}                , "O1" , "C1" , "C2" , "O5" },
        // Common sugar derivatives
        // Phosphate/sulfate
          { "[SP]1", "N[2]"   , "Phi"  , 180.0  ,  20.0  ,  20.0  , 1.0   , "permutation" , ""   , 1 , 1 , {"derivative"}       , {"monosaccharide"}                  , "O." , ".1" , "N." , "C."  },
          { "[SP]1", "N[2]"   , "Psi"  , -40.0  ,  20.0  ,  20.0  , 1.0   , "permutation" , ""   , 2 , 1 , {"derivative"}       , {"monosaccharide"}                  , ".1" , "N." , "C." , "H."  },
          { "[SP]1", "O[1-9]" , "Phi"  , 180.0  ,  20.0  ,  20.0  , 1.0   , "permutation" , ""   , 1 , 1 , {"derivative"}       , {"monosaccharide"}                  , "O." , ".1" , "O." , "C."  },
          { "[SP]1", "O[1-5]" , "Psi"  ,   0.0  ,  20.0  ,  20.0  , 1.0   , "permutation" , ""   , 2 , 1 , {"derivative"}       , {"monosaccharide"}                  , ".1" , "O." , "C." , "H."  },
          { "[SP]1", "O[6-9]" , "Psi"  , 180.0  ,  20.0  ,  20.0  , 1.0   , "permutation" , ""   , 2 , 1 , {"derivative"}       , {"monosaccharide"}                  , ".1" , "O." , "C." , "C."  },
          { "[SP]1", "O6"     , "Omg"  , -60.0  ,  20.0  ,  20.0  , 1.0   , "permutation" , "gg" , 3 , 1 , {"derivative"}       , {"monosaccharide"}                  , "O6" , "C6" , "C5" , "O5"  },
  //        These two were removed to stop them showing up as options the carb builder, but for some applications like grafting you may want them enabled. No functionality for that context yet tho.
  //        { "[SP]1", "O6"     , "Omg"  ,  60.0  ,  20.0  ,  20.0  , 1.0   , "permutation" , "gt" , 3 , 2 , {"derivative"}       , {"monosaccharide"}                  , "O6" , "C6" , "C5" , "O5"  },
  //        { "[SP]1", "O6"     , "Omg"  , 180.0  ,  20.0  ,  20.0  , 1.0   , "permutation" , "tg" , 3 , 3 , {"derivative"}       , {"gauche-effect=galacto"} , "O6" , "C6" , "C5" , "O5"  },
        // Ac ACX
          { "C1A"  , "O[1-9]" , "Phi"  , 180.0  ,  20.0  ,  20.0  , 1.0   , "permutation" , "t"  , 1 , 1 , {"derivative"}       , {"monosaccharide"}                  , "C2A", "C1A", "O." , "C."  },
          { "C1A"  , "O[1-5]" , "Psi"  ,   0.0  ,  20.0  ,  20.0  , 1.0   , "permutation" , "c"  , 2 , 1 , {"derivative"}       , {"monosaccharide"}                  , "C1A", "O." , "C." , "H."  },
          { "C1A"  , "O[6-9]" , "Psi"  , 180.0  ,  20.0  ,  20.0  , 1.0   , "permutation" , "-g" , 2 , 1 , {"derivative"}       , {"monosaccharide"}                  , "C1A", "O." , "C." , "C."  },
        // Me MEX
          { "CH3"  , "O[1-9]" , "Phi"  , 180.0  ,  20.0  ,  20.0  , 1.0   , "permutation" , "t"  , 1 , 1 , {"derivative"}       , {"monosaccharide"}                  , "CH3", "O." , "C." , "H."  },
          { "CH3"  , "O[1-5]" , "Psi"  ,   0.0  ,  20.0  ,  20.0  , 1.0   , "permutation" , "c"  , 2 , 1 , {"derivative"}       , {"monosaccharide"}                  , "CH3", "O." , "C." , "H."  },
          { "CH3"  , "O[6-9]" , "Psi"  , 180.0  ,  20.0  ,  20.0  , 1.0   , "permutation" , "-g" , 2 , 1 , {"derivative"}       , {"monosaccharide"}                  , "CH3", "O." , "C." , "C."  },
        // Any derivative linked to O8 of Sia. Copied from 2-8 linkage values.
          { ".*"   , "O8"     , "Omg8" ,-160.4  ,  20.0  ,  20.0  , 1.0   , "permutation" , ""   , 3 , 1 , {"derivative"}                        , {"ulosonate"}    , "O8" , "C8" , "C7" , "C6"  },
          { ".*"   , "O8"     , "Omg7" , -59.6  ,  20.0  ,  20.0  , 1.0   , "permutation" , ""   , 4 , 1 , {"derivative"}                        , {"ulosonate"}    , "C8" , "C7" , "C6" , "O6"  },
          { ".*"   , "O8"     , "Omg9" ,  73.4  ,  20.0  ,  20.0  , 1.0   , "permutation" , ""   , 5 , 1 , {"derivative"}                        , {"ulosonate"}    , "O9" , "C9" , "C8" , "O8"  },
          // Protein linkages
          // ASN // Values are from Petrescu et al 2004.
          { "C."   , "ND2"    , "Chi1" , 191.6  ,  14.4  ,  14.4  , 0.497 , "conformer"   , "A"  , 4 , 1 , {"none"}       , {"amino-acid"}            , "CG" , "CB" , "CA" , "N"   },
          { "C."   , "ND2"    , "Chi2" , 177.6  ,  43.0  ,  43.0  , 0.497 , "conformer"   , "A"  , 3 , 1 , {"none"}       , {"amino-acid"}            , "ND2", "CG" , "CB" , "CA"  },
          { "C."   , "ND2"    , "Psi"  , 177.3  ,  12.3  ,  12.3  , 0.497 , "conformer"   , "A"  , 2 , 1 , {"none"}       , {"amino-acid"}            , "C." , "ND2", "CG" , "CB"  },
          { "C1"   , "ND2"    , "Phi"  , 261.0  ,  21.3  ,  21.3  , 0.497 , "conformer"   , "A"  , 1 , 1 , {"none"}       , {"amino-acid"}            , "C." , "C." , "ND2", "CG"  },
          { "C."   , "ND2"    , "Chi1" ,  63.6  ,   8.9  ,   8.9  , 0.178 , "conformer"   , "B"  , 4 , 2 , {"none"}       , {"amino-acid"}            , "CG" , "CB" , "CA" , "N"   },
          { "C."   , "ND2"    , "Chi2" , 191.1  ,  31.6  ,  31.6  , 0.178 , "conformer"   , "B"  , 3 , 2 , {"none"}       , {"amino-acid"}            , "ND2", "CG" , "CB" , "CA"  },
          { "C."   , "ND2"    , "Psi"  , 178.5  ,  13.9  ,  13.9  , 0.178 , "conformer"   , "B"  , 2 , 2 , {"none"}       , {"amino-acid"}            , "C." , "ND2", "CG" , "CB"  },
          { "C1"   , "ND2"    , "Phi"  , 253.7  ,  21.5  ,  21.5  , 0.178 , "conformer"   , "B"  , 1 , 2 , {"none"}       , {"amino-acid"}            , "C." , "C." , "ND2", "CG"  },
          { "C."   , "ND2"    , "Chi1" , 290.6  ,  12.7  ,  12.7  , 0.235 , "conformer"   , "C"  , 4 , 3 , {"none"}       , {"amino-acid"}            , "CG" , "CB" , "CA" , "N"   },
          { "C."   , "ND2"    , "Chi2" , 152.9  ,  23.9  ,  23.9  , 0.235 , "conformer"   , "C"  , 3 , 3 , {"none"}       , {"amino-acid"}            , "ND2", "CG" , "CB" , "CA"  },
          { "C."   , "ND2"    , "Psi"  , 173.1  ,  12.2  ,  12.2  , 0.235 , "conformer"   , "C"  , 2 , 3 , {"none"}       , {"amino-acid"}            , "C." , "ND2", "CG" , "CB"  },
          { "C1"   , "ND2"    , "Phi"  , 268.0  ,  20.3  ,  20.3  , 0.235 , "conformer"   , "C"  , 1 , 3 , {"none"}       , {"amino-acid"}            , "C." , "C." , "ND2", "CG"  },
          { "C."   , "ND2"    , "Chi1" , 302.3  ,  11.5  ,  11.5  , 0.090 , "conformer"   , "D"  , 4 , 4 , {"none"}       , {"amino-acid"}            , "CG" , "CB" , "CA" , "N"   },
          { "C."   , "ND2"    , "Chi2" , 255.0  ,  28.8  ,  28.8  , 0.090 , "conformer"   , "D"  , 3 , 4 , {"none"}       , {"amino-acid"}            , "ND2", "CG" , "CB" , "CA"  },
          { "C."   , "ND2"    , "Psi"  , 178.1  ,  11.5  ,  11.5  , 0.090 , "conformer"   , "D"  , 2 , 4 , {"none"}       , {"amino-acid"}            , "C." , "ND2", "CG" , "CB"  },
          { "C1"   , "ND2"    , "Phi"  , 267.5  ,  23.9  ,  23.9  , 0.090 , "conformer"   , "D"  , 1 , 4 , {"none"}       , {"amino-acid"}            , "C." , "C." , "ND2", "CG"  },
          // THR // Values are from Lovell et al "PENULTIMATE ROTAMER LIBRARY"
          { "C."   , "OG1"    , "Chi1" ,  59.0  ,  10.0  ,  10.0  , 0.49  , "permutation" , "g"  , 3 , 1 , {"none"}       , {"amino-acid"}            , "OG1", "CB" , "CA" , "N"   },
          { "C."   , "OG1"    , "Chi1" ,-171.0  ,   6.0  ,   6.0  , 0.07  , "permutation" , "t"  , 3 , 2 , {"none"}       , {"amino-acid"}            , "OG1", "CB" , "CA" , "N"   },
          { "C."   , "OG1"    , "Chi1" , -61.0  ,   7.0  ,   7.0  , 0.43  , "permutation" , "-g" , 3 , 3 , {"none"}       , {"amino-acid"}            , "OG1", "CB" , "CA" , "N"   },
          { "C."   , "OG1"    , "Psi"  , -60.0  ,  60.0  ,  60.0  , 1.000 , "permutation" , "-g" , 2 , 1 , {"none"}       , {"amino-acid"}            , "C." , "OG1", "CB" , "CA"  },
          { "C."   , "OG1"    , "Phi"  , 180.0  ,  20.0  ,  20.0  , 1.000 , "permutation" , "t"  , 1 , 1 , {"none"}       , {"amino-acid"}            , "C." , "C." , "OG1", "CB"  },
           // SER // Values not checked
          { "C."   , "OG"     , "Chi1" ,  64.0  ,  10.0  ,  10.0  , 0.48  , "permutation" , "g"  , 3 , 1 , {"none"}       , {"amino-acid"}            , "OG" , "CB" , "CA" , "N"   },
          { "C."   , "OG"     , "Chi1" , 178.0  ,  11.0  ,  11.0  , 0.22  , "permutation" , "t"  , 3 , 2 , {"none"}       , {"amino-acid"}            , "OG" , "CB" , "CA" , "N"   },
          { "C."   , "OG"     , "Chi1" , -65.0  ,   9.0  ,   9.0  , 0.29  , "permutation" , "-g" , 3 , 3 , {"none"}       , {"amino-acid"}            , "OG" , "CB" , "CA" , "N"   },
          { "C."   , "OG"     , "Psi"  , -60.0  ,  20.0  ,  20.0  , 1.000 , "permutation" , "-g" , 2 , 1 , {"none"}       , {"amino-acid"}            , "C." , "OG" , "CB" , "CA"  },
          { "C."   , "OG"     , "Phi"  , 180.0  ,  20.0  ,  20.0  , 1.000 , "permutation" , "t"  , 1 , 1 , {"none"}       , {"amino-acid"}            , "C." , "C." , "OG1", "CB"  },
           // TYR // Values not checked
          { "C."   , "OH"     , "Chi1" , -60.0  ,  20.0  ,  20.0  , 1.000 , "permutation" , "-g" , 7 , 1 , {"none"}       , {"amino-acid"}            , "CG" , "CB" , "CA" , "N"   },
          { "C."   , "OH"     , "Chi1" ,  60.0  ,  20.0  ,  20.0  , 1.000 , "permutation" , "g"  , 7 , 2 , {"none"}       , {"amino-acid"}            , "CG" , "CB" , "CA" , "N"   },
          { "C."   , "OH"     , "Chi1" , 180.0  ,  20.0  ,  20.0  , 1.000 , "permutation" , "t"  , 7 , 3 , {"none"}       , {"amino-acid"}            , "CG" , "CB" , "CA" , "N"   },
          { "C."   , "OH"     , "Chi2" , -60.0  ,  20.0  ,  20.0  , 1.000 , "permutation" , "-g" , 6 , 1 , {"none"}       , {"amino-acid"}            , "CD1", "CG" , "CB" , "CA"  },
          { "C."   , "OH"     , "Psi"  , -60.0  ,  20.0  ,  20.0  , 1.000 , "permutation" , "-g" , 2 , 1 , {"none"}       , {"amino-acid"}            , "C." , "OH" , "CZ" , "CE1" },
          { "C."   , "OH"     , "Phi"  , 180.0  ,  20.0  ,  20.0  , 1.000 , "permutation" , "t"  , 1 , 1 , {"none"}       , {"amino-acid"}            , "C." , "C." , "OH ", "CZ"  },
          // ROH
          { "C1"   , "O1"     , "Phi"  , 180.0  ,  20.0  ,  20.0  , 1.0   , "permutation" , "g"  , 1 , 1 , {"none"}     , {"aglycon"}               , "C2" , "C1" , "O1" , "H1"  },
          { "C2"   , "O1"     , "Phi"  , 180.0  ,  20.0  ,  20.0  , 1.0   , "permutation" , "g"  , 1 , 1 , {"none"}     , {"aglycon"}               , "C2" , "C1" , "O1" , "H1"  },
          { "C1"   , "O"      , "Phi"  , 180.0  ,  20.0  ,  20.0  , 1.0   , "permutation" , "g"  , 1 , 1 , {"none"}     , {"aglycon"}               , "C2" , "C1" , "O1" , "H1"  },
          { "C2"   , "O"      , "Phi"  , 180.0  ,  20.0  ,  20.0  , 1.0   , "permutation" , "g"  , 1 , 1 , {"none"}     , {"aglycon"}               , "C2" , "C1" , "O1" , "H1"  },
      };
}

// clang-format on
//    Statistical analysis of the protein environment of N-glycosylation sites: implications for occupancy, structure,
//    and folding Andrei-J. Petrescu  Adina-L. Milac  Stefana M. Petrescu  Raymond A. Dwek Mark R. Wormald Glycobiology,
//    Volume 14, Issue 2, 1 February 2004, Pages 103–114,
