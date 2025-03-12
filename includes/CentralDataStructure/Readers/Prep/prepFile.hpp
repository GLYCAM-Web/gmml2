#ifndef INCLUDES_CENTRALDATASTRUCTURE_READERS_PREP_PREPFILE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_READERS_PREP_PREPFILE_HPP

#include "includes/CentralDataStructure/molecule.hpp"
#include "includes/CentralDataStructure/Readers/Prep/prepAtom.hpp"
#include "includes/CentralDataStructure/Readers/Prep/prepResidue.hpp"
#include <map>
#include <string>
#include <fstream>
#include <ostream>
#include <vector>

namespace prep
{
    class PrepFile : public cds::Molecule
    {
      public:
        //////////////////////////////////////////////////////////
        //                     TYPE DEFINITION                  //
        //////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////
        //                       Constructor                    //
        //////////////////////////////////////////////////////////
        PrepFile(const std::string& prep_file);
        PrepFile(const std::string& prep_file, const std::vector<std::string> queryNames);
        //////////////////////////////////////////////////////////
        //                           ACCESSOR                   //
        //////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////
        //                           MUTATOR                    //
        //////////////////////////////////////////////////////////
        void SetAtomConnectivities();
        void Generate3dStructures();
        //////////////////////////////////////////////////////////
        //                         FUNCTIONS                    //
        //////////////////////////////////////////////////////////
        void Write(const std::string& prep_file);
        void Write(std::ostream& out_stream);
        //////////////////////////////////////////////////////////
        //                     DISPLAY FUNCTIONS                //
        //////////////////////////////////////////////////////////
        std::string Print() const;

      private:
        //////////////////////////////////////////////////////////
        //                           ACCESSOR                   //
        //////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////
        //                           MUTATOR                    //
        //////////////////////////////////////////////////////////
        void ReadAllResidues(std::ifstream& in_file);
        void ReadQueryResidues(std::ifstream& in_file, const std::vector<std::string>& queryNames);
        // void ReadOnlyQueryResidues(std::ifstream &in_file, std::vector<std::string>& query_residue_names);
        //////////////////////////////////////////////////////////
        //                         ATTRIBUTES                   //
        //////////////////////////////////////////////////////////
    };
} // namespace prep
#endif
