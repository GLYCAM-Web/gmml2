#ifndef INCLUDES_CENTRALDATASTRUCTURE_READERS_PDB_PDBRESIDUE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_READERS_PDB_PDBRESIDUE_HPP

#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbAtom.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbResidueId.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbData.hpp"

#include <string>
#include <sstream>
#include <ostream>

namespace pdb
{
    class PdbResidue : public cds::Residue
    {
      public:
        //////////////////////////////////////////////////////////
        //                       CONSTRUCTOR                    //
        //////////////////////////////////////////////////////////
        PdbResidue(PdbData& data, size_t residueId, std::stringstream& singleResidueSecion, std::string firstLine);
        PdbResidue(PdbData&, size_t, const std::string residueName, const PdbResidue* referenceResidue);
        //////////////////////////////////////////////////////////
        //                       ACCESSOR                       //
        //////////////////////////////////////////////////////////
        ResidueId getId() const;

        inline const std::string& getInsertionCode() const
        {
            return insertionCode_;
        }

        inline const std::string& getChainId() const
        {
            return chainId_;
        }

        const std::string getNumberAndInsertionCode() const;

        inline bool HasTerCard() const
        {
            return hasTerCard_;
        }

        //////////////////////////////////////////////////////////
        //                       MUTATOR                        //
        //////////////////////////////////////////////////////////
        inline void AddTerCard()
        {
            hasTerCard_ = true;
        }

        inline void RemoveTerCard()
        {
            hasTerCard_ = false;
        }

        inline void setInsertionCode(const std::string s)
        {
            insertionCode_ = s;
        }

        inline void setChainId(const std::string s)
        {
            chainId_ = s;
        }

        cds::Atom* addPdbAtom(PdbData& data, size_t residueId, const std::string& line);
        cds::Atom* addPdbAtom(PdbData& data, size_t residueId, const std::string& name, const cds::Coordinate& c);
        void deletePdbAtom(PdbData& data, size_t residueId, cds::Atom* atom);

        //////////////////////////////////////////////////////////
        //                       FUNCTIONS                      //
        //////////////////////////////////////////////////////////
        void modifyNTerminal(PdbData& data, size_t residueId, const std::string& type);
        void modifyCTerminal(PdbData& data, size_t residueId, const std::string& type);

        //////////////////////////////////////////////////////////
        //                       DISPLAY FUNCTION               //
        //////////////////////////////////////////////////////////
        inline std::string printId() const
        {
            return this->getId().print();
        }

        void Print(std::ostream& out) const;

      private:
        //////////////////////////////////////////////////////////
        //                       ACCESSOR                       //
        //////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////
        //                       ATTRIBUTES                     //
        //////////////////////////////////////////////////////////
        std::string insertionCode_ = "";
        std::string chainId_       = "";
        bool hasTerCard_           = false;
    };
} // namespace pdb
#endif
