#ifndef INCLUDES_CENTRALDATASTRUCTURE_READERS_PDB_PDBCHAIN_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_READERS_PDB_PDBCHAIN_HPP
#include "includes/CentralDataStructure/Readers/Pdb/pdbResidue.hpp"
#include "includes/CentralDataStructure/molecule.hpp"
#include <string>
#include <iostream>

// #include <functional>
namespace pdb
{
    class PdbAtom;
    class PdbResidue;

    class PdbChain : public cds::Molecule
    {
      public:
        //////////////////////////////////////////////////////////
        //                       CONSTRUCTOR                    //
        //////////////////////////////////////////////////////////
        PdbChain(std::stringstream& stream_block, const std::string& chainId);
        //    PdbChain(PdbResidue pdbResidue);
        //    PdbChain(std::vector<PdbResidue> pdbResidues);
        //////////////////////////////////////////////////////////
        //                       ACCESSOR                       //
        //////////////////////////////////////////////////////////
        //    const std::string& GetChainId() const;
        //////////////////////////////////////////////////////////
        //                       MUTATOR                        //
        //////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////
        //                       FUNCTIONS                      //
        //////////////////////////////////////////////////////////
        void tagTerminalResidues();
        void InsertCap(const PdbResidue& refResidue, const std::string& type);
        void ModifyTerminal(const std::string& type, PdbResidue* terminalResidue);
        PdbResidue* getNTerminal();
        PdbResidue* getCTerminal();

      private:
        //////////////////////////////////////////////////////////
        //                       FUNCTIONS                      //
        //////////////////////////////////////////////////////////
        // std::string extractResidueId(const std::string &line);
        std::stringstream extractSingleResidueFromRecordSection(std::stringstream& pdbFileStream, std::string line);

        // PdbResidue& GetFirstResidue() const;
        //////////////////////////////////////////////////////////
        //                       ATTRIBUTES                     //
        //////////////////////////////////////////////////////////
        // std::vector<PdbResidue> pdbResidues_;
        // std::string chainId_;
    };
} // namespace pdb
#endif
