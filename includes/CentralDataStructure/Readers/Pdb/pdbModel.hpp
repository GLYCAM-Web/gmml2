#ifndef INCLUDES_CENTRALDATASTRUCTURE_READERS_PDB_PDBMODEL_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_READERS_PDB_PDBMODEL_HPP
#include "includes/CentralDataStructure/assembly.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbPreprocessorInputs.hpp"
#include "includes/CentralDataStructure/Parameters/parameterManager.hpp"
#include <vector>
#include <sstream>
#include <ostream>

namespace pdb
{
    class PdbAtom;
    class PdbResidue;
    class PdbChain;

    class PdbModel : public cds::Assembly
    {
      public:
        //////////////////////////////////////////////////////////
        //                       CONSTRUCTOR                    //
        //////////////////////////////////////////////////////////
        PdbModel();
        PdbModel(std::stringstream& stream_block);
        //////////////////////////////////////////////////////////
        //                       FUNCTIONS                      //
        //////////////////////////////////////////////////////////
        void ChangeResidueName(const std::string& selector, const std::string& newName);

        // Preprocessing functions
        void preProcessCysResidues(pdb::PreprocessorInformation& ppInfo);
        void preProcessHisResidues(pdb::PreprocessorInformation& ppInfo, const pdb::PreprocessorOptions& inputOptions);
        void preProcessChainTerminals(pdb::PreprocessorInformation& ppInfo,
                                      const pdb::PreprocessorOptions& inputOptions);
        void preProcessGapsUsingDistance(pdb::PreprocessorInformation& ppInfo,
                                         const pdb::PreprocessorOptions& inputOptions);
        void preProcessMissingUnrecognized(pdb::PreprocessorInformation& ppInfo,
                                           const cdsParameters::ParameterManager& parmManager);
        //        void bondAtomsByDistance();
        //////////////////////////////////////////////////////////
        //                       DISPLAY FUNCTION               //
        //////////////////////////////////////////////////////////
        void Print(std::ostream& out) const;
        void Write(std::ostream& stream) const;

      private:
        //////////////////////////////////////////////////////////
        //                       PRIVATE FUNCTIONS              //
        //////////////////////////////////////////////////////////
        std::string extractChainId(const std::string& line);
        std::stringstream extractSingleChainFromRecordSection(std::stringstream& stream_block, std::string line,
                                                              const std::string& initialChainID);
        void extractCoordinatesFromModel(std::stringstream& stream_block, std::string line);
    };
} // namespace pdb
#endif
