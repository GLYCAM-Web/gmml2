/* File: gmml2.i */
%module gmml2
%include <std_string.i>
%include <std_iostream.i>
%include<std_map.i>
%include<std_vector.i>

// #https://www.swig.org/Doc1.3/Customization.html#exception
%include exception.i
%exception
{
try	{
		$action
	} catch(const std::runtime_error& e) {
		SWIG_exception(SWIG_RuntimeError, e.what());
	} catch (const std::invalid_argument& e) {
   		SWIG_exception(SWIG_ValueError, e.what());
	} catch (const std::out_of_range& e) {
   		SWIG_exception(SWIG_IndexError, e.what());
	} catch(...) {
		SWIG_exception(SWIG_RuntimeError,"Unknown exception");
	}
}

%{
#define SWIG_FILE_WITH_INIT

#include "includes/CentralDataStructure/CondensedSequence/graphVizDotConfig.hpp"
#include "includes/CentralDataStructure/InternalPrograms/DrawGlycan/drawGlycan.hpp"
#include "includes/CentralDataStructure/InternalPrograms/Sequence/sequence.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbPreprocessorInputs.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbFile.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbResidueId.hpp"
#include "includes/CentralDataStructure/InternalPrograms/CarbohydrateBuilder/carbohydrateBuilder.hpp"
using namespace pdb;
%}

%inline %{
std::ostream & get_cout() { return std::cout; }
%}

///Naming conflicts///
//%rename(cds_PdbFile) pdb::PdbFile;
//%rename(cds_iPdbLineLength) pdb::iPdbLineLength;
//%rename(B_foo) B::foo;

//%include "includes/CodeUtils/constants.hpp"

// DrawGlycan
%include "includes/CentralDataStructure/CondensedSequence/graphVizDotConfig.hpp"
%include "includes/CentralDataStructure/InternalPrograms/DrawGlycan/drawGlycan.hpp"
%include "includes/CentralDataStructure/InternalPrograms/Sequence/sequence.hpp"
// MDPrep
%include "includes/CentralDataStructure/Readers/Pdb/pdbPreprocessorInputs.hpp"
%include "includes/CentralDataStructure/Readers/Pdb/pdbFile.hpp"
%include "includes/CentralDataStructure/Readers/Pdb/pdbResidueId.hpp"
// CarbohydrateBuilder
%include "includes/CentralDataStructure/InternalPrograms/CarbohydrateBuilder/carbohydrateBuilder.hpp"

///MD Prep///
%template (AtomInfoVector) std::vector<pdb::AtomInfo>;
%template (GapInAminoAcidChainVector) std::vector<pdb::GapInAminoAcidChain>;
%template (ResidueIdVector) std::vector<pdb::ResidueId>;
%template (ChainTerminalVector) std::vector<pdb::ChainTerminal>;
%template (DisulphideBondVector) std::vector<pdb::DisulphideBond>;
%template (hisSelectionPairVector) std::vector<std::pair<std::string,std::string>>;
