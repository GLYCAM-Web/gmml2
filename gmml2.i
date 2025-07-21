/* File: gmml2.i */
%module gmml2
%include <std_string.i>
%include <std_iostream.i>
%include <std_map.i>
%include <std_vector.i>

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

#include "include/sequence/graphViz.hpp"
#include "include/internalPrograms/DrawGlycan/drawGlycan.hpp"
#include "include/sequence/sequence.hpp"
#include "include/readers/Pdb/pdbPreprocessorInputs.hpp"
#include "include/readers/Pdb/pdbFile.hpp"
#include "include/readers/Pdb/pdbResidueId.hpp"
#include "include/internalPrograms/CarbohydrateBuilder/carbohydrateBuilder.hpp"
using namespace pdb;
using namespace gmml;
%}

%inline %{
std::ostream & get_cout() { return std::cout; }
%}

// DrawGlycan
%include "include/sequence/graphViz.hpp"
%include "include/internalPrograms/DrawGlycan/drawGlycan.hpp"
%include "include/sequence/sequence.hpp"
// MDPrep
%include "include/readers/Pdb/pdbPreprocessorInputs.hpp"
%include "include/readers/Pdb/pdbFile.hpp"
%include "include/readers/Pdb/pdbResidueId.hpp"
// CarbohydrateBuilder
%include "include/internalPrograms/CarbohydrateBuilder/carbohydrateBuilder.hpp"

///MD Prep///
%template (AtomInfoVector) std::vector<pdb::AtomInfo>;
%template (GapInAminoAcidChainVector) std::vector<pdb::GapInAminoAcidChain>;
%template (ResidueIdVector) std::vector<pdb::ResidueId>;
%template (ChainTerminalVector) std::vector<pdb::ChainTerminal>;
%template (DisulphideBondVector) std::vector<pdb::DisulphideBond>;
%template (hisSelectionPairVector) std::vector<std::pair<std::string,std::string>>;


///Carbohydrate Builder///
%template(string_vector) std::vector<std::string>;
%template(dihedral_options_vector) std::vector<gmml::DihedralOptions>;
%template(linkage_options_vector) std::vector<gmml::LinkageOptions>;
%template(single_rotamer_info_vector) std::vector<gmml::SingleRotamerInfo>;
