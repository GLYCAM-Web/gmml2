#ifndef INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_SEQUENCEMANIPULATOR_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_SEQUENCEMANIPULATOR_HPP

#include "includes/CentralDataStructure/CondensedSequence/parsedResidue.hpp"
#include "includes/CentralDataStructure/CondensedSequence/graphVizDotConfig.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/molecule.hpp"

#include <vector>
#include <string>

namespace cdsCondensedSequence
{
    std::vector<ParsedResidue*> parsedResiduesOrderedByConnectivity(std::vector<cds::Residue*> residues);
    void setIndexByConnectivity(std::vector<cds::Residue*> residues);
    void labelSequence(std::vector<cds::Residue*> residues);
    std::string printGraphViz(GraphVizDotConfig& configs, std::vector<cds::Residue*> residues);
    std::string printSequence(std::vector<cds::Residue*> residues, bool withLabels);
    std::string reorderSequence(cds::Molecule* molecule);
    ParsedResidue* terminalResidue(std::vector<cds::Residue*> residues);
} // namespace cdsCondensedSequence
#endif
