#ifndef INCLUDES_CENTRALDATASTRUCTURE_WRITERS_OFFWRITER_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_WRITERS_OFFWRITER_HPP

#include "includes/CentralDataStructure/residue.hpp"

#include <vector>
#include <string>
#include <ostream>

namespace cds
{
    struct AtomOffData
    {
        AtomOffData(int number_, std::string name_, std::string type_, int atomicNumber_, double charge_,
                    Coordinate coordinate_, std::vector<int> children_);
        int number;
        std::string name;
        std::string type;
        int atomicNumber;
        double charge;
        Coordinate coordinate;
        std::vector<int> children;
    };

    struct ResidueOffData
    {
        ResidueOffData(int number_, std::string name_, ResidueType type_, std::vector<AtomOffData> atoms_,
                       std::vector<int> connections_);
        int number;
        std::string name;
        ResidueType type;
        std::vector<AtomOffData> atoms;
        std::vector<int> atomsConnectedToOtherResidues;
    };

    std::string getOffType(const cds::ResidueType queryType);
    void WriteOffFileUnit(const std::vector<ResidueOffData>& residues, std::ostream& stream,
                          const std::string unitName);
    void WriteResiduesToOffFile(std::vector<cds::Residue*> residues, std::ostream& stream);
    void WriteToOffFile(const std::vector<Residue*>& residues, std::ostream& stream, const std::string unitName);
} // namespace cds
#endif
