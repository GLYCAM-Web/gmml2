#ifndef INCLUDES_CENTRALDATASTRUCTURE_WRITERS_PDBWRITER_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_WRITERS_PDBWRITER_HPP

#include "includes/CentralDataStructure/molecule.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/atom.hpp"

#include <string>
#include <ostream>

// Not all the options are available now, like e.g. writing molecule chain numbers, but I will add them as I need them
namespace cds
{
    struct AtomPdbData
    {
        Coordinate coordinate;
        int number;
        std::string name;
        std::string element;
        std::string recordName;
        std::string residueAlternativeLocation;
        int residueNumber;
        std::string residueName;
        std::string chainId;
        std::string insertionCode;
        double occupancy;
        double temperatureFactor;
    };

    struct ResiduePdbData
    {
        ResidueType type;
        std::vector<AtomPdbData> atoms;
    };

    AtomPdbData toAtomPdbData(const cds::Atom* atom, std::string recordName, std::string residueName, int residueNumber,
                              std::string chainId, std::string insertionCode, double occupancy,
                              double temperatureFactor);

    AtomPdbData toAtomPdbData(const cds::Atom* atom, std::string recordName, std::string residueName,
                              int residueNumber);

    std::vector<bool> residueTER(const std::vector<ResidueType>& types);
    std::vector<AtomPdbData> residuePdbAtoms(Residue* residue);
    std::vector<std::vector<AtomPdbData>> residuePdbAtoms(std::vector<Residue*>& residues);

    void writeAssemblyToPdb(std::ostream& stream, const std::vector<cds::Molecule*> molecules);
    void writeMoleculeToPdb(std::ostream& stream, const std::vector<bool>& residueTER,
                            const std::vector<std::vector<AtomPdbData>>& residueAtoms);
    void writeResidueToPdb(std::ostream& stream, const std::vector<AtomPdbData>& atoms);
    void writeAtomToPdb(std::ostream& stream, const AtomPdbData& data);
    void writeConectCards(std::ostream& stream, std::vector<std::pair<int, int>> atomsPairsConnectedToOtherResidues);
    void writeTrajectoryToPdb(std::ostream& stream, const std::vector<cds::Molecule*> molecules);

} // namespace cds
#endif /* INCLUDES_CENTRALDATASTRUCTURE_WRITERS_PDBWRITER_HPP_ */
