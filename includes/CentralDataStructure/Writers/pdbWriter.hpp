#ifndef INCLUDES_CENTRALDATASTRUCTURE_WRITERS_PDBWRITER_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_WRITERS_PDBWRITER_HPP

#include "includes/CentralDataStructure/molecule.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/atom.hpp"

#include <ostream>
#include <string>
#include <vector>

// Not all the options are available now, like e.g. writing molecule chain numbers, but I will add them as I need them
namespace cds
{
    struct AtomPdbData
    {
        AtomPdbData(const cds::Atom* atom, std::string recordName_, std::string residueName_, int residueNumber_,
                    std::string chainId_, std::string insertionCode_, double occupancy_, double temperatureFactor_);
        AtomPdbData(std::vector<cds::Atom*> atoms, std::vector<std::string> recordNames,
                    std::vector<std::string> residueNames, std::vector<int> residueNumbers);
        std::vector<Coordinate> coordinate;
        std::vector<int> number;
        std::vector<std::string> name;
        std::vector<std::string> element;
        std::vector<std::string> recordName;
        std::vector<std::string> residueAlternativeLocation;
        std::vector<int> residueNumber;
        std::vector<std::string> residueName;
        std::vector<std::string> chainId;
        std::vector<std::string> insertionCode;
        std::vector<double> occupancy;
        std::vector<double> temperatureFactor;
    };

    struct ResiduePdbData
    {
        ResiduePdbData(std::vector<std::vector<size_t>> atomIndices_, AtomPdbData atomData_)
            : atomIndices(atomIndices_), atomData(atomData_)
        {}

        std::vector<std::vector<size_t>> atomIndices;
        AtomPdbData atomData;
    };

    std::vector<bool> residueTER(const std::vector<ResidueType>& types);
    ResiduePdbData toResiduePdbData(std::vector<Residue*>& residues);

    void writeAssemblyToPdb(std::ostream& stream, const std::vector<cds::Molecule*> molecules);
    void writeMoleculeToPdb(std::ostream& stream, const std::vector<bool>& residueTER,
                            const std::vector<std::vector<size_t>>& indices, const AtomPdbData& atomData);
    void writeAtomToPdb(std::ostream& stream, const AtomPdbData& data, size_t index);
    void writeConectCards(std::ostream& stream, std::vector<std::pair<int, int>> atomsPairsConnectedToOtherResidues);
    void writeTrajectoryToPdb(std::ostream& stream, const std::vector<cds::Molecule*> molecules);

} // namespace cds
#endif /* INCLUDES_CENTRALDATASTRUCTURE_WRITERS_PDBWRITER_HPP_ */
