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
        AtomPdbData(std::vector<cds::Atom*>& atoms, std::vector<std::string> recordNames_,
                    std::vector<double> occupancies_, std::vector<double> temperatureFactors_);
        AtomPdbData(std::vector<cds::Atom*>& atoms_, std::vector<std::string> recordNames_);
        std::vector<cds::Atom*> atoms;
        std::vector<Coordinate> coordinates;
        std::vector<int> numbers;
        std::vector<std::string> names;
        std::vector<std::string> elements;
        std::vector<std::string> recordNames;
        std::vector<double> occupancies;
        std::vector<double> temperatureFactors;
    };

    struct ResiduePdbData
    {
        ResiduePdbData(std::vector<std::vector<size_t>> atomIndices_, std::vector<int> numbers_,
                       std::vector<std::string> names_, std::vector<std::string> chainIds_,
                       std::vector<std::string> insertionCodes_);

        std::vector<std::vector<size_t>> atomIndices;
        std::vector<int> numbers;
        std::vector<std::string> names;
        std::vector<std::string> chainIds;
        std::vector<std::string> insertionCodes;
    };

    struct PdbWriterData
    {
        ResiduePdbData residues;
        AtomPdbData atoms;
    };

    std::vector<bool> residueTER(const std::vector<ResidueType>& types);
    PdbWriterData toPdbWriterData(std::vector<Residue*>& residues);

    void writeMoleculeToPdb(std::ostream& stream, const std::vector<size_t>& residueIndices,
                            const std::vector<bool>& residueTER, const PdbWriterData& data);
    void writeAtomToPdb(std::ostream& stream, const ResiduePdbData& residues, size_t residueIndex,
                        const AtomPdbData& atoms, size_t atomIndex);
    void writeConectCards(std::ostream& stream, const std::vector<int>& atomNumbers,
                          std::vector<std::pair<size_t, size_t>> connectionIndices);
    void writeTrajectoryToPdb(std::ostream& stream, const std::vector<cds::Molecule*> molecules);

} // namespace cds
#endif
