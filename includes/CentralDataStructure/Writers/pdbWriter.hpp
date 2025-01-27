#ifndef INCLUDES_CENTRALDATASTRUCTURE_WRITERS_PDBWRITER_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_WRITERS_PDBWRITER_HPP

#include "includes/CentralDataStructure/FileFormats/pdbFileData.hpp"

#include "includes/CentralDataStructure/molecule.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/atom.hpp"

#include <ostream>
#include <string>
#include <vector>

// Not all the options are available now, like e.g. writing molecule chain numbers, but I will add them as I need them
namespace cds
{
    std::vector<bool> residueTER(const std::vector<ResidueType>& types);
    PdbFileAtomData toPdbFileAtomData(const std::vector<cds::Atom*>& atoms, std::vector<std::string> recordNames);
    PdbFileAtomData toPdbFileAtomData(const std::vector<cds::Atom*>& atoms, std::vector<std::string> recordNames,
                                      std::vector<double> occupancies, std::vector<double> temperatureFactors);
    PdbFileData toPdbFileData(std::vector<Residue*>& residues);

    void writeTrajectoryToPdb(std::ostream& stream, const std::vector<cds::Molecule*>& molecules);
    void WritePdb(std::ostream& stream, cds::Molecule* molecule);
} // namespace cds
#endif
