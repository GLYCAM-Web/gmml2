#ifndef INCLUDES_CENTRALDATASTRUCTURE_READERS_PDB_PDBDATA_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_READERS_PDB_PDBDATA_HPP

#include "includes/CentralDataStructure/cdsFunctions/graphInterface.hpp"

#include <string>
#include <vector>

namespace pdb
{
    struct AtomData
    {
        std::vector<std::string> recordNames;
        std::vector<double> occupancies;
        std::vector<double> temperatureFactors;
    };

    struct PdbData
    {
        AtomData atoms;
        std::vector<std::vector<size_t>> moleculeResidueOrder;
        cds::GraphIndexData indices;
    };
} // namespace pdb

#endif