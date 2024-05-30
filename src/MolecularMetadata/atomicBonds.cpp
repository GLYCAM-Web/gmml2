#include "includes/MolecularMetadata/atomicBonds.hpp"
#include <map>

using namespace atomicBonds;

namespace
{
    static const std::map<std::string, std::pair<double, double>> bondLengthMap = {
        {"CC", std::make_pair(1.22,  1.67)},
        {"CO", std::make_pair(1.08,  1.68)},
        {"OC", std::make_pair(1.08,  1.68)},
        {"CN", std::make_pair(1.26,  1.55)},
        {"NC", std::make_pair(1.26,  1.55)},
        {"OP", std::make_pair(1.35, 1.776)},
        {"PO", std::make_pair(1.35, 1.776)},
        {"OS", std::make_pair(1.43,  1.78)},
        {"SO", std::make_pair(1.43,  1.78)},
        {"NS", std::make_pair(1.62,  1.77)},
        {"SN", std::make_pair(1.62,  1.77)},
        {"CS", std::make_pair(1.50,  1.91)},
        {"SC", std::make_pair(1.50,  1.91)},
        {"SS", std::make_pair(1.50,  2.10)}
    };
}

std::pair<double, double> atomicBonds::getBondLengthByAtomType(const std::string& atom1Element,
                                                               const std::string& atom2Element)
{ // Using PDB bond length statistics provided by Chenghua on 2/5/19

    std::string bothAtoms = atom1Element + atom2Element;

    if (bondLengthMap.find(bothAtoms) != bondLengthMap.end())
    {
        std::pair<double, double> cutoffDistances = bondLengthMap.at(bothAtoms);
        return cutoffDistances;
    }
    else
    {
        // gmml::log(__LINE__, __FILE__,  gmml::INF, "Using default binding cutoff of 1.65");
        return std::make_pair(atomicBonds::minCutOff, atomicBonds::maxCutOff);
    }
}

double atomicBonds::getMaxBondLengthByAtomType(const std::string& atom1Element, const std::string& atom2Element)
{ // Using PDB bond length statistics provided by Chenghua on 2/5/19

    std::pair<double, double> result = atomicBonds::getBondLengthByAtomType(atom1Element, atom2Element);
    return result.second;
}
