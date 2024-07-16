#include "includes/MolecularMetadata/atomicBonds.hpp"

#include <array>
#include <vector>
#include <map>

namespace
{
    using MolecularMetadata::Element;
    typedef std::map<std::array<Element, 2>, std::pair<double, double>> BondLengthMap;

    BondLengthMap bondLengths()
    {
        Element C = Element::C;
        Element O = Element::O;
        Element N = Element::N;
        Element P = Element::P;
        Element S = Element::S;
        return {
            {{C, C},  {1.22, 1.67}},
            {{C, O},  {1.08, 1.68}},
            {{O, C},  {1.08, 1.68}},
            {{C, N},  {1.26, 1.55}},
            {{N, C},  {1.26, 1.55}},
            {{O, P}, {1.35, 1.776}},
            {{P, O}, {1.35, 1.776}},
            {{O, S},  {1.43, 1.78}},
            {{S, O},  {1.43, 1.78}},
            {{N, S},  {1.62, 1.77}},
            {{S, N},  {1.62, 1.77}},
            {{C, S},  {1.50, 1.91}},
            {{S, C},  {1.50, 1.91}},
            {{S, S},  {1.50, 2.10}}
        };
    }

    const BondLengthMap bondLengthMap = bondLengths();
} // namespace

std::pair<double, double> MolecularMetadata::getBondLengthByAtomType(Element atom1Element, Element atom2Element)
{ // Using PDB bond length statistics provided by Chenghua on 2/5/19

    std::array<Element, 2> bothAtoms = {atom1Element, atom2Element};

    if (bondLengthMap.find(bothAtoms) != bondLengthMap.end())
    {
        return bondLengthMap.at(bothAtoms);
    }
    else
    {
        // gmml::log(__LINE__, __FILE__,  gmml::INF, "Using default binding cutoff of 1.65");
        return {MolecularMetadata::minCutOff, MolecularMetadata::maxCutOff};
    }
}

double MolecularMetadata::getMaxBondLengthByAtomType(Element atom1Element, Element atom2Element)
{ // Using PDB bond length statistics provided by Chenghua on 2/5/19

    std::pair<double, double> result = getBondLengthByAtomType(atom1Element, atom2Element);
    return result.second;
}
