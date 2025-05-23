#include "includes/MolecularMetadata/elements.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/strings.hpp"

#include <vector>
#include <map>
#include <utility>

namespace
{
    using MolecularMetadata::Element;
    const std::vector<std::string> elementNames = {
        "",   "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na", "Mg", "Al", "Si", "P",  "S",
        "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As",
        "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
        "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho",
        "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po",
        "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md",
        "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
    };

    struct FlaggedDouble
    {
        double value;
        bool valid = false;
    };

    std::vector<FlaggedDouble> withValues(const std::vector<std::pair<Element, double>>& values)
    {
        std::vector<FlaggedDouble> result(Element::ElementCount, {0.0, false});
        for (auto& pair : values)
        {
            result[pair.first] = {pair.second, true};
        }
        return result;
    }

    std::vector<FlaggedDouble> chimeraRadii = withValues({
  // taken from the worst-case of chimera's united or all atom radii, whichever is larger
        { Element::H, 1.000},
        { Element::C, 1.880},
        { Element::N, 1.640},
        { Element::O, 1.500},
        { Element::F, 1.560},
        { Element::P, 1.871},
        { Element::S, 1.782},
        {Element::Cl, 1.735},
        {Element::Br, 1.978},
        { Element::I, 2.094}
    });

    std::vector<FlaggedDouble> amberRadii = withValues({
  // taken from the worst-case of chimera's united or all atom radii, whichever is larger
        { Element::H, 1.2000},
        { Element::C, 1.9080},
        { Element::N, 1.8240},
        { Element::O, 1.6612},
        { Element::F, 1.7500},
        { Element::P, 2.1000},
        { Element::S, 2.0000},
        {Element::Cl, 1.9480},
        {Element::Br, 2.2200},
        { Element::I, 2.3500}
    });

    std::vector<bool> isHeavy =
        codeUtils::indicesToBools(Element::ElementCount, {Element::C, Element::O, Element::N, Element::S, Element::P});

    std::vector<FlaggedDouble> lennardJonesEpsilons = withValues({
        { Element::H, 4.47789},
        { Element::C, 6.36953},
        { Element::N, 9.75379},
        { Element::O, 5.12647},
        { Element::F, 1.60592},
        { Element::P, 5.03050},
        { Element::S, 4.36927},
        {Element::Cl, 4.48328},
        {Element::Br, 1.97063},
        { Element::I, 1.53938}
    });

    std::vector<FlaggedDouble> lennardJonesSigmas = withValues({
        { Element::H, 0.552357},
        { Element::C, 1.354170},
        { Element::N, 1.265080},
        { Element::O, 1.175990},
        { Element::F, 1.015620},
        { Element::P, 1.906520},
        { Element::S, 1.870890},
        {Element::Cl, 1.817430},
        {Element::Br, 2.138160},
        { Element::I, 2.476700}
    });

    auto values = [](const std::vector<FlaggedDouble>& vec)
    {
        std::vector<double> result(vec.size(), 0.0);
        for (size_t n = 0; n < vec.size(); n++)
        {
            if (vec[n].valid)
            {
                result[n] = vec[n].value;
            }
        }
        return result;
    };

    auto bools = [](const std::vector<FlaggedDouble>& vec)
    {
        std::vector<bool> result(vec.size(), false);
        for (size_t n = 0; n < vec.size(); n++)
        {
            result[n] = vec[n].valid;
        }
        return result;
    };

    MolecularMetadata::PotentialTable defaultPotentialTable {
        {values(lennardJonesEpsilons), bools(lennardJonesEpsilons)},
        {  values(lennardJonesSigmas),   bools(lennardJonesSigmas)}
    };

    codeUtils::SparseVector<double> amberAtomRadii   = {values(amberRadii), bools(amberRadii)};
    codeUtils::SparseVector<double> chimeraAtomRadii = {values(chimeraRadii), bools(chimeraRadii)};
} // namespace

MolecularMetadata::Element MolecularMetadata::toElement(const std::string& str)
{
    auto it = std::find(elementNames.begin(), elementNames.end(), str);
    return (it == elementNames.end()) ? Element::Unknown : Element(it - elementNames.begin());
}

const std::string& MolecularMetadata::elementName(Element element)
{
    if (element >= elementNames.size())
    {
        throw std::runtime_error("Element goes beyond range of handled elements: " + std::to_string(element));
    }
    return elementNames[element];
}

const codeUtils::SparseVector<double>& MolecularMetadata::amberVanDerWaalsRadii()
{
    return amberAtomRadii;
}

const codeUtils::SparseVector<double>& MolecularMetadata::chimeraVanDerWaalsRadii()
{
    return chimeraAtomRadii;
}

const codeUtils::SparseVector<double>& MolecularMetadata::defaultVanDerWaalsRadii()
{
    return chimeraVanDerWaalsRadii();
}

bool MolecularMetadata::isHeavyElement(Element element)
{
    return isHeavy[element];
}

const MolecularMetadata::PotentialTable& MolecularMetadata::potentialTable()
{
    return defaultPotentialTable;
}

void MolecularMetadata::validateElementsInPotentialTable(const PotentialTable& potential,
                                                         const std::vector<Element>& vec)
{
    for (Element a : vec)
    {
        const std::string& name = elementNames[a];
        if (!potential.epsilon.hasValue[a])
        {
            throw std::runtime_error("Missing Lennard-Jones epsilon for " + name);
        }
        if (!potential.sigma.hasValue[a])
        {
            throw std::runtime_error("Missing Lennard-Jones sigma for " + name);
        }
    }
}

double MolecularMetadata::potentialWeight(const PotentialTable& table, Element a, Element b)
{
    double epsilon = 0.5 * (table.epsilon.values[a] + table.epsilon.values[b]);
    double sigma   = 0.5 * (table.sigma.values[a] + table.sigma.values[b]);
    double sigma2  = sigma * sigma;
    double sigma4  = sigma2 * sigma2;
    return epsilon * sigma4 * sigma4 * sigma4;
}

MolecularMetadata::Element MolecularMetadata::findElementAtomicNumber(const std::string& queryElement)
{
    Element result = toElement(queryElement);
    if (result != Element::Unknown)
    {
        return result;
    }
    else
    {
        std::string message = "Did not find this Element in the list of atomic Elements: " + queryElement;
        gmml::log(__LINE__, __FILE__, gmml::ERR, message);
        throw std::runtime_error(message);
    }
}
