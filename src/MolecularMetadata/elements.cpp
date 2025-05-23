#include "includes/MolecularMetadata/elements.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/parsing.hpp"
#include "includes/CodeUtils/strings.hpp"

#include <cmath>
#include <functional>
#include <vector>
#include <map>
#include <utility>
#include <optional>
#include <sstream>

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

    std::vector<FlaggedDouble> vdwRadii = withValues({
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
        { Element::H, 0.0157},
        { Element::C, 0.0860},
        { Element::N, 0.1700},
        { Element::O, 0.2100},
        { Element::F, 0.0610},
        { Element::P, 0.2000},
        { Element::S, 0.2500},
        {Element::Cl, 0.2650},
        {Element::Br, 0.4200},
        { Element::I, 0.5000}
    });

    std::vector<FlaggedDouble> mass = withValues({
        { Element::H,   1.0},
        { Element::C,  12.0},
        { Element::N,  14.0},
        { Element::O,  16.0},
        { Element::F,  19.0},
        { Element::P,  31.0},
        { Element::S,  32.1},
        {Element::Cl,  35.5},
        {Element::Se,  78.9},
        {Element::Br,  79.9},
        { Element::I, 126.9}
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

    codeUtils::SparseVector<double> atomicMass   = {values(mass), bools(mass)};
    codeUtils::SparseVector<double> atomVdwRadii = {values(vdwRadii), bools(vdwRadii)};
} // namespace

MolecularMetadata::ChemicalFormula MolecularMetadata::toFormula(const std::vector<Element>& elements)
{
    std::vector<int> count(ElementCount, 0);
    size_t nonzero = 0;
    for (Element element : elements)
    {
        nonzero += (count[element] == 0);
        count[element]++;
    }
    ChemicalFormula result;
    result.reserve(nonzero);
    for (size_t n = 0; n < count.size(); n++)
    {
        if (count[n] > 0)
        {
            result.push_back({Element(n), count[n]});
        }
    }
    return result;
}

// expects element names and numbers together, separated by space. Implicit 1 is not supported
// e.g "C3 H5 N1 O1 Se1"
MolecularMetadata::ChemicalFormula MolecularMetadata::parseFormula(const std::string& formula)
{
    auto isUpper = [](char c)
    {
        return std::isupper(c);
    };
    auto isDigit = [](char c)
    {
        return std::isdigit(c);
    };
    ChemicalFormula result;
    result.reserve(std::count_if(formula.begin(), formula.end(), isUpper));
    auto begin = formula.begin();
    auto from  = begin;
    while (from != formula.end())
    {
        auto until      = std::find(from, formula.end(), ' ');
        auto firstDigit = std::find_if(from, until, isDigit);
        MolecularMetadata::Element element =
            MolecularMetadata::toElement(formula.substr(from - begin, firstDigit - from));
        std::optional<int> count = codeUtils::parseInt(formula.substr(firstDigit - begin, until - firstDigit));
        result.push_back({element, count.value()});
        from = (until == formula.end()) ? formula.end() : until + 1;
    }
    return result;
};

MolecularMetadata::ChemicalFormula MolecularMetadata::formulaSum(const std::vector<ChemicalFormula>& formulas)
{
    std::vector<int> count(ElementCount, 0);
    for (auto& formula : formulas)
    {
        for (auto& a : formula)
        {
            count[a.first] += a.second;
        }
    }
    ChemicalFormula result;
    for (size_t n = 0; n < ElementCount; n++)
    {
        if (count[n] != 0)
        {
            result.push_back({Element(n), count[n]});
        }
    }
    return result;
}

std::string MolecularMetadata::toString(const ChemicalFormula& formula)
{
    std::ostringstream ss;
    for (auto& a : formula)
    {
        ss << elementName(a.first) << a.second;
    }
    return ss.str();
}

MolecularMetadata::Element MolecularMetadata::toElement(const std::string& str)
{
    auto it = std::find(elementNames.begin(), elementNames.end(), str);
    return (it == elementNames.end()) ? Element::Unknown : Element(it - elementNames.begin());
}

std::vector<bool> MolecularMetadata::foundElements(const std::vector<Element>& elements)
{
    std::function<size_t(const Element&)> toIndex = [](Element e)
    {
        return size_t(e);
    };
    return codeUtils::indicesToBools(ElementCount, codeUtils::vectorMap(toIndex, elements));
}

const std::string& MolecularMetadata::elementName(Element element)
{
    if (element >= elementNames.size())
    {
        throw std::runtime_error("Element goes beyond range of handled elements: " + std::to_string(element));
    }
    return elementNames[element];
}

codeUtils::SparseVector<double> MolecularMetadata::vanDerWaalsRadii()
{
    return atomVdwRadii;
}

bool MolecularMetadata::isHeavyElement(Element element)
{
    return isHeavy[element];
}

const codeUtils::SparseVector<double>& MolecularMetadata::elementMass()
{
    return atomicMass;
}

double MolecularMetadata::totalMass(const codeUtils::SparseVector<double>& mass, const ChemicalFormula& formula)
{
    double result = 0.0;
    for (auto& a : formula)
    {
        Element element = a.first;
        if (!mass.hasValue[element])
        {
            throw std::runtime_error("Error: missing atomic mass for element: " + elementName(element));
        }
        result += mass.values[element] * a.second;
    }
    return result;
}

MolecularMetadata::PotentialTable MolecularMetadata::potentialTable(const codeUtils::SparseVector<double>& radii,
                                                                    const std::vector<bool>& usedElements)
{
    std::vector<size_t> usedElementIndices = codeUtils::boolsToIndices(usedElements);
    size_t usedElementCount                = usedElementIndices.size();
    std::vector<size_t> indices(ElementCount, 0);
    size_t index = 0;
    for (size_t n : usedElementIndices)
    {
        indices[n] = index;
        index++;
    }
    const codeUtils::SparseVector<double>& epsilons = {values(lennardJonesEpsilons), bools(lennardJonesEpsilons)};
    std::vector<std::vector<PotentialFactor>> factors(usedElementCount,
                                                      std::vector<PotentialFactor>(usedElementCount, {0.0, 0.0}));
    for (size_t n = 0; n < usedElementCount; n++)
    {
        size_t indexN = usedElementIndices[n];
        for (size_t k = n; k < usedElementCount; k++)
        {
            size_t indexK = usedElementIndices[k];
            if (!radii.hasValue[indexK])
            {
                throw std::runtime_error("No radius found for element: " + elementName(Element(indexK)));
            }
            if (!epsilons.hasValue[indexK])
            {
                throw std::runtime_error("No lennard-jones epsilon found for element: " + elementName(Element(indexK)));
            }
            double epsilon = std::sqrt(epsilons.values[indexN] * epsilons.values[indexK]);
            double sigma   = radii.values[indexN] + radii.values[indexK];
            PotentialFactor factor {epsilon, sigma};
            factors[n][k] = factor;
            factors[k][n] = factor;
        }
    }
    return {indices, factors};
}

MolecularMetadata::PotentialFactor MolecularMetadata::potentialFactor(const PotentialTable& table, Element a, Element b)
{
    return table.factors[table.elementIndex[a]][table.elementIndex[b]];
}

double MolecularMetadata::lennardJonesPotential(const PotentialFactor& factor, double squaredDistance)
{
    double rmin  = factor.rmin;
    double pow2  = (rmin * rmin) / (std::numeric_limits<double>::epsilon() + squaredDistance);
    double pow4  = pow2 * pow2;
    double pow6  = pow4 * pow4;
    double pow12 = pow6 * pow6;
    return factor.epsilon * (pow12 - 2.0 * pow6);
}

double MolecularMetadata::distanceAtZerolennardJonesPotential(const PotentialFactor& factor)
{
    return std::pow(2.0, -1.0 / 6.0) * factor.rmin;
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
