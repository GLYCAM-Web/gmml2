#include "includes/MolecularMetadata/elements.hpp"
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

    struct radius
    {
        double value;
        bool valid = false;
    };

    std::vector<radius> vanDerWaalsRadii()
    {
        std::vector<radius> result;
        result.resize(Element::ElementCount);

        auto setRadius = [&](Element element, double value)
        {
            result[element] = {value, true};
        };

        // taken from the worst-case of chimera's united or all atom radii, whichever is larger
        setRadius(Element::C, 1.88);
        setRadius(Element::N, 1.64);
        setRadius(Element::O, 1.5);
        setRadius(Element::S, 1.782);
        setRadius(Element::H, 1.0);
        setRadius(Element::P, 1.871);
        setRadius(Element::F, 1.560);
        setRadius(Element::Cl, 1.735);
        setRadius(Element::Br, 1.978);
        setRadius(Element::I, 2.094);

        return result;
    }

    std::vector<radius> atomRadii = vanDerWaalsRadii();
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

double MolecularMetadata::vanDerWaalsRadius(Element element)
{
    const radius& result = atomRadii[element];
    if (!result.valid)
    {
        std::string message = "No valid radius for element: " + std::to_string(element);
        gmml::log(__LINE__, __FILE__, gmml::ERR, message);
        throw std::runtime_error(message);
    }
    return result.value;
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
