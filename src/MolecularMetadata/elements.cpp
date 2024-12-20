#include "includes/MolecularMetadata/elements.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/strings.hpp"

#include <vector>
#include <map>
#include <utility>

namespace
{
    using MolecularMetadata::Element;
    std::map<std::string, Element> elementsByName = {
        { "H",  Element::H},
        {"He", Element::He},
        {"Li", Element::Li},
        {"Be", Element::Be},
        { "B",  Element::B},
        { "C",  Element::C},
        { "N",  Element::N},
        { "O",  Element::O},
        { "F",  Element::F},
        {"Ne", Element::Ne},
        {"Na", Element::Na},
        {"Mg", Element::Mg},
        {"Al", Element::Al},
        {"Si", Element::Si},
        { "P",  Element::P},
        { "S",  Element::S},
        {"Cl", Element::Cl},
        {"Ar", Element::Ar},
        { "K",  Element::K},
        {"Ca", Element::Ca},
        {"Sc", Element::Sc},
        {"Ti", Element::Ti},
        { "V",  Element::V},
        {"Cr", Element::Cr},
        {"Mn", Element::Mn},
        {"Fe", Element::Fe},
        {"Co", Element::Co},
        {"Ni", Element::Ni},
        {"Cu", Element::Cu},
        {"Zn", Element::Zn},
        {"Ga", Element::Ga},
        {"Ge", Element::Ge},
        {"As", Element::As},
        {"Se", Element::Se},
        {"Br", Element::Br},
        {"Kr", Element::Kr},
        {"Rb", Element::Rb},
        {"Sr", Element::Sr},
        { "Y",  Element::Y},
        {"Zr", Element::Zr},
        {"Nb", Element::Nb},
        {"Mo", Element::Mo},
        {"Tc", Element::Tc},
        {"Ru", Element::Ru},
        {"Rh", Element::Rh},
        {"Pd", Element::Pd},
        {"Ag", Element::Ag},
        {"Cd", Element::Cd},
        {"In", Element::In},
        {"Sn", Element::Sn},
        {"Sb", Element::Sb},
        {"Te", Element::Te},
        { "I",  Element::I},
        {"Xe", Element::Xe},
        {"Cs", Element::Cs},
        {"Ba", Element::Ba},
        {"La", Element::La},
        {"Ce", Element::Ce},
        {"Pr", Element::Pr},
        {"Nd", Element::Nd},
        {"Pm", Element::Pm},
        {"Sm", Element::Sm},
        {"Eu", Element::Eu},
        {"Gd", Element::Gd},
        {"Tb", Element::Tb},
        {"Dy", Element::Dy},
        {"Ho", Element::Ho},
        {"Er", Element::Er},
        {"Tm", Element::Tm},
        {"Yb", Element::Yb},
        {"Lu", Element::Lu},
        {"Hf", Element::Hf},
        {"Ta", Element::Ta},
        { "W",  Element::W},
        {"Re", Element::Re},
        {"Os", Element::Os},
        {"Ir", Element::Ir},
        {"Pt", Element::Pt},
        {"Au", Element::Au},
        {"Hg", Element::Hg},
        {"Tl", Element::Tl},
        {"Pb", Element::Pb},
        {"Bi", Element::Bi},
        {"Po", Element::Po},
        {"At", Element::At},
        {"Rn", Element::Rn},
        {"Fr", Element::Fr},
        {"Ra", Element::Ra},
        {"Ac", Element::Ac},
        {"Th", Element::Th},
        {"Pa", Element::Pa},
        { "U",  Element::U},
        {"Np", Element::Np},
        {"Pu", Element::Pu},
        {"Am", Element::Am},
        {"Cm", Element::Cm},
        {"Bk", Element::Bk},
        {"Cf", Element::Cf},
        {"Es", Element::Es},
        {"Fm", Element::Fm},
        {"Md", Element::Md},
        {"No", Element::No},
        {"Lr", Element::Lr},
        {"Rf", Element::Rf},
        {"Db", Element::Db},
        {"Sg", Element::Sg},
        {"Bh", Element::Bh},
        {"Hs", Element::Hs},
        {"Mt", Element::Mt},
        {"Ds", Element::Ds},
        {"Rg", Element::Rg},
        {"Cn", Element::Cn},
        {"Nh", Element::Nh},
        {"Fl", Element::Fl},
        {"Mc", Element::Mc},
        {"Lv", Element::Lv},
        {"Ts", Element::Ts},
        {"Og", Element::Og}
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
    const auto& it = elementsByName.find(str);
    return (it != elementsByName.end()) ? it->second : Element::Unknown;
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
