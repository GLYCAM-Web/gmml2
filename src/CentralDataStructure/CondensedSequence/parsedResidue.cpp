#include "includes/CentralDataStructure/CondensedSequence/parsedResidue.hpp"
#include "includes/CodeUtils/casting.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/strings.hpp"
#include <sstream>

using cds::ResidueType;
using cdsCondensedSequence::ParsedResidue;
using cdsCondensedSequence::ParsedResidueComponents;

ParsedResidue::ParsedResidue(const ParsedResidueComponents& components_) : components(components_)
{
    this->setName(components_.fullString);
    this->SetType(components_.type);
}

std::string cdsCondensedSequence::getLink(cds::ResidueType type, const std::string& linkage)
{
    switch (type)
    {
        case (ResidueType::Sugar):
            return linkage.substr(linkage.size() - 1, 1);
        case (ResidueType::Derivative):
            return linkage.substr(linkage.size() - 1, 1);
        case (ResidueType::Deoxy):
            return linkage.substr(0, 1);
        default:
            return "0";
    }
}

std::vector<ParsedResidue*> ParsedResidue::GetChildren() const
{
    std::vector<ParsedResidue*> resRet;
    for (auto& currNodeRes : this->getChildren())
    {
        resRet.push_back(codeUtils::erratic_cast<ParsedResidue*>(currNodeRes));
    }
    return resRet;
}

std::vector<ParsedResidue*> ParsedResidue::GetParents() const
{
    std::vector<ParsedResidue*> resRet;
    for (auto& currNodeRes : this->getParents())
    {
        resRet.push_back(codeUtils::erratic_cast<ParsedResidue*>(currNodeRes));
    }
    return resRet;
}

ParsedResidue* ParsedResidue::GetParent() const
{
    std::vector<cds::Residue*> parents = this->getParents();
    if (parents.empty())
    {
        std::string message = "Bad situation: Parsed residue named " + this->getName() + " has no parent\n";
        gmml::log(__LINE__, __FILE__, gmml::ERR, message);
        throw std::runtime_error(message);
    }
    return codeUtils::erratic_cast<ParsedResidue*>(parents.front());
}

std::string ParsedResidue::GetChildLinkagesForGlycamResidueNaming() const
{
    const std::vector<ParsedResidue*> children = this->GetChildren();
    std::vector<std::string> nonDeoxyLinkNames;
    nonDeoxyLinkNames.reserve(children.size());
    for (auto child : children)
    {
        if (child->GetType() != ResidueType::Deoxy)
        { // For glycam residue name, e.g. 2YB, do not want deoxy linkages to impact the residue name.
            nonDeoxyLinkNames.push_back(child->GetLink());
        }
    }
    if (nonDeoxyLinkNames.empty())
    {
        return "Terminal";
    }
    else
    {
        return codeUtils::join(",", nonDeoxyLinkNames);
    }
}

std::string ParsedResidue::GetName() const
{
    return this->GetIsomer() + this->GetResidueName() + this->GetRingType() + this->GetResidueModifier() +
           this->GetRingShape();
}

std::string ParsedResidue::GetLinkageName() const
{ // Should only ever be zero or one inEdges in my current design.
    // e.g  Galb1-4[Glcb1-3]Manb1-OH. OH is the parent to Manb, the linkage is b1-.
    // Manb is parent (has out edges) to both Glc and Gal, but they only have one in edge each from Manb. The inedge
    // from Man to Gal has the linkage name 1-4.
    for (auto& linkage : this->getInEdges())
    {
        return linkage->getLabel();
    }
    if (this->getInEdges().empty() && this->GetType() == ResidueType::Sugar)
    { // ano-ano linkage as in Glca1-1Glcb, return the "b" of Glc
        return this->GetConfiguration();
    }
    return ""; // aglycones like ROH and OME will not have linkage info.
}
