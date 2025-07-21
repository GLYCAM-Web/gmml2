#include "include/sequence/parsedResidue.hpp"

#include "include/util/casting.hpp"
#include "include/util/containers.hpp"
#include "include/util/logging.hpp"
#include "include/util/strings.hpp"

#include <sstream>

namespace gmml
{
    ParsedResidue::ParsedResidue(const sequence::ParsedResidueComponents& components_) : components(components_)
    {
        this->setName(components_.fullString);
        this->SetType(components_.type);
    }

    std::vector<ParsedResidue*> ParsedResidue::GetChildren() const
    {
        std::vector<ParsedResidue*> resRet;
        for (auto& currNodeRes : this->getChildren())
        {
            resRet.push_back(util::erratic_cast<ParsedResidue*>(currNodeRes));
        }
        return resRet;
    }

    std::vector<ParsedResidue*> ParsedResidue::GetParents() const
    {
        std::vector<ParsedResidue*> resRet;
        for (auto& currNodeRes : this->getParents())
        {
            resRet.push_back(util::erratic_cast<ParsedResidue*>(currNodeRes));
        }
        return resRet;
    }

    ParsedResidue* ParsedResidue::GetParent() const
    {
        std::vector<Residue*> parents = this->getParents();
        if (parents.empty())
        {
            std::string message = "Bad situation: Parsed residue named " + this->getName() + " has no parent\n";
            util::log(__LINE__, __FILE__, util::ERR, message);
            throw std::runtime_error(message);
        }
        return util::erratic_cast<ParsedResidue*>(parents.front());
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
} // namespace gmml
