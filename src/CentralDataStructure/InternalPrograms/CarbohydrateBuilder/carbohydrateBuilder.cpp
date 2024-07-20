#include "includes/CentralDataStructure/InternalPrograms/CarbohydrateBuilder/carbohydrateBuilder.hpp"
#include "includes/CentralDataStructure/Selections/shaperSelections.hpp" // cdsSelections
#include "includes/CentralDataStructure/Shapers/residueLinkage.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralShape.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralAngles.hpp"
#include "includes/CodeUtils/logging.hpp"

#include <vector>
#include <string>
#include <sstream>
#include <stdexcept>

using cdsCondensedSequence::carbohydrateBuilder;
//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
carbohydrateBuilder::carbohydrateBuilder(std::string condensedSequence)

try : carbohydrate_(condensedSequence)
{}

// If a ctor throws, even if you catch it the standard guarantees another throw. So this is just to make a message.
catch (const std::string& exceptionMessage)
{
    gmml::log(__LINE__, __FILE__, gmml::ERR,
              "carbohydrateBuilder class caught this exception message: " + exceptionMessage);
    throw std::runtime_error(exceptionMessage);
}
catch (const std::runtime_error& error)
{
    gmml::log(__LINE__, __FILE__, gmml::ERR, "carbohydrateBuilder class caught a runtime error:");
    gmml::log(__LINE__, __FILE__, gmml::ERR, error.what());
    throw error;
}
catch (...)
{
    std::string message = "carbohydrateBuilder class caught a throw that was not anticipated. Please report how you "
                          "got to this to glycam@gmail.com.";
    gmml::log(__LINE__, __FILE__, gmml::ERR, message);
    throw std::runtime_error(message);
}

//////////////////////////////////////////////////////////
//                       FUNCTIONS                      //
//////////////////////////////////////////////////////////
void carbohydrateBuilder::GenerateSingle3DStructureDefaultFiles(std::string fileOutputDirectory,
                                                                std::string outputFileNaming)
{
    this->carbohydrate_.Generate3DStructureFiles(fileOutputDirectory, outputFileNaming);
    return;
}

void carbohydrateBuilder::GenerateSpecific3DStructure(cdsCondensedSequence::SingleRotamerInfoVector conformerInfo,
                                                      std::string fileOutputDirectory)
{
    //     std::string linkageIndex; // What Dan is calling linkageLabel. Internal index determined at C++ level and
    //     given to frontend to track. std::string linkageName; // Can be whatever the user wants it to be, default to
    //     same as index. std::string dihedralName; // omg / phi / psi / chi1 / chi2 std::string selectedRotamer; // gg
    //     / tg / g- etc std::string numericValue; // user entered 64 degrees. Could be a v2 feature.
    // With a conformer (aka rotamerSet), setting will be different as each rotatable_dihedral will be set to e.g. "A",
    // whereas for linkages with combinatorial rotamers (e,g, phi -g/t, omg gt/gg/tg), we need to set each dihedral as
    // specified, but maybe it will be ok to go through and find the value for "A" in each rotatable dihedral.. yeah
    // actually it should be fine. Leaving comment for time being.
    for (auto& rotamerInfo : conformerInfo)
    {
        std::stringstream ss;
        ss << "linkage: " << rotamerInfo.linkageIndex << " "
           << this->convertIncomingRotamerNamesToStandard(rotamerInfo.dihedralName) << " being set to "
           << rotamerInfo.selectedRotamer << std::endl;
        gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
        int currentLinkageIndex = std::stoi(rotamerInfo.linkageIndex);
        cds::ResidueLinkage* currentLinkage =
            cdsSelections::selectLinkageWithIndex(this->carbohydrate_.GetGlycosidicLinkages(), currentLinkageIndex);
        std::string standardDihedralName = this->convertIncomingRotamerNamesToStandard(rotamerInfo.dihedralName);
        cds::setSpecificShape(currentLinkage->rotatableDihedrals, standardDihedralName, rotamerInfo.selectedRotamer);
    }
    std::string fileName = "structure";
    this->carbohydrate_.ResolveOverlaps();
    this->carbohydrate_.Generate3DStructureFiles(fileOutputDirectory, fileName);
    return;
}

std::string carbohydrateBuilder::GetNumberOfShapes(bool likelyShapesOnly) const
{
    return this->carbohydrate_.GetNumberOfShapes(likelyShapesOnly);
}

// Commenting out for as not being used, and will be confusing later. The front-end calls a differnt function that will
// build a single, specific rotamer.
void carbohydrateBuilder::GenerateUpToNRotamers(int maxRotamers)
{
    std::vector<cds::ResidueLinkage> linkagesOrderedForPermutation =
        cdsSelections::SplitLinkagesIntoPermutants(this->carbohydrate_.GetGlycosidicLinkages());
    this->generateLinkagePermutationsRecursively(linkagesOrderedForPermutation.begin(),
                                                 linkagesOrderedForPermutation.end(), maxRotamers);
}

cdsCondensedSequence::LinkageOptionsVector carbohydrateBuilder::GenerateUserOptionsDataStruct()
{
    cdsCondensedSequence::LinkageOptionsVector userOptionsForSequence;
    // carbohydrate_.SetIndexByConnectivity();
    for (auto& linkage : this->carbohydrate_.GetGlycosidicLinkages())
    {
        // std::cout << "linko nameo: " << linkage.GetName() << std::endl;
        cdsCondensedSequence::DihedralOptionsVector possibleRotamers, likelyRotamers;
        std::vector<std::string> buffer;
        auto multipleRotamerDihedrals = cds::rotatableDihedralsWithMultipleRotamers(linkage.rotatableDihedrals);
        for (auto& rotatableDihedral : multipleRotamerDihedrals)
        {
            for (auto& metadata : rotatableDihedral.metadataVector)
            {
                buffer.push_back(metadata.rotamer_name_);
            }
            std::string dihedralName = cds::likelyName(rotatableDihedral.metadataVector);
            possibleRotamers.emplace_back(dihedralName, buffer);
            buffer.clear();
            for (auto& metadata : cds::likelyMetadata(rotatableDihedral.metadataVector))
            {
                buffer.push_back(metadata.rotamer_name_);
            }
            likelyRotamers.emplace_back(dihedralName, buffer);
            buffer.clear();
        }
        // If there are multiple rotamers for this linkage
        if (!multipleRotamerDihedrals.empty())
        { // Build struct in vector with emplace_back via constructor in struct
            auto& linkageResidues = linkage.link.residues;
            userOptionsForSequence.emplace_back(
                linkage.name, std::to_string(linkage.index), std::to_string(linkageResidues.first->getNumber()),
                std::to_string(linkageResidues.second->getNumber()), likelyRotamers, possibleRotamers);
        }
    }
    return userOptionsForSequence;
}

//////////////////////////////////////////////////////////
//                   PRIVATE FUNCTIONS                  //
//////////////////////////////////////////////////////////
// Adapted from resolve_overlaps.
// Goes deep and then as it falls out of the iteration it's setting and writing the shapes.
void carbohydrateBuilder::generateLinkagePermutationsRecursively(std::vector<cds::ResidueLinkage>::iterator linkage,
                                                                 std::vector<cds::ResidueLinkage>::iterator end,
                                                                 int maxRotamers, int rotamerCount)
{
    for (int shapeNumber = 0; shapeNumber < cds::numberOfShapes(linkage->rotatableDihedrals); ++shapeNumber)
    {
        ++rotamerCount;
        if (rotamerCount <= maxRotamers)
        {
            cds::setSpecificShapeUsingMetadata(linkage->rotatableDihedrals, shapeNumber);
            if (std::next(linkage) != end)
            {
                this->generateLinkagePermutationsRecursively(std::next(linkage), end, maxRotamers, rotamerCount);
            }
            // else // At the end
            // {
            // Check for issues? Resolve
            // Figure out name of file:
            // http://128.192.9.183/eln/gwscratch/2020/01/10/succinct-rotamer-set-labeling-for-sequences/ Write PDB file
            this->carbohydrate_.Generate3DStructureFiles(".", std::to_string(rotamerCount));
            // }
        }
    }
    return;
}

// In too much of a rush to do this properly, so I'll make it private and
// so dumb that you'll have to write a proper one. Yes you!
// Metadata belongs in metadata, not in code.
// Should have put this at the point that the name is compared in
// bool Rotatable_dihedral::SetSpecificShape(std::string dihedralName, std::string selectedRotamer)
// Or somehow that comparison should be moved into metadata. I wonder could you overload the
// == comparison operator for the metadata struct to check through this list when checking the name?
std::string carbohydrateBuilder::convertIncomingRotamerNamesToStandard(std::string incomingName)
{ // Lamda function to see if string is in the passed in vector
    auto isAinList = [](std::vector<std::string> names, std::string query)
    {
        if (std::find(names.begin(), names.end(), query) != names.end())
        {
            return true;
        }
        return false;
    }; // Thought about converting incomingName to lowercase first.
    if (isAinList({"Omega", "omega", "Omg", "omg", "OMG", "omga", "omg1", "Omg1", "o"}, incomingName))
    {
        return "Omg";
    }
    if (isAinList({"Phi", "phi", "PHI", "h"}, incomingName))
    {
        return "Phi";
    }
    if (isAinList({"Psi", "psi", "PSI", "s"}, incomingName))
    {
        return "Psi";
    }
    if (isAinList({"Chi1", "chi1", "CHI1", "c1"}, incomingName))
    {
        return "Chi1";
    }
    if (isAinList({"Chi2", "chi2", "CHI2", "c2"}, incomingName))
    {
        return "Chi2";
    }
    if (isAinList({"Omega7", "omega7", "OMG7", "omga7", "Omg7", "omg7", "o7"}, incomingName))
    {
        return "Omg7";
    }
    if (isAinList({"Omega8", "omega8", "OMG8", "omga8", "Omg8", "omg8", "o8"}, incomingName))
    {
        return "Omg8";
    }
    if (isAinList({"Omega9", "omega9", "OMG9", "omga9", "Omg9", "omg9", "o9"}, incomingName))
    {
        return "Omg9";
    }
    std::stringstream ss;
    ss << "Specified rotamer name: \"" << incomingName
       << "\", is not recognized in convertIncomingRotamerNamesToStandard.";
    throw ss.str();
}
