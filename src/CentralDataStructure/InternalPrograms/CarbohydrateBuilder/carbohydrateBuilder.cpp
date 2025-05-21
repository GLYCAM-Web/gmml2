#include "includes/CentralDataStructure/InternalPrograms/CarbohydrateBuilder/carbohydrateBuilder.hpp"
#include "includes/CentralDataStructure/Parameters/parameterManager.hpp"
#include "includes/CentralDataStructure/Selections/shaperSelections.hpp" // cdsSelections
#include "includes/CentralDataStructure/Shapers/residueLinkage.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralShape.hpp"
#include "includes/CentralDataStructure/Shapers/dihedralAngleSearch.hpp"
#include "includes/CentralDataStructure/CondensedSequence/sequenceParser.hpp"
#include "includes/MolecularMetadata/GLYCAM/dihedralangledata.hpp"
#include "includes/MolecularMetadata/elements.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/containerTypes.hpp"
#include "includes/CodeUtils/filesystem.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/version.h"

#include <vector>
#include <string>
#include <sstream>
#include <stdexcept>

using cdsCondensedSequence::carbohydrateBuilder;
//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
carbohydrateBuilder::carbohydrateBuilder(std::string condensedSequence)

try : parameters_(cdsParameters::loadParameters(codeUtils::gmmlHomeDirPath)),
    carbohydrate_(parameters_, MolecularMetadata::defaultVanDerWaalsRadii(),
                  cdsCondensedSequence::parseAndReorder(condensedSequence))
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
    this->Generate3DStructureFiles(fileOutputDirectory, outputFileNaming);
}

void carbohydrateBuilder::Generate3DStructureFiles(const std::string& fileOutputDirectory,
                                                   const std::string& outputFileNaming)
{
    std::vector<std::string> headerLines {"Produced by GMML (https://github.com/GLYCAM-Web/gmml2)  version " +
                                          std::string(GMML_VERSION)};
    this->carbohydrate_.Generate3DStructureFiles(fileOutputDirectory, outputFileNaming, headerLines);
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
    const codeUtils::SparseVector<double>& elementRadii         = MolecularMetadata::defaultVanDerWaalsRadii();
    const GlycamMetadata::DihedralAngleDataTable& metadataTable = GlycamMetadata::dihedralAngleDataTable();
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
        cds::setSpecificShape(metadataTable, currentLinkage->rotatableDihedrals, currentLinkage->dihedralMetadata,
                              standardDihedralName, rotamerInfo.selectedRotamer);
    }
    std::string fileName = "structure";
    this->carbohydrate_.ResolveOverlaps(elementRadii, metadataTable, cdsCondensedSequence::defaultSearchSettings);
    this->Generate3DStructureFiles(fileOutputDirectory, fileName);
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
    const GlycamMetadata::DihedralAngleDataTable metadataTable = GlycamMetadata::dihedralAngleDataTable();
    std::vector<cds::ResidueLinkage> linkagesOrderedForPermutation =
        cdsSelections::SplitLinkagesIntoPermutants(this->carbohydrate_.GetGlycosidicLinkages());
    this->generateLinkagePermutationsRecursively(metadataTable, linkagesOrderedForPermutation.begin(),
                                                 linkagesOrderedForPermutation.end(), maxRotamers);
}

cdsCondensedSequence::LinkageOptionsVector carbohydrateBuilder::GenerateUserOptionsDataStruct()
{
    const GlycamMetadata::DihedralAngleDataTable metadataTable = GlycamMetadata::dihedralAngleDataTable();
    cdsCondensedSequence::LinkageOptionsVector userOptionsForSequence;
    // carbohydrate_.SetIndexByConnectivity();
    for (auto& linkage : cds::nonDerivativeResidueLinkages(this->carbohydrate_.GetGlycosidicLinkages()))
    {
        // std::cout << "linko nameo: " << linkage.GetName() << std::endl;
        cdsCondensedSequence::DihedralOptionsVector possibleRotamers, likelyRotamers;
        std::vector<std::string> buffer;
        std::vector<size_t> multipleRotamerDihedralIndices =
            cds::rotatableDihedralsWithMultipleRotamers(linkage.dihedralMetadata);
        for (size_t index : multipleRotamerDihedralIndices)
        {
            auto& metadataVector = linkage.dihedralMetadata[index];
            for (auto& metadata : metadataVector)
            {
                buffer.push_back(metadataTable.entries[metadata].rotamer_name_);
            }
            std::vector<size_t> likely = GlycamMetadata::likelyMetadata(metadataTable, metadataVector);
            std::string dihedralName   = likely.empty() ? "Boo" : metadataTable.entries[likely[0]].dihedral_angle_name_;
            possibleRotamers.emplace_back(dihedralName, buffer);
            buffer.clear();
            for (auto& metadata : likely)
            {
                buffer.push_back(metadataTable.entries[metadata].rotamer_name_);
            }
            likelyRotamers.emplace_back(dihedralName, buffer);
            buffer.clear();
        }
        // If there are multiple rotamers for this linkage
        if (!multipleRotamerDihedralIndices.empty())
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
void carbohydrateBuilder::generateLinkagePermutationsRecursively(
    const GlycamMetadata::DihedralAngleDataTable& metadataTable, std::vector<cds::ResidueLinkage>::iterator linkage,
    std::vector<cds::ResidueLinkage>::iterator end, int maxRotamers, int rotamerCount)
{
    auto defaultAngle = [](const GlycamMetadata::DihedralAngleData metadata)
    {
        return metadata.default_angle;
    };

    for (size_t shapeNumber = 0;
         shapeNumber < cds::numberOfShapes(metadataTable, linkage->rotamerType, linkage->dihedralMetadata);
         ++shapeNumber)
    {
        ++rotamerCount;
        if (rotamerCount <= maxRotamers)
        {
            std::function<std::vector<size_t>(const GlycamMetadata::DihedralAngleDataTable&, const std::vector<size_t>)>
                specificShape =
                    [&shapeNumber](const GlycamMetadata::DihedralAngleDataTable&, const std::vector<size_t>&)
            {
                return std::vector<size_t> {shapeNumber};
            };
            cds::setShapeToPreference(*linkage,
                                      cds::linkageShapePreference(specificShape, defaultAngle, metadataTable,
                                                                  linkage->rotamerType, linkage->dihedralMetadata));
            if (std::next(linkage) != end)
            {
                this->generateLinkagePermutationsRecursively(metadataTable, std::next(linkage), end, maxRotamers,
                                                             rotamerCount);
            }
            // else // At the end
            // {
            // Check for issues? Resolve
            // Figure out name of file:
            // http://128.192.9.183/eln/gwscratch/2020/01/10/succinct-rotamer-set-labeling-for-sequences/ Write PDB file
            this->Generate3DStructureFiles(".", std::to_string(rotamerCount));
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
        return codeUtils::contains(names, query);
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
