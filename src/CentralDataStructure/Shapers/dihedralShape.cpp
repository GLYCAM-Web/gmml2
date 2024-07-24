#include "includes/CentralDataStructure/Shapers/dihedralShape.hpp"
#include "includes/CentralDataStructure/Geometry/rotationMatrix.hpp"
#include "includes/CentralDataStructure/Geometry/coordinate.hpp"
#include "includes/CentralDataStructure/atom.hpp"

using cds::RotatableDihedral;

void cds::setDihedralAngle(RotatableDihedral& dihedral, cds::AngleWithMetadata target)
{
    auto matrix = rotationTo(dihedralCoordinates(dihedral), constants::degree2Radian(target.value));
    matrix.rotateCoordinates(dihedral.movingCoordinates);
    dihedral.currentMetadata = target.metadata;
}

bool cds::setSpecificShape(RotatableDihedral& dihedral, std::string dihedralName, std::string selectedRotamer)
{
    auto& metadataVector = dihedral.metadataVector;
    if (dihedralName == metadataVector[0].dihedral_angle_name_)
    {
        for (auto& metadata : metadataVector)
        {
            if (metadata.rotamer_name_ == selectedRotamer)
            {
                setDihedralAngle(dihedral, {metadata.default_angle_value_, metadata});
                return true;
            }
        }
    }
    return false;
}

std::vector<cds::AngleWithMetadata> cds::currentShape(const std::vector<RotatableDihedral>& dihedrals)
{
    std::vector<AngleWithMetadata> result;
    result.reserve(dihedrals.size());

    for (auto& dihedral : dihedrals)
    {
        result.push_back({cds::angle(dihedralCoordinates(dihedral)), dihedral.currentMetadata});
    }

    return result;
}

std::vector<std::vector<cds::AngleWithMetadata>> cds::currentShape(const std::vector<ResidueLinkage>& linkages)
{
    std::vector<std::vector<AngleWithMetadata>> result;
    result.reserve(linkages.size());

    for (auto& a : linkages)
    {
        result.push_back(currentShape(a.rotatableDihedrals));
    }

    return result;
}

void cds::setShape(std::vector<ResidueLinkage>& linkages, const std::vector<std::vector<AngleWithMetadata>>& angles)
{
    for (size_t n = 0; n < linkages.size(); n++)
    {
        setShape(linkages[n].rotatableDihedrals, angles[n]);
    }
}

void cds::setShape(std::vector<cds::RotatableDihedral>& dihedrals, const std::vector<AngleWithMetadata>& angles)
{
    for (size_t n = 0; n < angles.size(); n++)
    {
        setDihedralAngle(dihedrals[n], angles[n]);
    }
}

void cds::setDefaultShapeUsingMetadata(std::vector<cds::RotatableDihedral>& dihedrals)
{
    for (auto& dihedral : dihedrals)
    {
        setDihedralAngle(dihedral, defaultAngle(dihedral.metadataVector[0]));
    }
}

void cds::setRandomShapeUsingMetadata(std::vector<RotatableDihedral>& dihedrals)
{
    auto rotamerType = dihedrals.at(0).metadataVector.at(0).rotamer_type_;
    if (rotamerType == gmml::MolecularMetadata::GLYCAM::RotamerType::permutation)
    {
        for (auto& entry : dihedrals)
        {
            auto angle = randomAngleEntryUsingMetadata(entry.metadataVector);
            setDihedralAngle(entry, angle);
        }
    }
    else if (rotamerType == gmml::MolecularMetadata::GLYCAM::RotamerType::conformer)
    {
        int numberOfConformers = dihedrals.at(0).metadataVector.size();
        std::uniform_int_distribution<> distr(0, (numberOfConformers - 1)); // define the range
        int randomlySelectedConformerNumber = distr(rng);
        for (auto& entry : dihedrals)
        {
            auto angle = randomDihedralAngleWithinMetadataRange(entry.metadataVector[randomlySelectedConformerNumber]);
            setDihedralAngle(entry, angle);
        }
    }
}

void cds::setRandomShapeUsingMetadata(std::vector<ResidueLinkage>& linkages)
{
    for (auto& linkage : linkages)
    {
        setRandomShapeUsingMetadata(linkage.rotatableDihedrals);
    }
}

void cds::setSpecificShapeUsingMetadata(std::vector<RotatableDihedral>& dihedrals, int shapeNumber)
{
    for (auto& entry : dihedrals)
    {
        setDihedralAngle(entry, defaultAngle(entry.metadataVector[shapeNumber]));
    }
}

void cds::setSpecificShape(std::vector<RotatableDihedral>& dihedrals, std::string dihedralName,
                           std::string selectedRotamer)
{
    for (auto& rotatableDihedral : dihedrals)
    {
        // This will call RotatableDihedrals that don't have dihedralName (phi,psi), and nothing will happen. Hmmm.
        if (setSpecificShape(rotatableDihedral, dihedralName, selectedRotamer))
        {
            return; // Return once you manage to set a shape.
        }
    }
    std::string errorMessage = "Did not set " + dihedralName + " to " + selectedRotamer +
                               " as requested in ResidueLinkage::SetSpecificShape()";
    gmml::log(__LINE__, __FILE__, gmml::ERR, errorMessage);
    throw std::runtime_error(errorMessage);
}
