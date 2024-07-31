#include "includes/CentralDataStructure/Shapers/dihedralShape.hpp"
#include "includes/CentralDataStructure/Geometry/rotationMatrix.hpp"
#include "includes/CentralDataStructure/Geometry/coordinate.hpp"
#include "includes/CentralDataStructure/atom.hpp"

using cds::RotatableDihedral;

void cds::setDihedralAngle(RotatableDihedral& dihedral, cds::AngleWithMetadata target)
{
    auto matrix = rotationTo(dihedralCoordinates(dihedral), constants::degree2Radian(target.value));
    matrix.rotateCoordinates(dihedral.movingCoordinates);
    dihedral.currentMetadataIndex = target.metadataIndex;
}

bool cds::setSpecificShape(RotatableDihedral& dihedral, std::string dihedralName, std::string selectedRotamer)
{
    auto& metadataVector = dihedral.metadataVector;
    if (dihedralName == metadataVector[0].dihedral_angle_name_)
    {
        for (size_t n = 0; n < metadataVector.size(); n++)
        {
            auto& metadata = metadataVector[n];
            if (metadata.rotamer_name_ == selectedRotamer)
            {
                setDihedralAngle(dihedral, {metadata.default_angle_value_, metadata.default_angle_value_, n});
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
        auto& metadataIndex   = dihedral.currentMetadataIndex;
        auto& currentMetadata = dihedral.metadataVector[metadataIndex];
        result.push_back(
            {cds::angle(dihedralCoordinates(dihedral)), currentMetadata.default_angle_value_, metadataIndex});
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
        double defaultAngle = dihedral.metadataVector[0].default_angle_value_;
        setDihedralAngle(dihedral, {defaultAngle, defaultAngle, 0});
    }
}

void cds::setRandomShapeUsingMetadata(MetadataDistribution randomMetadata, AngleDistribution randomAngle,
                                      gmml::MolecularMetadata::GLYCAM::RotamerType rotamerType,
                                      std::vector<RotatableDihedral>& dihedrals)
{
    auto randomAngleWithMetadata = [&](size_t index, const DihedralAngleData& metadata)
    {
        double angle = randomAngle(metadata);
        return AngleWithMetadata {angle, metadata.default_angle_value_, index};
    };

    if (rotamerType == gmml::MolecularMetadata::GLYCAM::RotamerType::permutation)
    {
        for (auto& entry : dihedrals)
        {
            size_t index  = randomMetadata(entry.metadataVector);
            auto metadata = entry.metadataVector[index];
            setDihedralAngle(entry, randomAngleWithMetadata(index, metadata));
        }
    }
    else if (rotamerType == gmml::MolecularMetadata::GLYCAM::RotamerType::conformer)
    {
        size_t conformerIndex = randomMetadata(dihedrals.at(0).metadataVector);
        for (auto& entry : dihedrals)
        {
            auto metadata = entry.metadataVector[conformerIndex];
            setDihedralAngle(entry, randomAngleWithMetadata(conformerIndex, metadata));
        }
    }
}

void cds::setRandomShapeUsingMetadata(MetadataDistribution randomMetadata, AngleDistribution randomAngle,
                                      std::vector<ResidueLinkage>& linkages)
{
    for (auto& linkage : linkages)
    {
        setRandomShapeUsingMetadata(randomMetadata, randomAngle, linkage.rotamerType, linkage.rotatableDihedrals);
    }
}

void cds::setSpecificShapeUsingMetadata(std::vector<RotatableDihedral>& dihedrals, size_t shapeNumber)
{
    for (auto& entry : dihedrals)
    {
        auto& metadata = entry.metadataVector[shapeNumber];
        setDihedralAngle(entry, {metadata.default_angle_value_, metadata.default_angle_value_, shapeNumber});
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
