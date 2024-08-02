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

bool cds::setSpecificShape(RotatableDihedral& dihedral, const DihedralAngleDataVector& metadataVector,
                           std::string dihedralName, std::string selectedRotamer)
{
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

std::vector<cds::AngleWithMetadata> cds::currentShape(const std::vector<RotatableDihedral>& dihedrals,
                                                      const std::vector<DihedralAngleDataVector>& metadata)
{
    std::vector<AngleWithMetadata> result;
    result.reserve(dihedrals.size());

    for (size_t n = 0; n < dihedrals.size(); n++)
    {
        auto& dihedral        = dihedrals[n];
        auto& metadataIndex   = dihedral.currentMetadataIndex;
        auto& currentMetadata = metadata[n][metadataIndex];
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
        result.push_back(currentShape(a.rotatableDihedrals, a.dihedralMetadata));
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

void cds::setDefaultShapeUsingMetadata(std::vector<cds::RotatableDihedral>& dihedrals,
                                       const std::vector<DihedralAngleDataVector>& metadata)
{
    for (size_t n = 0; n < dihedrals.size(); n++)
    {
        double defaultAngle = metadata[n][0].default_angle_value_;
        setDihedralAngle(dihedrals[n], {defaultAngle, defaultAngle, 0});
    }
}

void cds::setRandomShapeUsingMetadata(MetadataDistribution randomMetadata, AngleDistribution randomAngle,
                                      gmml::MolecularMetadata::GLYCAM::RotamerType rotamerType,
                                      std::vector<RotatableDihedral>& dihedrals,
                                      const std::vector<DihedralAngleDataVector>& metadata)
{
    auto randomAngleWithMetadata = [&](size_t index, const DihedralAngleDataVector& metadataVector)
    {
        auto& metadata = metadataVector[index];
        double angle   = randomAngle(metadata);
        return AngleWithMetadata {angle, metadata.default_angle_value_, index};
    };

    if (rotamerType == gmml::MolecularMetadata::GLYCAM::RotamerType::permutation)
    {
        for (size_t n = 0; n < dihedrals.size(); n++)
        {
            size_t index = randomMetadata(metadata[n]);
            setDihedralAngle(dihedrals[n], randomAngleWithMetadata(index, metadata[n]));
        }
    }
    else if (rotamerType == gmml::MolecularMetadata::GLYCAM::RotamerType::conformer)
    {
        size_t conformerIndex = randomMetadata(metadata[0]);
        for (size_t n = 0; n < dihedrals.size(); n++)
        {
            setDihedralAngle(dihedrals[n], randomAngleWithMetadata(conformerIndex, metadata[n]));
        }
    }
}

void cds::setRandomShapeUsingMetadata(MetadataDistribution randomMetadata, AngleDistribution randomAngle,
                                      std::vector<ResidueLinkage>& linkages)
{
    for (auto& linkage : linkages)
    {
        setRandomShapeUsingMetadata(randomMetadata, randomAngle, linkage.rotamerType, linkage.rotatableDihedrals,
                                    linkage.dihedralMetadata);
    }
}

void cds::setSpecificShapeUsingMetadata(std::vector<RotatableDihedral>& dihedrals,
                                        const std::vector<DihedralAngleDataVector>& metadata, size_t shapeNumber)
{
    for (size_t n = 0; n < dihedrals.size(); n++)
    {
        double defaultAngle = metadata[n][shapeNumber].default_angle_value_;
        setDihedralAngle(dihedrals[n], {defaultAngle, defaultAngle, shapeNumber});
    }
}

void cds::setSpecificShape(std::vector<RotatableDihedral>& dihedrals,
                           const std::vector<DihedralAngleDataVector>& metadata, std::string dihedralName,
                           std::string selectedRotamer)
{
    for (size_t n = 0; n < dihedrals.size(); n++)
    {
        // This will call RotatableDihedrals that don't have dihedralName (phi,psi), and nothing will happen. Hmmm.
        if (setSpecificShape(dihedrals[n], metadata[n], dihedralName, selectedRotamer))
        {
            return; // Return once you manage to set a shape.
        }
    }
    std::string errorMessage = "Did not set " + dihedralName + " to " + selectedRotamer +
                               " as requested in ResidueLinkage::SetSpecificShape()";
    gmml::log(__LINE__, __FILE__, gmml::ERR, errorMessage);
    throw std::runtime_error(errorMessage);
}
