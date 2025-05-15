#include "includes/CentralDataStructure/Editors/superimposition.hpp"

#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
#include <eigen3/Eigen/Geometry>

namespace
{
    void ReplaceCoordinatesFromMatrix(std::vector<cds::Coordinate>& coordinates, const Eigen::Matrix3Xd& matrix)
    {
        for (size_t k = 0; k < coordinates.size(); k++)
        {
            coordinates[k] = {matrix(0, k), matrix(1, k), matrix(2, k)};
        }
    }

    double sign(double a)
    {
        return a < 0 ? -1.0 : 1.0;
    }

    Eigen::Affine3d Find3DAffineTransform(Eigen::Matrix3Xd in, Eigen::Matrix3Xd out)
    {
        // Default output
        Eigen::Affine3d A;
        A.linear()      = Eigen::Matrix3d::Identity(3, 3);
        A.translation() = Eigen::Vector3d::Zero();
        if (in.cols() != out.cols())
        {
            throw "Find3DAffineTransform(): input data mis-match";
        }
        // First find the scale, by finding the ratio of sums of some distances,
        // then bring the datasets to the same scale.
        double dist_in = 0, dist_out = 0;
        for (int col = 0; col < in.cols() - 1; col++)
        {
            dist_in  += (in.col(col + 1) - in.col(col)).norm();
            dist_out += (out.col(col + 1) - out.col(col)).norm();
        }
        if (dist_in <= 0 || dist_out <= 0)
        {
            return A;
        }
        double scale            = dist_out / dist_in;
        out                     /= scale;
        // Find the centroids then shift to the origin
        Eigen::Vector3d in_ctr  = Eigen::Vector3d::Zero();
        Eigen::Vector3d out_ctr = Eigen::Vector3d::Zero();
        for (int col = 0; col < in.cols(); col++)
        {
            in_ctr  += in.col(col);
            out_ctr += out.col(col);
        }
        in_ctr  /= in.cols();
        out_ctr /= out.cols();
        for (int col = 0; col < in.cols(); col++)
        {
            in.col(col)  -= in_ctr;
            out.col(col) -= out_ctr;
        }
        // SVD
        Eigen::MatrixXd Cov = in * out.transpose();
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(Cov, Eigen::ComputeThinU | Eigen::ComputeThinV);
        // Find the rotation
        double d          = sign((svd.matrixV() * svd.matrixU().transpose()).determinant());
        Eigen::Matrix3d I = Eigen::Matrix3d::Identity(3, 3);
        I(2, 2)           = d;
        Eigen::Matrix3d R = svd.matrixV() * I * svd.matrixU().transpose();
        // The final transform
        A.linear()        = scale * R;
        A.translation()   = scale * (out_ctr - R * in_ctr);
        return A;
    }
} // namespace

Eigen::Matrix3Xd cds::generateMatrix(const std::vector<cds::Coordinate>& coordinates)
{
    Eigen::Matrix3Xd result(3, coordinates.size());
    for (size_t k = 0; k < coordinates.size(); k++)
    {
        cds::Coordinate coord = coordinates[k];
        for (size_t n = 0; n < 3; n++)
        {
            result(n, k) = coord.nth(n);
        }
    }
    return result;
}

std::vector<cds::Coordinate> cds::matrixCoordinates(const Eigen::Matrix3Xd& matrix)
{
    std::vector<cds::Coordinate> result;
    Eigen::Index columns = matrix.cols();
    result.reserve(columns);
    for (Eigen::Index k = 0; k < columns; k++)
    {
        result.push_back({matrix(0, k), matrix(1, k), matrix(2, k)});
    }
    return result;
}

cds::AffineTransform cds::affineTransform(const std::vector<cds::Coordinate>& target,
                                          const std::vector<cds::Coordinate>& moving)
{
    Eigen::Matrix3Xd targetMatrix = generateMatrix(target);
    Eigen::Matrix3Xd movingMatrix = generateMatrix(moving);
    Eigen::Affine3d transform     = Find3DAffineTransform(movingMatrix, targetMatrix);
    return {movingMatrix, targetMatrix, transform};
}

void cds::Superimpose(std::vector<Coordinate>& moving, const std::vector<Coordinate>& target)
{
    Eigen::Matrix3Xd targetMatrix = generateMatrix(target);
    Eigen::Matrix3Xd movingMatrix = generateMatrix(moving);

    // Figure out how to move assembly moving onto target
    Eigen::Affine3d Affine = Find3DAffineTransform(movingMatrix, targetMatrix);

    // Create a matirx containing the moved co-ordinates of assembly moving.
    Eigen::Matrix3Xd movedMatrix = (Affine * movingMatrix);

    // Replace co-ordinates of moving with moved co-ordinates
    ReplaceCoordinatesFromMatrix(moving, movedMatrix);
}

void cds::Superimpose(std::vector<Coordinate>& moving, const std::vector<Coordinate>& target,
                      std::vector<Coordinate>& alsoMoving)
{
    Eigen::Matrix3Xd targetMatrix     = generateMatrix(target);
    Eigen::Matrix3Xd movingMatrix     = generateMatrix(moving);
    Eigen::Matrix3Xd alsoMovingMatrix = generateMatrix(alsoMoving);

    // Figure out how to move assembly moving onto target
    Eigen::Affine3d Affine = Find3DAffineTransform(movingMatrix, targetMatrix);

    // Create a matrix containing the moved co-ordinates of assembly moving and alsoMoving.
    Eigen::Matrix3Xd movedMatrix     = (Affine * movingMatrix);
    Eigen::Matrix3Xd alsoMovedMatrix = (Affine * alsoMovingMatrix);

    // Replace co-ordinates of moving with moved co-ordinates
    ReplaceCoordinatesFromMatrix(moving, movedMatrix);
    ReplaceCoordinatesFromMatrix(alsoMoving, alsoMovedMatrix);
}
