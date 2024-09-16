#include "includes/CentralDataStructure/Editors/superimposition.hpp"

#include "includes/CentralDataStructure/Geometry/types.hpp"
#include "includes/CodeUtils/references.hpp"
#include <eigen3/Eigen/Geometry>

namespace
{
    void GenerateMatrixFromCoordinates(const std::vector<cds::CoordinateReference>& coordinates,
                                       Eigen::Matrix3Xd& matrix)
    {
        int col = 0; // Column index for matrix
        for (auto& coordinate : coordinates)
        {
            auto coord = coordinate.get();
            for (size_t n = 0; n < 3; n++)
            {
                matrix(n, col) = coord.nth(n);
            }
            ++col;
        }
    }

    void ReplaceCoordinatesFromMatrix(std::vector<cds::CoordinateReference>& coordinates,
                                      const Eigen::Matrix3Xd& matrix)
    {
        int col = 0; // Column index for matrix
        for (auto& coordinate : coordinates)
        {
            coordinate.set({matrix(0, col), matrix(1, col), matrix(2, col)});
            ++col;
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

void cds::Superimpose(std::vector<CoordinateReference>& moving, const std::vector<CoordinateReference>& target)
{
    Eigen::Matrix3Xd movingMatrix(3, moving.size()), targetMatrix(3, target.size());

    // Create Matrices containing co-ordinates of moving and target
    GenerateMatrixFromCoordinates(moving, movingMatrix);
    GenerateMatrixFromCoordinates(target, targetMatrix);

    // Figure out how to move assembly moving onto target
    Eigen::Affine3d Affine = Find3DAffineTransform(movingMatrix, targetMatrix);

    // Create a matirx containing the moved co-ordinates of assembly moving.
    Eigen::Matrix3Xd movedMatrix = (Affine * movingMatrix);

    // Replace co-ordinates of moving with moved co-ordinates
    ReplaceCoordinatesFromMatrix(moving, movedMatrix);
}

void cds::Superimpose(std::vector<CoordinateReference>& moving, const std::vector<CoordinateReference>& target,
                      std::vector<CoordinateReference>& alsoMoving)
{
    Eigen::Matrix3Xd movingMatrix(3, moving.size()), targetMatrix(3, target.size());
    Eigen::Matrix3Xd alsoMovingMatrix(3, alsoMoving.size()); // separate from above line for clarity

    // Create a matrices containing co-ordinates of assembly moving and target
    GenerateMatrixFromCoordinates(moving, movingMatrix);
    GenerateMatrixFromCoordinates(target, targetMatrix);
    GenerateMatrixFromCoordinates(alsoMoving, alsoMovingMatrix);

    // Figure out how to move assembly moving onto target
    Eigen::Affine3d Affine = Find3DAffineTransform(movingMatrix, targetMatrix);

    // Create a matrix containing the moved co-ordinates of assembly moving and alsoMoving.
    Eigen::Matrix3Xd movedMatrix     = (Affine * movingMatrix);
    Eigen::Matrix3Xd alsoMovedMatrix = (Affine * alsoMovingMatrix);

    // Replace co-ordinates of moving with moved co-ordinates
    ReplaceCoordinatesFromMatrix(moving, movedMatrix);
    ReplaceCoordinatesFromMatrix(alsoMoving, alsoMovedMatrix);
}
