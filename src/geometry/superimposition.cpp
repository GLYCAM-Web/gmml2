#include "include/geometry/superimposition.hpp"

#include "include/geometry/geometryTypes.hpp"

#include <eigen3/Eigen/Geometry>
#include <stdexcept>
#include <vector>

namespace gmml
{
    namespace
    {
        double sign(double a) { return a < 0 ? -1.0 : 1.0; }
    } // namespace

    Eigen::Affine3d affineTransform(const Eigen::Matrix3Xd& target, const Eigen::Matrix3Xd& initial)
    {
        Eigen::Matrix3Xd in = initial;
        Eigen::Matrix3Xd out = target;
        // Default output
        Eigen::Affine3d A;
        A.linear() = Eigen::Matrix3d::Identity(3, 3);
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
            dist_in += (in.col(col + 1) - in.col(col)).norm();
            dist_out += (out.col(col + 1) - out.col(col)).norm();
        }
        if (dist_in <= 0 || dist_out <= 0)
        {
            return A;
        }
        double scale = dist_out / dist_in;
        out /= scale;
        // Find the centroids then shift to the origin
        Eigen::Vector3d in_ctr = Eigen::Vector3d::Zero();
        Eigen::Vector3d out_ctr = Eigen::Vector3d::Zero();
        for (int col = 0; col < in.cols(); col++)
        {
            in_ctr += in.col(col);
            out_ctr += out.col(col);
        }
        in_ctr /= in.cols();
        out_ctr /= out.cols();
        for (int col = 0; col < in.cols(); col++)
        {
            in.col(col) -= in_ctr;
            out.col(col) -= out_ctr;
        }
        // SVD
        Eigen::MatrixXd Cov = in * out.transpose();
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(Cov, Eigen::ComputeThinU | Eigen::ComputeThinV);
        // Find the rotation
        double d = sign((svd.matrixV() * svd.matrixU().transpose()).determinant());
        Eigen::Matrix3d I = Eigen::Matrix3d::Identity(3, 3);
        I(2, 2) = d;
        Eigen::Matrix3d R = svd.matrixV() * I * svd.matrixU().transpose();
        // The final transform
        A.linear() = scale * R;
        A.translation() = scale * (out_ctr - R * in_ctr);
        return A;
    }

    Eigen::Matrix3Xd generateMatrix(const std::vector<Coordinate>& coordinates)
    {
        Eigen::Matrix3Xd result(3, coordinates.size());
        for (size_t k = 0; k < coordinates.size(); k++)
        {
            Coordinate coord = coordinates[k];
            for (size_t n = 0; n < 3; n++)
            {
                result(n, k) = coord.nth(n);
            }
        }
        return result;
    }

    std::vector<Coordinate> matrixCoordinates(const Eigen::Matrix3Xd& matrix)
    {
        std::vector<Coordinate> result;
        Eigen::Index columns = matrix.cols();
        result.reserve(columns);
        for (Eigen::Index k = 0; k < columns; k++)
        {
            result.push_back({matrix(0, k), matrix(1, k), matrix(2, k)});
        }
        return result;
    }
} // namespace gmml
