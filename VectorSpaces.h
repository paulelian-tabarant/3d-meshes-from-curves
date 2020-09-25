#pragma once

#include <Eigen/Dense>
#include <vector>

/**
 * @brief Miscellaneous conversion functions between vector spaces
 * 
 */
namespace VectorSpaces {
    const double EPSILON_DOUBLE { std::numeric_limits<double>::epsilon() };
    const double MAX_DOUBLE     { std::numeric_limits<double>::max()     };

    /**
     * @brief Convert (X,Y) planar coordinates into the 3D space
     * 
     * @param x 
     * @param y 
     * @param plane The containing plane in a 3x4 matrix with (O,X,Y,Z) as columns
     * @return Eigen::RowVector3d 
     */
    Eigen::RowVector3d planarTo3DRow(double x, double y, Eigen::Matrix<double, 3, 4> plane); 

    /**
     * @brief Linear interpolation of a normalized value
     * between minVal and maxVal
     * 
     * @param normVal 
     * @param minVal 
     * @param maxVal 
     * @return double 
     */
    double rescale(double normVal, double minVal, double maxVal);

    /**
     * @brief Bijection from [vMin, vMax]
     * to [0, 1] unit interval
     * 
     * @param val 
     * @param vMin 
     * @param vMax 
     * @return double 
     */
    double normalize(double val, double vMin, double vMax);

    /**
     * @brief Compute an orthogonal projection of a 3D points
     * on a plane defined in the same coordinates system
     * 
     * @param point 
     * @param plane 
     * @return Eigen::Vector2d 
     */
    Eigen::Vector2d orthoProjection(Eigen::Vector3d point, Eigen::Matrix<double, 3, 4> plane);

    /**
     * @brief Compute barycentric coordinates of a point
     * inside a 2D triangle
     * 
     * @param p 
     * @param a 
     * @param b 
     * @param c 
     * @return Eigen::Vector3d 
     */
    Eigen::Vector3d barycentric2D(Eigen::Vector2d p, Eigen::Vector2d a, Eigen::Vector2d b, Eigen::Vector2d c);

    /**
     * @brief Get a set of curves as a vector of vectors on the given plane
     */
    std::vector<std::vector<Eigen::Vector3d>> curvesProjection3d(std::vector<std::vector<Eigen::Vector3d>>& C,
                            Eigen::Matrix<double, 3, 4> plane);

    /**
     * @brief Compute the normal value from a quaternion, thanks to a 
     * reference vector
     */
    Eigen::Vector3d quaternionToNormal(Eigen::Quaterniond q);

    /**
     * @brief Convert 3D positions from an elliptic frame to a unit sphere
     * coordinate system, based on given ellipse dimensions
     */
    Eigen::RowVector3d ellipseToUnitSphere(Eigen::RowVector3d pos, Eigen::RowVector3d c, Eigen::Matrix3d rot, Eigen::RowVector3d dim);

    /**
     * @brief Project 3D position to an octahedron defined by its position (center), its orientation
     * (rot) and its dimensions (dim)
     */
    Eigen::RowVector3d projectToOctahedronNearest(Eigen::RowVector3d pos, Eigen::RowVector3d c, Eigen::Matrix3d rot, Eigen::RowVector3d dim);

    /**
     * @brief Project 3D position to an octahedron defined by its position (center),
     * its orientation (rot) and its dimensions (dim) as normalized coordinates :
     * position is converted in a "unit octahedron" coordinates system
     */
    Eigen::RowVector3d projectToUnitOctahedronNearest(Eigen::RowVector3d pos, Eigen::RowVector3d c, Eigen::Matrix3d rot, Eigen::RowVector3d dim);

    /**
     * @brief Project 3D position inside a 3D triangle, resulting in 2D coordinates
     * defined by the two vectors t0t1 and t0t2
     */
    Eigen::Vector2d trianglePlanarCoords(Eigen::RowVector3d pos, Eigen::RowVector3d t0, Eigen::RowVector3d t1, Eigen::RowVector3d t2);
}
