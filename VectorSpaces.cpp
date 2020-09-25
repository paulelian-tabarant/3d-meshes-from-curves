#include "VectorSpaces.h"
#include <iostream>

#include <Eigen/Dense>

Eigen::RowVector3d VectorSpaces::planarTo3DRow(double x, double y, Eigen::Matrix<double, 3, 4> plane) {
    Eigen::Vector3d
        O = plane.col(0),
        U = plane.col(1),
        V = plane.col(2);

        return (O + U * x + V * y).transpose();
}

double VectorSpaces::rescale(double normVal, double minVal, double maxVal) {
    return minVal * (1 - normVal) + maxVal * normVal;
}

double VectorSpaces::normalize(double val, double minVal, double maxVal) {
    return (val - minVal) / (maxVal - minVal);
}

Eigen::Vector2d VectorSpaces::orthoProjection(Eigen::Vector3d point, Eigen::Matrix<double, 3, 4> plane) {
    Eigen::Vector3d
    O = plane.col(0),
    X = plane.col(1),
    Y = plane.col(2);
    double x = (point - O).dot(X),
           y = (point - O).dot(Y);
    return Eigen::Vector2d(x, y);
}

Eigen::Vector3d VectorSpaces::barycentric2D(Eigen::Vector2d p, Eigen::Vector2d a, Eigen::Vector2d b, Eigen::Vector2d c) {
    // Transcribed from Christer Ericson's Real-Time Collision Detection book
    Eigen::Vector2d v0 = b - a, v1 = c - a, v2 = p - a;
    double d00 = v0.squaredNorm(),
           d01 = v0.dot(v1),
           d11 = v1.squaredNorm(),
           d20 = v2.dot(v0),
           d21 = v2.dot(v1);
    double denom = d00 * d11 - d01 * d01;
    double alpha = (d11 * d20 - d01 * d21) / denom,
           beta  = (d00 * d21 - d01 * d20) / denom,
           gamma = 1 - alpha - beta;
    return Eigen::Vector3d(alpha, beta, gamma);
}

std::vector<std::vector<Eigen::Vector3d>> VectorSpaces::curvesProjection3d(
        std::vector<std::vector<Eigen::Vector3d>>& C, 
        Eigen::Matrix<double, 3, 4> plane) {

    std::vector<std::vector<Eigen::Vector3d>> C_p(C.size());
    for (size_t ic = 0; ic < C.size(); ic++) {
        for (size_t ip = 0; ip < C[ic].size(); ip++) {
            // Compute 3d projection of each coordinate of the source curves on the plane
            Eigen::Vector2d pos2d = VectorSpaces::orthoProjection(C[ic][ip], plane);
            Eigen::RowVector3d pos3d_p = VectorSpaces::planarTo3DRow(pos2d(0), pos2d(1), plane);
            C_p[ic].push_back(pos3d_p);
        }
    }
    return C_p;
}

Eigen::Vector3d VectorSpaces::quaternionToNormal(Eigen::Quaterniond q) {
    // Quaternions are computed from the y axis on Unity
    Eigen::Vector3d k = Eigen::Vector3d(0, 1, 0); 
    return q._transformVector(k);   
}

Eigen::RowVector3d VectorSpaces::ellipseToUnitSphere(Eigen::RowVector3d pos, Eigen::RowVector3d c, Eigen::Matrix3d rot, Eigen::RowVector3d dim) {
    Eigen::RowVector3d newPos;

    // Convert coordinates in the ellipse local space
    newPos = (pos - c).transpose();
    newPos *= rot;
    newPos(0) /= dim(0);
    newPos(1) /= dim(1);
    newPos(2) /= dim(2);

    return newPos;
}


Eigen::RowVector3d VectorSpaces::projectToOctahedronNearest(Eigen::RowVector3d pos, Eigen::RowVector3d center, Eigen::Matrix3d rot, Eigen::RowVector3d dim) {
    Eigen::RowVector3d newPos;

    // Convert coordinates back from the ellipse local space
    double a =  dim(0) > EPSILON_DOUBLE ? 1./dim(0) : MAX_DOUBLE;
    double b =  dim(1) > EPSILON_DOUBLE ? 1./dim(1) : MAX_DOUBLE;
    double c =  dim(2) > EPSILON_DOUBLE ? 1./dim(2) : MAX_DOUBLE;
    newPos = (pos-center).transpose();
    newPos *= rot;

    double Alpha = (1. - ( fabs(newPos(0))*a + fabs(newPos(1))*b + fabs(newPos(2))*c ) ) / (a*a+b*b+c*c) ;
    newPos(0) = newPos(0) > 0 ? newPos(0) + Alpha * a : newPos(0) - Alpha * a;
    newPos(1) = newPos(1) > 0 ? newPos(1) + Alpha * b : newPos(1) - Alpha * b;
    newPos(2) = newPos(2) > 0 ? newPos(2) + Alpha * c : newPos(2) - Alpha * c;

//    printf("dans l'octaedre si ce nombre vaut 0: %lf\n", fabs(newPos(0))*a + fabs(newPos(1))*b + fabs(newPos(2))*c - 1.);
    newPos *= rot.transpose();
    newPos += center.transpose();

    return newPos;
}

Eigen::RowVector3d VectorSpaces::projectToUnitOctahedronNearest(Eigen::RowVector3d pos, Eigen::RowVector3d center, Eigen::Matrix3d rot, Eigen::RowVector3d dim) {
    Eigen::RowVector3d newPos;

    // Convert coordinates back from the ellipse local space
    double a = dim(0) > EPSILON_DOUBLE ? 1./dim(0) : MAX_DOUBLE;
    double b = dim(1) > EPSILON_DOUBLE ? 1./dim(1) : MAX_DOUBLE;
    double c = dim(2) > EPSILON_DOUBLE ? 1./dim(2) : MAX_DOUBLE;
    newPos = (pos-center).transpose();
    newPos *= rot;

    double Alpha = (1. - ( fabs(newPos(0))*a + fabs(newPos(1))*b + fabs(newPos(2))*c ) ) / (a*a+b*b+c*c) ;
    newPos(0) = newPos(0) > 0 ? newPos(0) + Alpha * a : newPos(0) - Alpha * a;
    newPos(1) = newPos(1) > 0 ? newPos(1) + Alpha * b : newPos(1) - Alpha * b;
    newPos(2) = newPos(2) > 0 ? newPos(2) + Alpha * c : newPos(2) - Alpha * c;

    newPos(0) *= a;
    newPos(1) *= b;
    newPos(2) *= c;

    return newPos;
}

Eigen::Vector2d VectorSpaces::trianglePlanarCoords(Eigen::RowVector3d pos, Eigen::RowVector3d t0, Eigen::RowVector3d t1, Eigen::RowVector3d t2) {
    Eigen::RowVector3d t0p = pos - t0,
                       t0t1 = t1 - t0,
                       t0t2 = t2 - t0;
    
    double x = t0p.dot(t0t1) / t0t1.squaredNorm(),
           y = t0p.dot(t0t2) / t0t2.squaredNorm();

    return Eigen::Vector2d(x, y);
}
