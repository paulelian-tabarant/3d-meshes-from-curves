#define _USE_MATH_DEFINES

#include "ApproxSurface.h"
#include "VectorSpaces.h"

#include <cmath>
#include <vector>
#include <map>
#include <algorithm>
#include <time.h>

#include <Eigen/SVD>
#include <Eigen/SparseQR>
#include <Eigen/OrderingMethods>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Geometry>

#include <iostream>

using namespace Eigen;

Eigen::Vector3d ApproxSurface::centroid(const std::vector<Coord> &coords)
{
    Vector3d c(Vector3d::Zero());
    int cpt = 0;
    for (const Coord &coord : coords) {
        c += coord.pos;
        cpt++;
    }
    c /= (double)cpt;
    return c;
}

void ApproxSurface::fittingPlane() 
{
    #ifndef TEST_FITTING_PLANE

        std::vector<Coord> coords; 
        curvesSrc->getCoords(coords, computeAsTape);

        Vector3d c = centroid(coords);
        P.col(0) = c;

        MatrixXd M(3, coords.size());

        // Store each point position 
        int cpt = 0;
        for (int ic = 0; ic < positions.size(); ic++) {
            for (int ip = 0; ip < positions[ic].size(); ip++) {
                M.col(cpt) = positions[ic][ip] - c;
                cpt++;
            }
        }

        // SVD decomposition of the M matrix using BDC method
        // (Bidiagonal Divide & Conquer)
        Eigen::BDCSVD<MatrixXd> svd(M, ComputeThinU);

        Matrix<double, 3, 3> U = svd.matrixU();
        // Last column gives the fitting plane normal vector
        P.col(3) = U.col(2);

        // If fitting plane directions should be updated
        if (P.isApprox(Matrix<double,3,4>::Zero()) || !constantDirs) {
            P.col(1) = U.col(0);
            P.col(2) = U.col(1);
        }

    #else

        std::vector<Coord> coords;
        curvesSrc->getCoords(coords, computeAsTape);

        Vector3d nAvg = Vector3d::Zero(), c = Vector3d::Zero();
        int nPos = 0;
        for (Coord &coord : coords) {
            c += coord.pos;
            nAvg += coord.normal;
            nPos++;
        }
        c /= nPos;
        nAvg.normalize();
        P.col(0) = c;
        P.col(3) = nAvg;
    
        // If plane directions should be updated
        if (computeFittingPlane || P.col(1).isApprox(Vector3d::Zero()) || !constantDirs) {
            // Algorithm from Hughes and Moller (1999)
            Vector3d u, v;
            if (fabs(nAvg(0)) > fabs(nAvg(2))) u = Vector3d(-nAvg(1), nAvg(0), 0.0);
            else u = Vector3d(0.0, -nAvg(2), nAvg(1));
            u.normalize();
            v = nAvg.cross(u);

            MatrixXd M(2, nPos);
            int cpt = 0;
            for (Coord &coord : coords) {
                M(0, cpt) = (coord.pos - c).dot(u);
                M(1, cpt) = (coord.pos - c).dot(v);
                cpt++;
            }
            BDCSVD<MatrixXd> svd(M, ComputeThinU);
            Matrix<double, 2, 2> U = svd.matrixU();
            P.col(1) = U(0,0) * u + U(1,0) * v;
            P.col(2) = U(0,1) * u + U(1,1) * v;

            computeFittingPlane = false;
        }

    #endif

    // Store P data for further computations
    Vector4d bounds = curvesSrc->getProjectionBounds(computeAsTape, P);
    xMin = bounds(0);
    yMin = bounds(1);
    xMax = bounds(2);
    yMax = bounds(3);
}
 
void ApproxSurface::setConstantFittingPlane(const bool &constant) {
    // When fitting plane is firstly set to constant, the former must be
    // computed once at least
    if (constant && !constantDirs) {
        fittingPlane();
        computeFittingPlane = true;
    }

    constantDirs = constant;
}

std::string ApproxSurface::triangleToMapStr(unsigned int i0, unsigned int i1, unsigned int i2) {
    if (i0 > i1)
        std::swap(i0, i1);
    if (i0 > i2)
        std::swap(i0, i2);
    if (i1 > i2)
        std::swap(i1, i2);
    
    std::ostringstream tStr;
    tStr << i0 << " " << i1 << " " << i2;
    return tStr.str();
}

void ApproxSurface::initTriangleIndicesData() 
{
    rowIndexOfTriangle.clear();

    // Filling map for Triangle -> index associations
    for (unsigned int i = 0; i < T.rows(); i++) {
        RowVector3i t = T.row(i);
        std::string key = triangleToMapStr(t(0), t(1), t(2));
        rowIndexOfTriangle[key] = i;
    }
}

void ApproxSurface::planarMesh(const unsigned int &M, const unsigned int &N) 
{
    V_plane = MatrixXd((M+2) * (N+1) + (M+1) / 2, 3);
    borderV.clear();

    // 6 T for N vertices locally
    T = MatrixXi((M/2 + 1) * (N*4 + 2), 3);

    // Getting bounds for optimal triangulation on the P
    // Bounds = {minX, minY, maxX, maxY}T
    unsigned int vIndex = 0;
    double y0 = -(1.f/M) / 2;
    // Vertices positions
    for (unsigned int i = 0; i < M+2; i++) {
        double y = y0 + i / (double)M;
        y = VectorSpaces::rescale(y, yMin, yMax);

        // Start current line from 0 or 1 half-resolution step
        bool even = i % 2 == 0;
        double x0 = even ? 0 : -(1.f/N) / 2;
        unsigned int Nsup = even ? N+1 : N+2;
        for (unsigned int j = 0; j < Nsup; j++) {
            double x = x0 + j / (double)N;
            x = VectorSpaces::rescale(x, xMin, xMax);
            V_plane.row(vIndex) = VectorSpaces::planarTo3DRow(x, y, P);

            if (i > 0 && i < M+1 && (j == 0 || j == Nsup-1))
                borderV.push_back(vIndex);

            vIndex++;
        }
    }

    // Triangles : we draw a tiling of hexagons based on the previous vertices
    t1[0] =    N+1, t1[1] =    N+2;
    t2[0] =    N+2, t2[1] =      1;
    t3[0] =      1, t3[1] = -(N+1);
    t4[0] = -(N+1), t4[1] = -(N+2);
    t5[0] = -(N+2), t5[1] =     -1;
    t6[0] =     -1, t6[1] =    N+1;
    // vi : store the first hexagon center index of the current iteration
    h0 = N+2;
    hstep = 2*(N+1) + 1;
    unsigned int tIndex = 0, vi = h0;
    for (unsigned int i = 1; i < M+1; i += 2) {
        // We draw a complete hexagon at first iteration of each row
        T.row(tIndex)     = RowVector3i(vi, vi+t5[0], vi+t5[1]);
        T.row(tIndex + 1) = RowVector3i(vi, vi+t6[0], vi+t6[1]);
        tIndex += 2;
        // Next, we only draw the 4 right-most T of the hexagon
        for (unsigned int j = 0; j < N; j++) {
            unsigned int vj = vi + j;
            T.row(tIndex)     = RowVector3i(vj, vj+t1[0], vj+t1[1]);
            T.row(tIndex + 1) = RowVector3i(vj, vj+t2[0], vj+t2[1]);
            T.row(tIndex + 2) = RowVector3i(vj, vj+t3[0], vj+t3[1]);
            T.row(tIndex + 3) = RowVector3i(vj, vj+t4[0], vj+t4[1]);
            tIndex += 4;
        }
        vi += hstep;
    }

    // Store resolution into local attributes
    mPlane = M;
    nPlane = N;
    // Store hexagon angles into local attributes
    Vector2d c = VectorSpaces::orthoProjection(V_plane.row(h0), P),
             p0 = VectorSpaces::orthoProjection(V_plane.row(h0 + t3[0]), P),
             p1 = VectorSpaces::orthoProjection(V_plane.row(h0 + t4[0]), P);
    Vector2d cp0 = p0 - c, cp1 = p1 - c;
    theta1 = acos(cp0.dot(cp1) / (cp0.norm() * cp1.norm()));
    theta2 = M_PI - 2 * theta1;

    // Init map for retrieving triangle row index from a triplet of indexes
    initTriangleIndicesData();
}

void ApproxSurface::addBorderTriangles() 
{
    // Should only be called if there is at least 1 triangle to add on each side
    if (borderV.size() < 6) return;

    unsigned int iBorder = T.rows();
    T.conservativeResize(T.rows() + borderV.size() / 2 - 1, 3);
    unsigned int tIndex = 0;
    for (size_t i = 0; i < borderV.size() - 4; i += 4) {
        int i0 = borderV[i], i1 = borderV[i+2], i2 = borderV[i+4];
        T.row(iBorder + tIndex) = RowVector3i(i0, i2, i1);
        i0 = borderV[i+1], i1 = borderV[i+3], i2 = borderV[i+5];
        T.row(iBorder + tIndex + 1) = RowVector3i(i0, i2, i1);
        tIndex += 2;
    }
}

void ApproxSurface::planarMesh(const unsigned int &N) 
{
    // We need a fitting plane before computing planar triangulation
    fittingPlane();

    double w = fabs(xMax - xMin),
           l = fabs(yMax - yMin);
    int M = ceil((l/w) * N);
    if (M < 0 || M > 1e6) {
        std::cerr << "Error: bad resolution value for height (M), w = " << w << std::endl;
		M = 2;
	}
    if (M % 2 == 0) 
        M--;

    planarMesh(M, N);
}

unsigned int ApproxSurface::planeBoundingTriangle(const double &x, const double &y) 
{
    double u = VectorSpaces::normalize(x, xMin, xMax),
           v = VectorSpaces::normalize(y, yMin, yMax);
    double yi = v * mPlane + 0.5;
    double xj = u * nPlane + 0.5;
    // j = 0 for first hexagon on the line
    unsigned int j = round(xj);
    // Find i such as the (i+1)th hexagon line is the nearest from the point
    unsigned int i = (int)yi;
    if (i % 2 == 0) 
        i++;
    i /= 2;
    // Index of the nearest hexagon center in V_plane matrix
    unsigned int hcIndex = h0 + i*hstep + (j - 1);

    // If (x,y) falls outside the fitting plane limits
    if (hcIndex < 0 || hcIndex >= V_plane.rows()) return -1;

    Vector3d hcPos3D = V_plane.row(hcIndex).transpose();
    Vector2d hcPos2D = VectorSpaces::orthoProjection(hcPos3D, P);

    // Find the bounding triangle inside the bounding hexagon
    // y axis is inverted from the triangular mesh point of view
    double angle = atan2(hcPos2D(1) - y, x - hcPos2D(0));
    unsigned int *t;
    if (angle > 0) {
        if (angle < theta1)
            t = t3;
        else if (angle < theta1 + theta2)
            t = t4;
        else
            t = t5;
    }
    else {
        if (angle > -theta1)
            t = t2;
        else if (angle > -(theta1 + theta2))
            t = t1;
        else
            t = t6;
    }

    std::string key = triangleToMapStr(hcIndex, hcIndex + t[0], hcIndex + t[1]);
    return rowIndexOfTriangle[key];
}

void ApproxSurface::sparseLSsolve(const SparseMatrix<double> &A, const MatrixXd &b, MatrixXd& X, const bool &factorize) 
{
    #ifndef LSCG
        static SparseQR<SparseMatrix<double>, COLAMDOrdering<int>> solver;
        // Solve the constraints system with least-squares
        if (factorize) {
            A.makeCompressed();
            solver.compute(A);
        }
        X = solver.solve(b);

    #else

        LeastSquaresConjugateGradient<SparseMatrix<double>> solver;
        solver.compute(A);
        solver.setTolerance(0.0001);
        X = solver.solve(b);

    #endif
}

void ApproxSurface::trianglePlanarCoords(const Vector3i &abc, Vector2d& a2d, Vector2d& b2d, Vector2d& c2d) 
{
    Vector3d a = V_plane.row(abc(0)).transpose(),
             b = V_plane.row(abc(1)).transpose(),
             c = V_plane.row(abc(2)).transpose();
    a2d = VectorSpaces::orthoProjection(a, P);
    b2d = VectorSpaces::orthoProjection(b, P);
    c2d = VectorSpaces::orthoProjection(c, P);
}

void ApproxSurface::cotangentsMatrix(const MatrixXd &vertices, const MatrixXi &triangles, SparseMatrix<double> &cotangents)
{
    cotangents = SparseMatrix<double>(vertices.rows(), vertices.rows());
    std::vector<Triplet<double>> nonZeroCoords;
    nonZeroCoords.reserve(triangles.size() * 3 * 4);

    for (unsigned int it = 0; it < triangles.rows(); it++) {
        RowVector3i triangle = triangles.row(it);
        for (unsigned int iv = 0; iv < 3; iv++) {
            int i0 {triangle[iv]}, i1 {triangle[(iv + 1) % 3]}, i2 {triangle[(iv + 2) % 3]};
            RowVector3d v0 = vertices.row(i0);
            RowVector3d v1 = vertices.row(i1);
            RowVector3d v2 = vertices.row(i2);
            RowVector3d e01 = v1 - v0, e02 = v2 - v0;
            double normProduct = e01.norm() * e02.norm();
            double cos0 = e01.dot(e02) / normProduct;
            // Consider the sine angle to be always positive (no oriented angle in laplacian formula)
            double sin0 = e01.cross(e02).norm() / normProduct;
            double cot0 = cos0 / sin0;

            nonZeroCoords.push_back(Triplet<double>(i1, i2, +cot0));
            nonZeroCoords.push_back(Triplet<double>(i2, i1, +cot0));
            nonZeroCoords.push_back(Triplet<double>(i1, i1, -cot0));
            nonZeroCoords.push_back(Triplet<double>(i2, i2, -cot0));
        }
    }

    cotangents.setFromTriplets(nonZeroCoords.begin(), nonZeroCoords.end());
}

void ApproxSurface::buildFromPlaneBarycentric() 
{

    unsigned int nPos = curvesSrc->totalSize(computeAsTape);
    // Matrices initial sizing & data acquisition
    Ac = SparseMatrix<double>(nPos, V_plane.rows());
    bc = MatrixXd(nPos, 3);

    std::vector<Curve> curves;
    curvesSrc->getCurves(curves, computeAsTape);

    // Non-zero coefficients are stored as a list of triplets <i, j, value>
    std::vector<Triplet<double>> nonZeroCoords;
    nonZeroCoords.reserve(3 * nPos);

    unsigned int pIndex = 0;
    // Compute each projected point barycentric coordinates and build the corresponding matrix
    for (Curve &curve : curves) {
        // Get curve weight for surface construction
        double wi = curve.weight;

        for (Coord &coord : curve.points) {
            Vector3d p = coord.pos;
            double wp = coord.weight;

            // Look for the bounding 2D triangle of the current projected point
            Vector2d xy = VectorSpaces::orthoProjection(p, P);

            unsigned int iT = planeBoundingTriangle(xy(0), xy(1));
            Vector3i bTriangle = T.row(iT).transpose();
            Vector2d a2d, b2d, c2d;
            trianglePlanarCoords(bTriangle, a2d, b2d, c2d);

            // Compute its barycentric coordinates into the 2D P, with curve weight factor
            Vector3d bary = VectorSpaces::barycentric2D(xy, a2d, b2d, c2d);

            // Add corresponding coefficients to the system matrix
            nonZeroCoords.push_back(Triplet<double>(pIndex, bTriangle(0), wi * wp * bary(0)));
            nonZeroCoords.push_back(Triplet<double>(pIndex, bTriangle(1), wi * wp * bary(1)));
            nonZeroCoords.push_back(Triplet<double>(pIndex, bTriangle(2), wi * wp * bary(2)));

            bc.row(pIndex) = wi * wp * p.transpose(); 

            pIndex++;
        }
    }
    // Generate coefficients matrix A and use QR decomposition to solve the problem
    Ac.setFromTriplets(nonZeroCoords.begin(), nonZeroCoords.end());
}

void ApproxSurface::sparseVerticalConcat(const SparseMatrix<double> &M1, const SparseMatrix<double> &M2, SparseMatrix<double>& M) 
{
    if (M1.cols() != M2.cols()) {
        std::cerr << "Concatenated matrices must have the same number of columns."
        << std::endl;
        return;
    }

    M = SparseMatrix<double>(M1.rows() + M2.rows(), M1.cols());
    // Create a list of triplets storing non zero coordinates of the matrix A
    std::vector<Triplet<double>> triplets;
    triplets.reserve(M1.nonZeros() + M2.nonZeros());

    // Fill with M1 part
    for (unsigned int k = 0; k < M1.outerSize(); ++k) {
        for (SparseMatrix<double>::InnerIterator it(M1, k); it; ++it) {
            triplets.push_back(Triplet<double>(it.row(), it.col(), it.value()));
        }
    }
    // Fill with M2 part
    for (unsigned int k = 0; k < M2.outerSize(); ++k) {
        for (SparseMatrix<double>::InnerIterator it(M2, k); it; ++it) {
            triplets.push_back(Triplet<double>(M1.rows() + it.row(), it.col(), it.value()));
        }
    }
    M.setFromTriplets(triplets.begin(), triplets.end());
}

void ApproxSurface::computeFromPlane() 
{
    // Build barycentric constraints in order to use Ac & bc
    buildFromPlaneBarycentric();

    // Compute cotangent coefficients into the 1st sparse matrix
    cotangentsMatrix(V_plane, T, L);

    // Add laplacian constraints to barycentric constraints (with weights)
    sparseVerticalConcat(Ac, L, Acl);
    SparseMatrix<double> Acl_w;
    sparseVerticalConcat(sqrt(1 / smoothing) * Ac, L, Acl_w);

    // Right member of the equation, which is only a zeros Nx3 matrix for the
    // laplacian part
    bcl = MatrixXd(bc.rows() + L.rows(), 3);
    bcl << 	bc,
			MatrixXd::Zero(L.rows(), 3);
    MatrixXd bcl_w(bc.rows() + L.rows(), 3);
    bcl_w << sqrt(1 / smoothing) * bc,
             MatrixXd::Zero(L.rows(), 3);

    // Solve the constraints system with least-squares
	sparseLSsolve(Acl_w, bcl_w, Vcl, true);
}

unsigned int ApproxSurface::i0LineNorth(const unsigned int &nLine) 
{
    if (nLine == 0) return 0;

    return 2 * (nLine * (nLine - 1)) + 1;
}

unsigned int ApproxSurface::i0LineSouth(const unsigned int &nV, const unsigned int &nLine) 
{
    return nV - 2 * (nLine * (nLine + 1)) - 1;
}

unsigned int ApproxSurface::hemisphereOffset(const unsigned int &nLine, const unsigned int &hIndex) 
{
    return hIndex * nLine;
}

void ApproxSurface::octahedronMesh(const unsigned int &nRes) 
{
    RowVector3d up			( 0.0,  0.0,  1.0),
                down		( 0.0,  0.0, -1.0),
                right		( 1.0,  0.0,  0.0),
                forward		( 0.0,  1.0,  0.0),
                left		(-1.0,  0.0,  0.0),
                backward	( 0.0, -1.0,  0.0);

    RowVector3d octaDirs[4] = {right, forward, left, backward}; 

    unsigned int nV = 4 * nRes * nRes + 2;
    unsigned int nT = 8 * nRes * nRes;

    MatrixXd V_unitOctahedron {MatrixXd::Zero(nV, 3)};
    T_octahedron = MatrixXi::Zero(nT, 3);
    unsigned int iV = 0, iT = 0;
    V_unitOctahedron.row(iV++) = up;

    // Up side of the octahedron
    for (unsigned int i = 1; i < nRes; i++) {
        for (unsigned int iDir = 0; iDir < 4; iDir++) {
            RowVector3d lineStart = octaDirs[iDir] - up;
            RowVector3d p0 = up + (i / (double) nRes) * lineStart;
            RowVector3d step = (1 / (double) nRes) * (octaDirs[(iDir+1) % 4]- octaDirs[iDir]);
            for (unsigned int j = 0; j < i; j++) {
                V_unitOctahedron.row(iV++) = p0 + j * step;
            }

        }
    }
    // Down side of the octahedron
    for (unsigned int i = 0; i < nRes; i++) {
        for (unsigned int iDir = 0; iDir < 4; iDir++) {
            RowVector3d lineStart = down - octaDirs[iDir];
            RowVector3d p0 = octaDirs[iDir] + (i / (double) nRes) * lineStart;
            RowVector3d step = (1 / (double) nRes) * (octaDirs[(iDir+1) % 4] - octaDirs[iDir]);
            for (unsigned int j = 0; j < nRes - i; j++) {
                V_unitOctahedron.row(iV++) = p0 + j * step;
            }
        }
    }
    V_unitOctahedron.row(iV++) = down;

    // Triangles
    for (unsigned int i = 1; i < nRes + 1; i++) {
        for (unsigned int iDir = 0; iDir < 4; iDir++) {
            unsigned int nPatterns = i - 1;
            unsigned int i0 = i0LineNorth(i) + hemisphereOffset(i, iDir);
            unsigned int i0_prev = i0LineNorth(i-1) + hemisphereOffset(i-1, iDir);

            unsigned int i1 = i0_prev, i2 = i0;
            for (unsigned int j = 0; j < nPatterns; j++) {
                T_octahedron.row(iT++) = RowVector3i(i1, i2, i2+1);

                if (iDir < 3 || j < nPatterns - 1) 
                    T_octahedron.row(iT++) = RowVector3i(i1, i2+1, i1+1);
                else
                    T_octahedron.row(iT++) = RowVector3i(i1, i2+1, i0LineNorth(i-1));

                i1++;
                i2++;
            }
            if (iDir < 3)
                T_octahedron.row(iT++) = RowVector3i(i1, i2, i2+1); 
            else
                T_octahedron.row(iT++) = RowVector3i(i0LineNorth(i-1), i2, i0LineNorth(i));
        }
    }
    for (unsigned int i = 1; i < nRes + 1; i++) {
        for (unsigned int iDir = 0; iDir < 4; iDir++) {
            unsigned int nPatterns = i - 1;
            unsigned int i0 = i0LineSouth(nV, i) + hemisphereOffset(i, iDir);
            unsigned int i0_prev = i0LineSouth(nV, i-1) + hemisphereOffset(i-1, iDir);

            unsigned int i1 = i0_prev, i2 = i0;
            for (unsigned int j = 0; j < nPatterns; j++) {
                T_octahedron.row(iT++) = RowVector3i(i2, i1, i2+1);

                if (iDir < 3 || j < nPatterns - 1) 
                    T_octahedron.row(iT++) = RowVector3i(i2+1, i1, i1+1);
                else
                    T_octahedron.row(iT++) = RowVector3i(i2+1, i1, i0LineSouth(nV, i-1));

                i1++;
                i2++;
            }
            if (iDir < 3)
                T_octahedron.row(iT++) = RowVector3i(i2, i1, i2+1); 
            else
                T_octahedron.row(iT++) = RowVector3i(i2, i0LineSouth(nV, i-1), i0LineSouth(nV, i));
        }
    }

    MatrixXd V_unitSphere = V_unitOctahedron.rowwise().normalized();
    nOctahedron = nRes;

    // Set octahedron dimensions according to input points
    fittingOctahedron();
    V_octahedron = MatrixXd(V_unitOctahedron.rows(), 3);
    V_octahedron.col(0) = V_unitOctahedron.col(0) * octaDim(0);
    V_octahedron.col(1) = V_unitOctahedron.col(1) * octaDim(1);
    V_octahedron.col(2) = V_unitOctahedron.col(2) * octaDim(2);
    V_octahedron = (V_octahedron * octaRot.inverse()).rowwise() + octaPos.transpose();

    V_ellipse = MatrixXd(V_unitSphere.rows(), 3);
    V_ellipse.col(0) = V_unitSphere.col(0) * octaDim(0);
    V_ellipse.col(1) = V_unitSphere.col(1) * octaDim(1);
    V_ellipse.col(2) = V_unitSphere.col(2) * octaDim(2);
    V_ellipse = (V_unitSphere * octaRot.inverse()).rowwise() + octaPos.transpose();
}

void ApproxSurface::fittingOctahedron() 
{
    // SVD decomposition
    std::vector<Coord> coords; 
    curvesSrc->getCoords(coords, computeAsTape);

    Vector3d c = centroid(coords);
    P.col(0) = c;

    MatrixXd M(3, coords.size());

    int cpt = 0;
    // Store each point position 
    for (Coord &coord : coords) {
        M.col(cpt) = coord.pos - c;
        cpt++;
    }

    // SVD decomposition of the M matrix using BDC method
    // (Bidiagonal Divide & Conquer)
    Eigen::BDCSVD<MatrixXd> svd(M, ComputeThinU);

    Matrix3d U = svd.matrixU();
    Eigen::Vector3d limits = curvesSrc->getFrameLimits(computeAsTape, c, U);

    octaPos = c;
    octaDim = limits;
    octaRot = U;
}

Vector3i ApproxSurface::octaBoundingTriangle(const Vector3d &pos)
{
	unsigned int iT;
	const unsigned int &N = nOctahedron;
	unsigned int nT = 8 * N * N;

	double alpha = fabs(pos(0)) + fabs(pos(1)) + fabs(pos(2));
	double x = fabs(pos(0) / alpha);
	double y = fabs(pos(1) / alpha);
	double z = fabs(pos(2) / alpha);

	unsigned int xInt = floor(x * N); 
	unsigned int yInt = floor(y * N); 
	unsigned int zInt = floor(z * N);

	unsigned int iDeltaT;
	if (pos(2) >= 0.0) {
		iDeltaT = 4 * (N - zInt - 1) * (N - zInt - 1);
	} else {
		iDeltaT = nT / 2 +  4 * (N - zInt - 1) * (N - zInt - 1); 
	}

	unsigned int lineTrianglesNb = 2 * (N - zInt - 1) + 1;
	unsigned int offsetTriangle;
	if ((pos(0)>=0) && (pos(1)>=0)) {
		iT = iDeltaT + 2 * yInt;
		offsetTriangle = 1;
	} else if ((pos(0)<0) && (pos(1)>=0)) {
		iT = iDeltaT + lineTrianglesNb + (lineTrianglesNb - 1 - 2 * yInt);
		offsetTriangle = -1;
	} else if ((pos(0)<0) && (pos(1)<0)) {
		iT = iDeltaT + 2 * lineTrianglesNb + 2 * yInt;
		offsetTriangle = 1;
	} else {
		iT = iDeltaT + 3 * lineTrianglesNb + (lineTrianglesNb - 2 * yInt - 1);
		offsetTriangle = -1;
	}

	if (xInt + yInt + zInt == N) {
		if (zInt == N) {
			iT = 0;
		} else if (yInt == N) {
			iT = iDeltaT + lineTrianglesNb - 1;
		}
	} 
	else if (xInt + yInt + zInt != N - 1) {
		iT += offsetTriangle;
	}

	return Eigen::Vector3i(T_octahedron.row(iT)(0), T_octahedron.row(iT)(1), T_octahedron.row(iT)(2));
}

// TODO: implement point weights consideration
void ApproxSurface::buildFromOctahedronBarycentric() 
{
    unsigned int nPos = curvesSrc->totalSize(computeAsTape);
    Ac = SparseMatrix<double>(nPos, V_octahedron.rows());
    bc = MatrixXd(nPos, 3);

    std::vector<Curve> curves;
    curvesSrc->getCurves(curves, computeAsTape);

    std::vector<Triplet<double>> nonZeroCoords;
    nonZeroCoords.reserve(3 * nPos);

    unsigned int cpt = 0;
    // Loop on each curve position 
    for (Curve &curve : curves) {
        for (Coord &coord : curve.points) {
            RowVector3d pos = coord.pos.transpose();
            RowVector3d posUnitOcta = VectorSpaces::projectToUnitOctahedronNearest(pos, octaPos, octaRot, octaDim);
            Vector3i triangle = octaBoundingTriangle(posUnitOcta);
            RowVector3d posOcta = VectorSpaces::projectToOctahedronNearest(pos, octaPos, octaRot, octaDim);

            RowVector3d v0 = V_octahedron.row(triangle(0)),
                        v1 = V_octahedron.row(triangle(1)),
                        v2 = V_octahedron.row(triangle(2));
            
            Vector2d pos2D = VectorSpaces::trianglePlanarCoords(posOcta, v0, v1, v2);
            Vector3d posBary = VectorSpaces::barycentric2D(pos2D, Vector2d(0,0), Vector2d(1,0), Vector2d(0,1));

            nonZeroCoords.push_back(Triplet<double>(cpt, triangle(0), posBary(0)));
            nonZeroCoords.push_back(Triplet<double>(cpt, triangle(1), posBary(1)));
            nonZeroCoords.push_back(Triplet<double>(cpt, triangle(2), posBary(2)));

            bc.row(cpt) = pos;

            cpt++;
        }
    }

    Ac.setFromTriplets(nonZeroCoords.begin(), nonZeroCoords.end());
}

void ApproxSurface::computeFromOctahedron() 
{
    buildFromOctahedronBarycentric();

    T = T_octahedron;

	// Compute cotangent matrix from ellipse instead of octahedron (smooth result)
    cotangentsMatrix(V_ellipse, T, L);
    sparseVerticalConcat(Ac, L, Acl);
    SparseMatrix<double> Acl_w;
    sparseVerticalConcat(sqrt(1 / smoothing) * Ac, L, Acl_w);

    bcl = MatrixXd(bc.rows() + L.rows(), 3);
    MatrixXd bcl_w(bc.rows() + L.rows(), 3);
    bcl_w << sqrt(1 / smoothing) * bc,
             MatrixXd::Zero(L.rows(), 3);

	sparseLSsolve(Acl_w, bcl_w, Vcl, true);
}

void ApproxSurface::compute(const unsigned int &res, BaseShape baseShape) 
{
    resetMesh();
    // We need at least 3 coordinates to compute a surface
    if (curvesSrc->totalSize(false) < 3) return;

    switch (baseShape)
    {
    case PLANE:
        planarMesh(res);
        computeFromPlane();
        addBorderTriangles();
        break;
    case OCTAHEDRON:
        octahedronMesh(res);
        computeFromOctahedron();
        break;
    default:
        break;
    }
}

void ApproxSurface::resetMesh() 
{
    T = MatrixXi();
    Vcl = MatrixXd();

	// Reset triangle indices map
	rowIndexOfTriangle.clear();
}

void ApproxSurface::computeNormals(MatrixXd &N)
{
    N = MatrixXd::Zero(Vcl.rows(), 3);
    for (unsigned int i = 0; i < T.rows(); i++) {
        RowVector3d v0 {Vcl.row(T.row(i)[0])},
                    v1 {Vcl.row(T.row(i)[1])},
                    v2 {Vcl.row(T.row(i)[2])};

        RowVector3d normal = (v1 - v0).cross(v2 - v0);
        normal.normalize();
        for (unsigned int j = 0; j < 3; j++) {
            unsigned int j0 {j}, j1 {(j + 1) % 3}, j2 {(j + 2) % 3};
            RowVector3d e01 {Vcl.row(T.row(i)[j1]) - Vcl.row(T.row(i)[j0])};
            RowVector3d e02 {Vcl.row(T.row(i)[j2]) - Vcl.row(T.row(i)[j0])};

            // Weight normal contribution to the vertex according to the incident angle
            double angle = acos(e01.dot(e02) / (e01.norm() * e02.norm()));
            N.row(T.row(i)[j]) += angle * normal;
        }
    }

    N.rowwise().normalize();
}
