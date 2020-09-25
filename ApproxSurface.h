#pragma once

#include "CurvesSet.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <map>
#include <memory>

#define TEST_FITTING_PLANE ;
#define LSCG ;

class CurvesSet;

/**
 * @brief 
 * Generate and handle operation on approximated 3D surfaces 
 * based on one instance of the CurvesSet class
 */
class ApproxSurface 
{
    public:
        enum BaseShape {
            PLANE = 0,
            OCTAHEDRON = 1
        };

    private:
        std::shared_ptr<CurvesSet> curvesSrc;

        /**
         * @brief Tell if current input curves should be considered
         * as tape or just raw curves (width == 0)
         */
        const bool computeAsTape;
		const double smoothing;

        /**
         * @brief The curves set fitting plane, 
         *  if already computed via planarMesh(M, N)
         * as well as its 2d boundaries on local coordinates system
         */
        Eigen::Matrix<double, 3, 4> P; 
        double xMin, yMin, xMax, yMax;
        bool constantDirs {false};
        // When fitting plane is set to constant, the former must be
        // computed once at least
        bool computeFittingPlane {false};

        std::vector<unsigned int> borderV;

        // Triangles : we draw a tiling of hexagons based on the previous vertices
        unsigned int t1[2], t2[2], t3[2], t4[2], t5[2], t6[2];
        double theta1, theta2;

        /**
         * @brief Store first hexagon center index in the list
         * of drawn vertices, & step value between this center 
         * and the center of the same hexagon on the next row
         */
        unsigned int h0, hstep;

		Eigen::MatrixXd V_plane;

        /**
         * @brief Approximated surface triangles,
         * valid as well for the planar mesh as for
         * the actual approximated surface
         */
        Eigen::MatrixXi T;

        /**
         * @brief If any triangulation computed,
         * current mesh resolution
         * on width (M) and length (N)
         */
        unsigned int mPlane, nPlane;

        /**
         * @brief Store the laplacian cotangent matrix from
         * the initial planar mesh {pos2D, triangles}.
         */
        Eigen::SparseMatrix<double> L;

        /**
         * @brief Matrix storing coefficients for solving
         * the barycentric coordinates problem, & adding 
         * laplacian constraints in a second step
         */
        Eigen::SparseMatrix<double> Ac, Acl;

        /**
         * @brief Store the right-hand member for solving
         * the barycentric coordinates problem, & adding
         * laplacian constraints in a second step
         */
        Eigen::MatrixXd bc, bcl;

        /**
         * @brief Positions of vertices of the real 
         * approximated surface according to resp. 
         * barycentric and laplacian + barycentric constraints.
         */
        Eigen::MatrixXd Vcl;
        /**
         * @brief Map associating each index triplet with the index
         * of the corresponding row inside T
         */
        std::map<std::string, int> rowIndexOfTriangle;

        // To be documented
        unsigned int nOctahedron;
        Eigen::MatrixXd V_octahedron, V_ellipse;
        Eigen::MatrixXi T_octahedron;


		// Octahedron dimensions
        Eigen::Vector3d octaPos, octaDim;
        Eigen::Matrix3d octaRot;

private:
        /**
         * @brief Compute the centroid value from a set of 3D points
         */
        Eigen::Vector3d centroid(const std::vector<Coord> &coords);

        /**
         * @brief Compute the best fitting plane from the current
         * curves set using a SVD decomposition
         */
        void fittingPlane();

        /**
         * @brief Compute planar triangulation with a given
         * N edges resolution per row, adapting length to have
         * as-regular-as possible hexagons
         */
        void planarMesh(const unsigned int &N);

        /**
         * @brief Compute a triangulation of the fitting plane
         * @param M width resolution in triangles
         * @param N length resolution in triangles
         */
        void planarMesh(const unsigned int &M, const unsigned int &N);

        /**
         * @brief Convert a triplet of integers representing a 
         * triangle inside T topology to a suitable string for
         * unique mapping of this triplet
         */
        std::string triangleToMapStr(unsigned int i0, unsigned int i1, unsigned int i2);

        /**
         * @brief Get the index of the 3 vertices defining
         * the bounding triangle of a couple (x,y) planar 
         * coordinates 
         */
        unsigned int planeBoundingTriangle(const double &x, const double &y);

        /**
         * @brief Concatenate two sparse matrices vertically 
         * into last matrix argument
         */
        void sparseVerticalConcat(
            const Eigen::SparseMatrix<double> &M1,
            const Eigen::SparseMatrix<double> &M2,
            Eigen::SparseMatrix<double>& M);

        /**
         * @brief Solve a linear system composed of a SparseMatrix A
         * and a right member b which is a dense matrix, in the least-
         * squares sense. Last parameter should indicate whether we 
         * need to factorize the A matrix (true) or if last computation
         * is still relevant (false).
         */
        void sparseLSsolve(const Eigen::SparseMatrix<double> &A,
            const Eigen::MatrixXd &b,
            Eigen::MatrixXd &X,
            const bool &factorize);

        /**
         * @brief Compute coordinates of a triangle defined by vertices index, in the fitting plane 
         * 2d coordinates system 
         */
        void trianglePlanarCoords(const Eigen::Vector3i &abc, Eigen::Vector2d& a2d, Eigen::Vector2d& b2d, Eigen::Vector2d& c2d);

        void cotangentsMatrix(const Eigen::MatrixXd &vertices, const Eigen::MatrixXi &triangles, Eigen::SparseMatrix<double> &cotangents);

        /**
         * @brief Compute each point barycentric coordinates
         * inside the triangular planar mesh, store them into
         * a sparse matrix for further surface approximation 
         */
        void buildFromPlaneBarycentric();

        void computeFromPlane();
        void computeFromOctahedron();

        /**
         * @brief Add triangles to the planar mesh to obtain a perfectly
         * rectangular plane. 
         */
        void addBorderTriangles();

        /**
         * @brief Build the map associating a triplet of vertices to
         * the index of corresponding triangle row index inside the 
         * mesh matrix (T_tape or T_raw)
         */
        void initTriangleIndicesData();

		// Octahedron mesh indices facilities
		unsigned int i0LineNorth(const unsigned int &nLine);
		unsigned int i0LineSouth(const unsigned int &nV, const unsigned int &nLine);
		unsigned int hemisphereOffset(const unsigned int &nLine, const unsigned int &hIndex);

		/**
		 * @brief Build an octahedron-shaped (normalized) mesh from a resolution
		 * corresponding to the number of triangle rows on each face
		 */
        void octahedronMesh(const unsigned int &nRes);

		/**
		 * @brief Compute a bounding octahedron from input positions (centroid,
		 * orientation and dimensions) stored inside private fields 
		 */
        void fittingOctahedron();

		/**
		 * @brief Determine the index of the bounding triangle of an input 3D position
		 * within the current octahedron (which should be already computed).
		 * Input position should be given from a "unit octahedron" coordinate system.
		 */
        Eigen::Vector3i octaBoundingTriangle(const Eigen::Vector3d &pos);

		/**
		 * @brief Build the coefficients matrix necessary for the computation of an
		 * output geometry, from a bounding octahedron (which should have been computed before).
		 * Each non-zero coordinate on a row, for example on column i, corresponds to a barycentric
		 * coordinate associated to one vertex of its bounding triangle.
		 */
        void buildFromOctahedronBarycentric();

		/**
		 * @brief Remove output geometry data, stored inside T_raw and Vcl
		 */
        void resetMesh();

    public:
        /**
         * @brief Construct a new Approx Surface object
         * @param curvesSrc
         * A CurvesSet instance which the approx surface should be
         * based on
         */
        ApproxSurface(const std::shared_ptr<CurvesSet> &curvesSrc, const bool &tape, const double &smoothing)
			: curvesSrc(curvesSrc), computeAsTape(tape), smoothing(smoothing) { }

        /**
         * @brief Retrieve the last computed fitting plane
         * @return Matrix<double, 3, 4>& 
         */
        void getFittingPlane(Eigen::Matrix<double, 3, 4>& plane) { plane = this->P; }
        Eigen::Vector4d getPlaneLimits() { return Eigen::Vector4d(xMin, yMin, xMax, yMax); }

        /**
         * @brief Get triangles topology of approximated
         * surface. Valid for both planar & 3D triangulation
         * @return The set of triangles identified by their
         * vertices indices inside the positions matrix
         */
        void getTriangles(Eigen::MatrixXi& T) { T = this->T; };

        /**
         * @brief Compute approximative surface based on the current curves set.
         * @param res width resolution in terms of number of edges
         * @param w laplacian smooting (> 1 if you want smoothing to have more influence
         * than approximating positions)
         */
        void compute(const unsigned int &res, BaseShape baseShape); 

        /**
         * @brief Get the surface vertices positions 
         */
        void getVertices(Eigen::MatrixXd& V) { V = this->Vcl; }

        /**
         * @brief Tell if the fitting plane directions should be recomputed
         * each time compute() is called, or left unmodified
         */
        void setConstantFittingPlane(const bool &constant);

        void computeNormals(Eigen::MatrixXd &N);
};