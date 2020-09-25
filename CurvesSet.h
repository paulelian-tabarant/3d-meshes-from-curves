#pragma once

#define MAXSIZE 256

#include <vector>
#include <Eigen/Dense>
#include <string>

/**
 * @brief 
 * Stores a unique 3D curve coordinate as a triplet of
 * one time value,
 * three position coefficients 
 * & four rotation (quaternion) coefficients
 */
struct Coord {
    Eigen::Vector3d pos;
    Eigen::Vector3d normal;
    double weight;

    Coord(const Eigen::Vector3d &p, const Eigen::Vector3d &n, const double &w = 1.0) :
        pos(p), normal(n), weight(w) {}

    Coord() : pos(Eigen::Vector3d::Zero()), 
              normal(Eigen::Vector3d(0.0, 1.0, 0.0)),
              weight(1.0) {}
};

struct Curve {
    std::vector<Coord> points;
    double weight;
    double width;

    Curve(const double &weight, const double &width) :
        weight(weight), width(width) {}
};

/**
 * @brief 
 * Class handling extraction and management of 3D curves data
 * Inside files, each curve should be separated by a line jump.
 * Each coordinate (time, position, quaternion rotation)
 * must be written on an individual line, and formatted as follows :
 * time X Y Z w x y z
 */
class CurvesSet 
{
    private:
        // Original raw curves data (time, pos, rot)
        std::vector<Curve> curves;
        std::vector<Curve> tapeCurves;
        bool tapeUpdated {false};
		double curvesWidth {0.0f};
    
    private:
        void computeTapeCurves();

    public:
        /**
         * @brief Construct a new Curves Set object
         * from data stored inside the file at path.
         */ // not maintained anymore
        //CurvesSet(std::string path, int decimation);

        CurvesSet(const std::vector<Curve>& curves) : curves(curves) { }

        CurvesSet() {}

        /**
         * @brief Clear curves input data
         */
        void clear() { curves.clear(); }

        void startNewCurve(const double &weight, const double &width);
		void startNewCurve(const double &weight) { startNewCurve(weight, curvesWidth); }
		void addPoint(const Coord &c);
        void addPoint(const Eigen::Vector3d &p, const Eigen::Vector3d &n, const double &w) { addPoint(Coord(p, n, w)); }

        /**
         * @brief Add new points to the current curve,
         * or a newly created curve if necessary
         */
        void addPoints(const std::vector<Coord>& points);

        /**
         * @brief Remove last curve data from the curves set
         */
        void removeLastCurve();

        /**
         * @brief 
         * Assign different weights to curves for surface generation.
         * If specified vector is smaller than nb of curves, other weights
         * will be considered as 1.
         */
        //void applyWeightsToCurves(Eigen::VectorXd weights);

        /**
         * @brief Set last curve width inside the input set 
         */
        void setCurvesWidth(const double &width) { curvesWidth = width; }

        /**
         * @brief Get the total number of 3D points
         * inside the curves set
         * @param tape if the size must consider "tape" representation
         * of the curves set to retrieve the size (lateral points)
         */
        unsigned int totalSize(const bool &tape = {false});

        /**
         * @brief Get the number of curves in the set 
         */
        unsigned int size() { return curves.size(); }

        /**
         * @brief 
         * Get bounds for projection of all 3D points on the
         * plan coordinates system
         * @param plane The plane space (origin, X, Y, Z in columns)
         * @return Eigen::Vector4d minX, minY, maxX, maxY on this plane
         */
        Eigen::Vector4d getProjectionBounds(const bool &tape, const Eigen::Matrix<double, 3, 4> &plane);

        /**
         * @brief Get x, y and z absolute limits (from origin), given the
         * current positions of curves.
         * @return Eigen::Vector3d dx, dx and dz from the origin
        */
        Eigen::Vector3d getFrameLimits(const bool &tape, const Eigen::Vector3d &O, const Eigen::Matrix3d &frame);

        /**
         * @brief Get curves coordinates stored in independent vectors
         */
        void getCoords(std::vector<Coord> &coords, const bool &tape);

        void getCurves(std::vector<Curve> &curves, const bool &tape);
};