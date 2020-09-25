#include "CurvesSet.h"
#include "VectorSpaces.h"

#include <fstream>
#include <cmath>

#if 0
// Caution : not maintained anymore
CurvesSet::CurvesSet(string path, int decimation) : curves(vector<vector<Coord>>()) 
{ 
    // Create a stream & a string to read curves coordinates line by line
    ifstream ifs(path);
    string coordStr;
    int curvesCpt = 0, pointsCpt = 0;

    // Scan the input file until the end
    while (getline(ifs, coordStr)) {
        istringstream iss(coordStr);

        if (coordStr.length() <= 1) {
            // Curve width value is stored as 1st line of a curve block
            double width;
            iss >> width;
            curvesWidth.push_back(width);
            // Line jump : create a new curve
            Curve newCurve(1.0, width);
            curves.push_back(newCurve);
            curvesCpt++;
            continue;
        }

        // If we want to under-sample the curves
        pointsCpt++;
        if (pointsCpt % decimation != 0) continue;

        // Add new coordinate to current curve
        Eigen::Vector3d p;
        Eigen::Quaterniond r;
        double w;

        // time, position & rotation reading
        iss >> p(0) >> p(1) >> p(2)
            >> r.w() >> r.x() >> r.y() >> r.z()
            >> w;

        Coord c(p, r, w);
        curves[curvesCpt].points.push_back(c);
    }

    ifs.close();
}
#endif

void CurvesSet::startNewCurve(const double &weight, const double &width)
{
    if (curves.size() > 0 && curves.back().points.size() == 0)
        return;
    
    curves.push_back(Curve(weight, width));
}

void CurvesSet::addPoint(const Coord &c)
{
    if (curves.size() == 0)
        return;
    
    curves.back().points.push_back(c);
    tapeUpdated = false;
}

void CurvesSet::addPoints(const std::vector<Coord>& points)
{
    if (curves.size() == 0)
        return;

    for (Coord point : points) {
        curves.back().points.push_back(point);
    }

    tapeUpdated = false;
}

void CurvesSet::computeTapeCurves()
{
    tapeCurves.clear();

    for (Curve &curve : curves) {
        Curve tapeCurve(curve.weight, curve.width);
        for (size_t ip = 0; ip < curve.points.size(); ip++) {
            // Last point of the curve : assume having same segment value as previous point
            Coord c = curve.points[ip];
            Eigen::Vector3d p0 = c.pos;
            Eigen::Vector3d n = c.normal.normalized();
            // Segment computation
            Eigen::Vector3d vd;
            if (ip < curve.points.size() - 1) {
                Eigen::Vector3d p1 = curve.points[ip + 1].pos;
                vd = n.cross(p1 - p0).normalized();
            }
            else {
                Eigen::Vector3d pPrev = curve.points[ip - 1].pos;
                vd = n.cross(p0 - pPrev).normalized();
            }

            // Create two-sided ribbon from each coord
            Coord cl = c, cr = c;
            cl.pos = p0 - curve.width/2 * vd;
            cr.pos = p0 + curve.width/2 * vd;

            tapeCurve.points.push_back(cl);
            tapeCurve.points.push_back(c);
            tapeCurve.points.push_back(cr);
        }
        tapeCurves.push_back(tapeCurve);
    }

    tapeUpdated = true;
}

void CurvesSet::removeLastCurve() 
{
    if (curves.size() == 0) return;

    curves.pop_back();

    tapeUpdated = false;
}

unsigned int CurvesSet::totalSize(const bool &tape)
{
    unsigned int count {0};
    for (Curve& curve : curves)
        count += curve.points.size();

    return tape ? 3 * count : count;
}

Eigen::Vector4d CurvesSet::getProjectionBounds(const bool &tape, const Eigen::Matrix<double, 3, 4> &plane)
{
    // As we often compute from the centroid, not necessary to use sup and inf bounds
    double minX = 0, minY = 0, maxX = 0, maxY = 0; 
    if (tape && !tapeUpdated)
        computeTapeCurves();

    std::vector<Curve> &curCurves = tape ? tapeCurves : curves;
    
    for (Curve &curve : curCurves) {
        // Compute X and Y distance of each point from the origin
        for (Coord &coord : curve.points) {
            Eigen::Vector2d xy = VectorSpaces::orthoProjection(coord.pos, plane);
            double x = xy(0),
                   y = xy(1);
            minX = x < minX ? x : minX;
            minY = y < minY ? y : minY;
            maxX = x > maxX ? x : maxX;
            maxY = y > maxY ? y : maxY;
        }
    }
    return Eigen::Vector4d(minX, minY, maxX, maxY);
}

Eigen::Vector3d CurvesSet::getFrameLimits(const bool &tape, const Eigen::Vector3d &O, const Eigen::Matrix3d &frame)
{
    Eigen::Vector3d limits;

    if (tape && !tapeUpdated)
        computeTapeCurves();

    std::vector<Curve> &curCurves = tape ? tapeCurves : curves;

    for (Curve &curve : curCurves) {
        for (Coord &coord : curve.points) {
            Eigen::Vector3d framePos = (coord.pos - O).transpose() * frame;
            double x = fabs(framePos(0)),
                   y = fabs(framePos(1)),
                   z = fabs(framePos(2));
            limits(0) = x > limits(0) ? x : limits(0);
            limits(1) = y > limits(1) ? y : limits(1);
            limits(2) = z > limits(2) ? z : limits(2);
        }
    }
    return limits;
}

#if 0
void CurvesSet::applyWeightsToCurves(Eigen::VectorXd weights) 
{
    // If weights vector is greater than number of curves,
    // we trim it. If smaller, we fill it with ones until
    // we reach the curves number size.
	size_t size0 = weights.size();
    Eigen::VectorXd weights0 = weights;
    weights = Eigen::VectorXd::Constant(size(), 1.0);
    if (size0 < curves.size())
        weights.head(size0) = weights0;
    else
        weights = weights0.head(curves.size());

    // Normalize weights in order not to influence global
    // surface computation (<c, alpha*w> = totalSize)
    Eigen::VectorXd c(curves.size());
    for (size_t i = 0; i < curves.size(); ++i) {
        c[i] = (curves[i]).points.size();
    }
    double alpha = totalSize() / c.dot(weights);
    for (size_t i = 0; i < curves.size(); ++i) {
        (curves[i]).weight = alpha * weights[i];
    }
}
#endif

void CurvesSet::getCoords(std::vector<Coord> &coords, const bool &tape) 
{
    if (tape && !tapeUpdated)
        computeTapeCurves();

    std::vector<Curve>& curCurves = tape ? tapeCurves : curves;

    for (Curve &curve : curCurves) {
        for (Coord &c : curve.points) {
            coords.push_back(c);
        }
    }
}

void CurvesSet::getCurves(std::vector<Curve> &curCurves, const bool &tape)
{
    if (tape && !tapeUpdated)
        computeTapeCurves();

    curCurves = tape ? tapeCurves : curves;
}

