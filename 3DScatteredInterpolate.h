#ifndef _3D_SCATTERED_INTERPOLATE_H_
#define _3D_SCATTERED_INTERPOLATE_H_

#include <vector>
#include <Eigen/Dense>
#include "libnabo/nabo/nabo.h"

namespace _3DScatteredInterpolate
{
    struct Point
    {
        double x;
        double y;
        double z;
        
        Point(double _x, double _y, double _z) : x(_x), y(_y), z(_z)
        {

        }
    };

    enum InterpolateMethod
    {
        INTERPOLATE_METHOD_LINEAR = 0,
        INTERPOLATE_METHOD_NEAREST
    };

    class ScatteredInterpolate
    {
    public:
        ScatteredInterpolate(std::vector<Point> &_sample_points, std::vector<double> &_values, unsigned int _neighbor_size = 10,
                            InterpolateMethod interpolate_method = INTERPOLATE_METHOD_LINEAR, bool _is_do_extrapolate = false);
        ~ScatteredInterpolate();

        /*
        * values of sample points, user could replace values freely after sample points are recored,
        * without initializing again
        */
        std::vector<double> values;

        double operator () (Point query_point);
        std::vector<double> operator () (std::vector<Point> &query_points);

    private:
        // kdtree
        Nabo::NNSearchF* nns;

        unsigned int neighborSize;

        Eigen::MatrixXf samplePoints;
        bool isExtraPolationAvailable;
        InterpolateMethod interpolateMethodToUse;
        
        double maxDistanceBetweenNeighbors;

        double doLinearInterpolate(Point query_point);
        double findNearestValue(Point query_point);
    };

}
#endif