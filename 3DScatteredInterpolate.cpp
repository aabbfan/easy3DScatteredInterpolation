#include "3DScatteredInterpolate.h"
#include <cstdio>
#include <cfloat>

using namespace _3DScatteredInterpolate;

ScatteredInterpolate::ScatteredInterpolate(std::vector<Point> &_sample_points, std::vector<double> &_values,
                                           unsigned int _neighbor_size, InterpolateMethod interpolate_method,
                                           bool _is_do_extrapolate)
                                           : values(_values), interpolateMethodToUse(interpolate_method), isExtraPolationAvailable(_is_do_extrapolate),
                                             neighborSize(_neighbor_size)
{
    if (_sample_points.size() < _neighbor_size)
    {
        fprintf(stderr, "Error! too few input sample points.\n");
        exit(0);
    }

    /* 
     * initialize sample points
     * libabo requires points matrix with (3 * N) dimension
     */
    this->samplePoints = Eigen::MatrixXf(3, _sample_points.size());

    for (int i = 0; i < _sample_points.size(); i++)
    {
        this->samplePoints(0, i) = _sample_points[i].x;
        this->samplePoints(1, i) = _sample_points[i].y;
        this->samplePoints(2, i) = _sample_points[i].z;
    }

    // initialize kdtree
    this->nns = Nabo::NNSearchF::createKDTreeLinearHeap(this->samplePoints);

    // initialize max distance parameter
    this->maxDistanceBetweenNeighbors = DBL_MIN;

    #pragma omp parallel for
    for (int i = 0; i < this->samplePoints.cols(); i++)
    {
        // find the nearest neighbor of this point
        Eigen::VectorXi indices(1);
        Eigen::VectorXf dists(1);

        Eigen::Vector3f point(this->samplePoints(0, i), this->samplePoints(1, i), this->samplePoints(2, i));
        this->nns->knn(point, indices, dists, 1, 0); // only one point

        if (dists[0] > this->maxDistanceBetweenNeighbors)
        {
            #pragma omp critical
            this->maxDistanceBetweenNeighbors = dists[0];
        }
    }
}

ScatteredInterpolate::~ScatteredInterpolate()
{
    delete nns;
}

double ScatteredInterpolate::doLinearInterpolate(Point query_point)
{
    Eigen::Vector3f point(query_point.x, query_point.y, query_point.z);

    Eigen::VectorXi indices(this->neighborSize);
    Eigen::VectorXf dists(this->neighborSize);
    this->nns->knn(point, indices, dists, this->neighborSize, 0, Nabo::NNSearchF::ALLOW_SELF_MATCH);

    // inverse distance weighting method
    double dist_sum = 0;
    double result = 0;
    for (int i = 0; i < indices.size(); i++)
    {
        dist_sum += 1 / dists[i] / dists[i];
        result += values[indices[i]] / dists[i] / dists[i];
    }
    result /= dist_sum;

    return result;
}

double ScatteredInterpolate::findNearestValue(Point query_point)
{
    // find the nearest one
    Eigen::Vector3f point(query_point.x, query_point.y, query_point.z);

    Eigen::VectorXi indices(1);
    Eigen::VectorXf dists(1);

    // nearest
    this->nns->knn(point, indices, dists, 1, 0, Nabo::NNSearchF::ALLOW_SELF_MATCH);

    return this->values[indices[0]];
}

double ScatteredInterpolate::operator () (Point query_point)
{
    // size of values and samplePoints should be the same
    if (this->values.size() != this->samplePoints.cols())
    {
        fprintf(stderr, "Error! the size of sample points and values should be the same.\n");
        return NAN;
    }

    // check size of sample points
    if (this->samplePoints.cols() == 0) return NAN;

    bool is_inside_domain = true;    
    if (this->samplePoints.cols() < this->neighborSize)
        is_inside_domain = false;

    // find the nearest one
    Eigen::Vector3f point(query_point.x, query_point.y, query_point.z);

    Eigen::VectorXi indices(1);
    Eigen::VectorXf dists(1);

    // nearest
    this->nns->knn(point, indices, dists, 1, 0, Nabo::NNSearchF::ALLOW_SELF_MATCH);

    // distance should below global maximum value
    if (dists[0] > this->maxDistanceBetweenNeighbors)
        is_inside_domain = false;

    // do interpolate
    if (is_inside_domain)
    {
        if (this->interpolateMethodToUse == INTERPOLATE_METHOD_LINEAR)
            return doLinearInterpolate(query_point);
        else
            return findNearestValue(query_point);
    }

    // do extrapolate
    if (isExtraPolationAvailable)
        return findNearestValue(query_point);

    return NAN;
}

std::vector<double> ScatteredInterpolate::operator () (std::vector<Point> &query_points)
{
    std::vector<double> result = std::vector<double> (query_points.size());

    #pragma omp parallel for
    for (int i = 0; i < query_points.size(); i++)
    {
        result[i] = (*this)(query_points[i]);
    }

    return result;
}
