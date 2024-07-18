#include "3DScatteredInterpolate.h"
#include <iostream>
#include <cmath>

double sampleFunction(double x, double y, double z)
{
    return std::sin(x) + 2 * y + std::cos(z);
}

int main()
{
    std::vector<_3DScatteredInterpolate::Point> sample_points;
    std::vector<double> sample_values;
    
    for (double x = 0; x < 10; x = x + 0.1)
    {
        for (double y = 0; y < 10; y = y + 0.1)
        {
            for (double z = 0; z < 10; z = z + 0.1)
            {
                _3DScatteredInterpolate::Point tmp_point(x, y, z);
                double tmp_value = sampleFunction(x, y, z);

                sample_points.push_back(tmp_point);
                sample_values.push_back(tmp_value);
            }
        }
    }

    _3DScatteredInterpolate::ScatteredInterpolate interpolate(sample_points, sample_values);
    std::cout << "sample function: sin(x) + 2*y + cos(z)" << std::endl;
    std::cout << "interpolating result for point(1.13, 2.51, 3.21) :" << std::endl;
    
    _3DScatteredInterpolate::Point query_point(1.13, 2.51, 3.21);
    std::cout << interpolate(query_point) << std::endl;

    std::cout << "sin(1.13) + 2*(2.51) + cos(3.21):" << std::endl;
    std::cout << sampleFunction(1.13, 2.51, 3.21) << std::endl;

    return 0;
}